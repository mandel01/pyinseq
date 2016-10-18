#!/usr/bin/env python

'''Main script for running the pyinseq package.'''
import argparse
import csv
import glob
import logging
import os
import numpy as np
import pandas as pd
import regex as re
import screed
import sys
import yaml
from shutil import copyfile
from collections import OrderedDict
from .demultiplex import demultiplex_fastq, write_reads
from .gbkconvert import gbk2fna, gbk2ftt
from .mapReads import bowtie_build, bowtie_map, parse_bowtie
from .processMapping import mapGenes, buildGeneTable
from .utils import convert_to_filename, createExperimentDirectories

logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s %(message)s')
logger = logging.getLogger(__name__)


def parseArgs(args):
    '''Parse command line arguments.'''
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',
                        help='input Illumina reads file or folder',
                        required=True)
    parser.add_argument('-s', '--samples',
                        help='sample list with barcodes. \
                        If not provided then entire folder provided for --input is analyzed',
                        required=False)
    parser.add_argument('-e', '--experiment',
                        help='experiment name (no spaces or special characters)',
                        required=True)
    parser.add_argument('-g', '--genome',
                        help='genome in GenBank format (one concatenated file for multiple contigs/chromosomes)',
                        required=True)
    parser.add_argument('-d', '--disruption',
                        help='fraction of gene disrupted (0.0 - 1.0)',
                        default=1.0)
    parser.add_argument('--nobarcodes',
                        help='barcodes have already been removed from the samples; \
                        -i should list the directory with filenames (.fastq.gz) \
                        corresponding to the sample names',
                        action='store_true',
                        default=False)
    parser.add_argument('--demultiplex',
                        help='demultiplex initial file into separate files by barcode',
                        action='store_true',
                        default=False)
    parser.add_argument('--compress',
                        help='compress (gzip) demultiplexed samples',
                        action='store_true',
                        default=False)
    parser.add_argument('--keepall',
                        help='keep all intermediate files generated',
                        action='store_true',
                        default=False)
    return parser.parse_args(args)


class cd:
    '''Context manager to change to the specified directory then back.'''
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


class Settings():
    '''Instantiate to set up settings for the experiment'''
    def __init__(self, experiment_name):
        # standard
        self.experiment = convert_to_filename(experiment_name)
        self.path = 'results/{}/'.format(self.experiment)
        self.genome_path = self.path + 'genome_lookup/'
        self.genome_index_path = self.path + 'genome_lookup/genome'
        self.raw_path = self.path + 'raw_data/'
        self.samples_yaml = self.path + 'samples.yml'
        self.summary_yaml = self.path + 'summary.yml'
        # may be modified
        self.keepall = False
        self.barcode_length = 4

    # Set up directories?


def set_paths(experiment_name):
    experiment = convert_to_filename(experiment_name)
    samples_yaml = 'results/{}/samples.yml'.format(experiment)
    summary_yaml = 'results/{}/summary.yml'.format(experiment)
    path = {'experiment': experiment,
            'samples_yaml': samples_yaml,
            'summary_yaml': summary_yaml}
    return path


def set_disruption(d):
    '''Check that gene disrution is 0.0 to 1.0; otherwise set to 1.0'''
    if d < 0.0 or d > 1.0:
        logger.error('Disruption value provided ({0}) is not in range 0.0 to 1.0; proceeding with default value of 1.0'.format(d))
        d = 1.0
    return d


def tab_delimited_samples_to_dict(sample_file):
    '''Read sample names, barcodes from tab-delimited into an OrderedDict.'''
    samplesDict = OrderedDict()
    with open(sample_file, 'r', newline='') as csvfile:
        for line in csv.reader(csvfile, delimiter='\t'):
            if not line[0].startswith('#'):  # ignore comment lines in original file
                # sample > filename-acceptable string
                # barcode > uppercase
                sample = convert_to_filename(line[0])
                barcode = line[1].upper()
                if sample not in samplesDict and barcode not in samplesDict.values():
                        samplesDict[sample] = {'barcode': barcode}
                else:
                    raise IOError('Error: duplicate sample {0} barcode {1}'.format(sample, barcode))
    return samplesDict


def yaml_samples_to_dict(sample_file):
    '''Read sample names, barcodes from yaml into an OrderedDict.'''
    with open(sample_file, 'r') as f:
        samplesDict = yaml.load(sample_file)
    return samplesDict


def directory_of_samples_to_dict(directory):
    '''Read sample names from a directory of .gz files into an OrderedDict.'''
    samplesDict = OrderedDict()
    for gzfile in list_files(directory):
        # TODO(convert internal periods to underscore? use regex?)
        # extract file name before any periods
        f = (os.path.splitext(os.path.basename(gzfile))[0].split('.')[0])
        samplesDict[f] = {}
    return samplesDict


def list_files(folder, ext='gz'):
    '''Return list of .gz files from the specified folder'''
    with cd(folder):
        return [f for f in glob.glob('*.{}'.format(ext))]


def extract_chromosome_sequence(fastq_file):
    '''
    Return a dictionary of all insertions in a file

    {barcode1:
        {chrom_seq1: count,
         chrom_seq2: count},
     barcode2:
        {chrom_seq1: count,
         chrom_seq2: count}}
    '''
    insertionDict = {}
    pattern = re.compile('''
    ^                               # beginning of string
    ([ACGT]{4})                     # group(1) = barcode, any 4-bp of mixed ACGT
    ([NACGT][ACGT]{13,14}(?:TA))    # group(2) = 16-17 bp of chromosomal sequence
                                    # first bp can be N
                                    # last two must be TA for transposon
    ACAGGTTG                        # flanking transposon sequence
    ''', re.VERBOSE)
    logger.debug('fastq_file: {0}'.format(fastq_file))
    with screed.open(fastq_file) as seqfile:
        for read in seqfile:
            m = re.search(pattern, read.sequence)
            try:
                barcode, chrom_seq = m.group(1), m.group(2)
                try:
                    insertionDict.setdefault(barcode, {chrom_seq: 0})[chrom_seq] += 1
                # insertionDict[barcode] exists but not {chrom_seq: ...} nested dict
                except(KeyError):
                    insertionDict[barcode][chrom_seq] = 1
            # no barcode/chrom_seq in the read
            except:
                pass
    logger.debug('insertionDict: {0}'.format(insertionDict))
    return insertionDict


def map_raw_reads(settings, samplesDict, insertionDict):
    '''
    Map sequences with bowtie.

    Iterate through each barcode in insertionDict and map reads using bowtie.
    Map only reads for barcodes in the samples dictionary (i.e., barcodes that are
    in the experiment)

    Input dictionary is a nested dictionary defined in extract_chromosome_sequence():

    {barcode1:
        {chrom_seq1: count,
         chrom_seq2: count},
     barcode2:
        {chrom_seq1: count,
         chrom_seq2: count}}
    '''
    for sample in samplesDict:
        barcode = samplesDict[sample]['barcode']
        reads = ','.join(insertionDict[barcode].keys())
        bowtie_output_file = settings.path + sample + '_bowtie.txt'
        logger.info('bowtie mapping for sample: {0}: {1}'.format(sample, samplesDict[sample]))
        # logger.debug('reads for bowtie mapping: {0}'.format(reads))
        # logger.debug('bowtie_output_file: {0}'.format(bowtie_output_file))
        bowtie_map(settings.genome_index_path, reads, bowtie_output_file)


def yamldump(d, f):
    '''Write dictionary d as yaml file f.yml'''
    with open(f + '.yml', 'w') as fo:
        fo.write(yaml.dump(d, default_flow_style=False))


def write_individual_insertion_count_dictionaries(insertionDict, samplesDict, settings):
    '''Write all insertion counts to individual yaml files'''
    for sample in samplesDict:
        barcode = samplesDict[sample]['barcode']
        yamldump(insertionDict[barcode], settings.path + sample + '_insertion_counts')


def read_individual_insertion_count_dictionary(sample, settings):
    '''Return pandas dataframe of insertion counts from one yaml file'''
    insertion_count_file = settings.path + sample + '_insertion_counts.yml'
    with open(insertion_count_file, 'r') as f:
        insertion_counts = yaml.load(f)
    df_insertion_counts = pd.DataFrame.from_records(insertion_counts, index=['counts']).transpose()
    df_insertion_counts['sample'] = sample
    return df_insertion_counts


def write_sites_files(dataframe_object, output_file):
    df = dataframe_object.groupby(['contig', 'insertion_nucleotide', 'sample',
                                   'orientation']).sum().unstack(level=['sample', 'orientation'])
    # join the first ('counts') and third ('+' or '-') column header values
    df.columns = [''.join(col[::2]) for col in df.columns.values]
    df.to_csv(output_file, sep='\t')


def read_sites_file(sample, settings):
    '''Return pandas dataframe of insertion counts'''
    df = pd.read_csv(settings.path + sample + '_sites.txt', sep='\t')
    return df


def reverse_complement(sequence):
    '''Return the DNA Reverse Complement'''
    complement = ''.maketrans('GATCRYSWKMBVDHN', 'CTAGYRSWMKVBHDN')
    return sequence.translate(complement)[::-1]


def original_sequence(row):
    '''Return the original sequence from the bowtie sequence (pandas)'''
    if row['orientation'] == '-':
        return reverse_complement(row['sequence'])
    return row['sequence']


def insertion_nucleotide(orientation, bowtie_nucleotide, read_length):
    '''Given bowtie mapping data provide the transposon insertion nucleotide'''
    return bowtie_nucleotide + read_length - 1 if orientation == '+' else bowtie_nucleotide + 1


def normalize_counts_per_million(df_sites):
    '''Normalize plus, minus, and total sites all to total sites and multiply by a million'''
    cpm = lambda x: 1E6 * x / df_sites['total'].sum()
    df_sites[['counts+', 'counts-', 'total']] = df_sites[['plus', 'minus', 'total']].applymap(cpm)
    return df_sites


def process_bowtie_results(settings, samplesDict, insertionDict):
    '''Read each bowtie result file into a dataframe, align with read data

       For each bowtie results file (orientation, replication, nucleotide, sequence):
       - Read into a pandas dataframe
       - Align to the insertion data (number of reads)
       - Write to a csv/txt file
       - Concatenate the results from all of the samples

    '''
    heads = ['orientation', 'contig', 'bowtie_nucleotide', 'sequence', 'mismatch']

    # dataframe of mapped reads of all barcodes
    df_insertions = pd.DataFrame.from_dict(insertionDict)
    logger.debug('df_insertions\n{0}'.format(df_insertions))

    # results for all of the samples
    df_all_samples = ''
    summary_sites_table = 'results/{experiment}/summary_sites_table.txt'.format(
        experiment=settings.experiment)

    for sample in samplesDict:
        # ### MAP TO SITES ### #
        bowtie_results_file = 'results/{experiment}/{sample}_bowtie.txt'.format(
            experiment=settings.experiment,
            sample=sample)

        # dataframe of bowtie results
        df_bt = pd.read_csv(bowtie_results_file, delimiter='\t', names=heads, usecols=[0, 1, 2, 3])
        # Set DNA sequence as the same as the original submission (original or reverse complement)
        df_bt['original_sequence'] = df_bt.apply(lambda row: original_sequence(row), axis=1)
        df_bt['insertion_nucleotide'] = \
            df_bt.apply(lambda row: insertion_nucleotide(row['orientation'],
                                                         row['bowtie_nucleotide'],
                                                         len(row['sequence'])), axis=1)
        # index by the DNA sequence that will be found in the results file
        df_bowtie_results_indexed = df_bt.set_index(['original_sequence'])
        df_insertion_counts = read_individual_insertion_count_dictionary(sample, settings)
        # load bowtie results for the sample
        df_sample = pd.concat([df_insertion_counts[['sample', 'counts']],
                              df_bowtie_results_indexed[['contig', 'orientation', 'insertion_nucleotide']]],
                             axis=1, join='inner')

        sites_summary_file = 'results/{experiment}/{sample}_sites.txt'.format(
            experiment=settings.experiment,
            sample=sample)
        # Write the sites file as tab-delimited csv
        write_sites_files(df_sample, sites_summary_file)
        # Aggregate sites data
        if not isinstance(df_all_samples, pd.DataFrame):
            df_all_samples = df_sample.copy()
        else:
            df_all_samples = df_all_samples.merge(df_sample, left_index=True, right_index=True, how='left')
    # Write all samples
    df_all_samples.to_csv(summary_sites_table, sep='\t')


def process_gene_counts(settings, samplesDict, disruption):
    '''Loop through all samples and map to genes'''
    # Load genome features (.ftt) as dataframe
    df_ftt = pd.read_csv(settings.genome_index_path + '.ftt', sep='\t')
    df_all_samples = ''
    summary_genes_table = 'results/{experiment}/summary_gene_table.txt'.format(
        experiment=settings.experiment)
    # Loop through sample files
    for sample in samplesDict:
        genes_summary_file = 'results/{experiment}/{sample}_genes.txt'.format(
            experiment=settings.experiment,
            sample=sample)
        df_sites = read_sites_file(sample, settings)
        df_sample = map_counted_insertions_to_genes(df_sites, df_ftt, sample, genes_summary_file, disruption)
        # Aggregate genes data
        if not isinstance(df_all_samples, pd.DataFrame):
            df_all_samples = df_sample.copy()
        else:
            df_all_samples = df_all_samples.merge(df_sample, left_index=True, right_index=True, how='left')
    # Write all samples
    df_all_samples.to_csv(summary_genes_table, sep='\t')


def map_counted_insertions_to_genes(df_sites, df_ftt, sample, genes_summary_file, disruption):
    '''Map insertions in df_sample to the genes in df_genome_ftt
       Insertions in multiple genes get counted in both.
    '''
    df_sites['total'] = df_sites['counts+'] + df_sites['counts-']
    df_sites.columns = ['Locus', 'nt', 'plus', 'minus', 'total']
    df_merged = df_ftt.merge(normalize_counts_per_million(df_sites), on=['Locus'], how='left')
    df_merged['count'] = df_merged.apply(lambda x: x['total'] if insertion_falls_in_gene(
        x['Location_Start'], x['Location_End'], x['Strand'], x['nt'], disruption) else 0, axis=1)
    df_aggregate = df_merged[list(df_ftt.columns) + ['count']].groupby(list(df_ftt.columns)).sum()
    df_aggregate_resorted = df_aggregate.sort_index(level=list(df_ftt.columns))
    df_aggregate_resorted.index.names = ['Contig', 'Start', 'End', 'Strand', 'Length', 'PID',
                                         'Gene', 'Synonym', 'Code', 'COG', 'Product']
    df_aggregate_resorted.columns = [sample]
    df_aggregate_resorted.to_csv(genes_summary_file, sep='\t')
    return df_aggregate_resorted


def insertion_falls_in_gene(start, end, strand, insertion, disruption):
    '''return True if the insertion disrupts the 5'-end of the gene defined by disruption (0.0-1.0)'''
    if insertion >= start:
        if insertion <= end:
            # 0.0 = 5'end ; 1.0 = 3'end
            # TODO(Should featureEnd have +1 added?)
            if strand == '+':
                three_primeness = (insertion - start) / (end - start)
            if strand == '-':
                three_primeness = (end - insertion) / (end - start)
            return three_primeness <= disruption


def parse_genbank_setup_bowtie(gbkfile, organism, genomeDir):
    logger.info('Preparing nucleotide fasta file from GenBank file to use in bowtie mapping.\n' \
        '  GenBank source file: {}'.format(gbkfile))
    gbk2fna(gbkfile, organism, genomeDir)
    logger.info('Preparing feature table file from GenBank file to use in gene mapping.\n' \
        '  GenBank source file: {}'.format(gbkfile))
    gbk2ftt(gbkfile, organism, genomeDir)
    # Change directory, build bowtie indexes, change directory back
    with cd(genomeDir):
        # TODO(pass settings object)
        logger.info('Building bowtie index files in {0}'.format('''need to pass settings.genome'''))
        bowtie_build(organism)


def pipeline_analysis():

    print('\n===================='
          '\n*     Analysis     *'
          '\n====================\n')

    # write samples.yml with data for each sample
    print('Writing file with summary data for each sample:\n  {}'.format(Settings.samples_yaml))
    print(yaml.dump(Settings.samplesDict, default_flow_style=False))
    with open(Settings.samples_yaml, 'w') as fo:
        fo.write(yaml.dump(Settings.samplesDict, default_flow_style=False))

    # write summary.yml with data for entire experiment
    print('Writing file with overall summary information:\n  {}'.format(Settings.summary_yaml))
    print(yaml.dump(Settings.summaryDict, default_flow_style=False))
    with open(Settings.summary_yaml, 'w') as fo:
        fo.write(yaml.dump(Settings.summaryDict, default_flow_style=False))


def main(args):
    logger.info('Process command line arguments')
    args = parseArgs(args)
    # Initialize the settings object
    settings = Settings(args.experiment)
    # Keep intermediate files
    settings.keepall = args.keepall
    gbkfile = args.genome
    reads = args.input
    disruption = set_disruption(float(args.disruption))
    # Organism reference files called 'genome.fna' etc
    organism = 'genome'
    # sample names and paths
    samples = args.samples
    barcodes_present = not args.nobarcodes
    if samples:
        samplesDict = tab_delimited_samples_to_dict(samples)
    else:
        reads = os.path.abspath(reads)
        samplesDict = directory_of_samples_to_dict(samples)
    logger.debug('samplesDict: {0}'.format(samplesDict))

    # --- SET UP DIRECTORIES --- #
    createExperimentDirectories(settings.experiment)
    # put 'organism' in summaryDict

    # --- SET UP BOWTIE --- #
    parse_genbank_setup_bowtie(gbkfile, organism, settings.genome_path)

    # --- DEMULTIPLEX, TRIM, AND MAP WITH BOWTIE --- #
    logger.info('Demultiplex reads')
    demultiplex_fastq(reads, samplesDict, settings)

    logger.info('Process INSeq samples')
    insertionDict = extract_chromosome_sequence(reads)
    logger.info('Writing insertion count files')
    write_individual_insertion_count_dictionaries(insertionDict, samplesDict, settings)
    logger.info('Mapping with bowtie')
    map_raw_reads(settings, samplesDict, insertionDict)

    # --- PROCESS BOWTIE MAPPINGS --- #
    logger.info('Process bowtie results')
    process_bowtie_results(settings, samplesDict, insertionDict)
    logger.info('Map to genes')
    process_gene_counts(settings, samplesDict, disruption)
    exit()

    ### TODO(REDO Settings...) ###
    if not samples:
        Settings.summaryDict['total reads'] = 0
        for sample in Settings.samplesDict:
            print(Settings.samplesDict[sample])
            Settings.summaryDict['total reads'] += Settings.samplesDict[sample]['reads_with_bc']

    # --- ANALYSIS OF RESULTS --- #
    pipeline_analysis()

    # --- CONFIRM ETION --- #
    print('\n===================='
          '\n*       Done       *'
          '\n====================\n')


if __name__ == '__main__':
    main()

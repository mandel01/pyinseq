#!/usr/bin/env python3
"""
Demultiplexes a FASTQ file into multiple files by 5' end barcode.

Output path includes the experiment and sample name:
(pysinseq)/results/{experiment}/{sample}.fastq

"""

import collections
import csv
import gzip
import logging
import regex as re
import screed
import sys
from .config import transposon_end
from .utils import convert_to_filename

logger = logging.getLogger(__name__)


def demultiplex_fastq(reads, samplesDict, settings):
    '''Demultiplex a fastq input file by 5' barcode into separate files.

       Use regex to identify the chromosome slice and save this in the read
       record as 'trim', e.g., read['trim'] = (4, 21)
       Save raw reads into '{experiment}/raw_data/{sampleName}.fastq'
       Save trimmed reads into '{experiment}/{sampleName}_trimmed.fastq'
    '''
    # Dictionary of lists to hold FASTQ reads until they are written to files
    # Keys are barcodes
    demultiplex_dict = {}
    for sample in samplesDict:
        bc = samplesDict[sample]['barcode']
        demultiplex_dict[bc] = []
    demultiplex_dict['other'] = []  # unassigned barcodes
    # count of reads
    nreads = 0
    nreads_dict = demultiplex_dict.copy()  # Think about collecting info about number of reads with/out good tn, etc...
    # For each line in the FASTQ file:
    #   Assign to barcode (fastq record into the dictionary)
    #   Cache by barcode; write to the appropriate output files (untrimmed and trimmed)
    pattern = re.compile('''
    ^                               # beginning of string
    ([NACGT][ACGT]{{13,14}}(?:TA))  # group(1) = 16-17 bp of chromosomal sequence
                                    # first bp can be N
                                    # last two must be TA for transposon
    ({left}|{right})                # group(2) flanking transposon sequence (left, right)
    '''.format(left=transposon_end['left'], right=transposon_end['right']), re.VERBOSE)
    with screed.open(reads) as seqfile:
        for read in seqfile:
            barcode = read.sequence[0:4]
            if barcode in demultiplex_dict:
                # barcode is in sample list == of interest
                m = re.search(pattern, read.sequence[4:])
                try:
                    # good transposon sequence
                    chrom_seq, tn_side = m.group(1), m.group(2)
                    read['trim'] = m.span(1)  # trim slice
                    read['tn_side'] = 'left' if tn_side == transposon_end['left'] else 'right'
                except AttributeError:
                    # no intact transposon; demultiplex but do not trim/map
                    read['trim'] = None
                demultiplex_dict[barcode].append(read)
            else:
                # barcode is not in sample list >> '_other.fastq' file
                demultiplex_dict['other'].append(read)
            # Every 5E6 sequences write and clear the dictionary
            nreads += 1
            if nreads % 5E6 == 0:
                logger.info('Demultiplexed {:,} samples'.format(nreads))
                write_reads(demultiplex_dict, samplesDict, settings)
                write_trimmed_reads(demultiplex_dict, samplesDict, settings)
                # Clear the dictionary after writing to file
                for sampleName in demultiplex_dict:
                    demultiplex_dict[sampleName] = []
    write_reads(demultiplex_dict, samplesDict, settings)
    write_trimmed_reads(demultiplex_dict, samplesDict, settings)
    logger.info('Total records demultiplexed: {:,}'.format(nreads))
    return nreads


def write_reads(demultiplex_dict, samplesDict, settings):
    '''Write the fastq data to the correct (demultiplexed) file.'''
    # inverting the value: key pairs to barcode: sample, and also adding 'other': '_other'
    barcode_dict = {'other': '_other'}
    for sample in samplesDict:
        barcode_dict[samplesDict[sample]['barcode']] = sample
    for barcode in demultiplex_dict:
        if demultiplex_dict[barcode]:
            with open('{path}raw_data/{sample}.fastq'.format(
                path=settings.path,
                sample=barcode_dict[barcode]), 'a') as fo:
                for read in demultiplex_dict[barcode]:
                    fo.write('@{n}\n{s}\n+\n{q}\n'.format(
                        n=read.name,
                        s=read.sequence,
                        q=read.quality))


def write_trimmed_reads(demultiplex_dict, samplesDict, settings):
    '''Write the fastq data to the correct (demultiplexed) file.'''
    # inverting the value: key pairs to barcode: sample. Exclude 'other' here
    barcode_dict = {}
    for sample in samplesDict:
        barcode_dict[samplesDict[sample]['barcode']] = sample
    for barcode in demultiplex_dict:
        if barcode != 'other':
            # Write all of the TnL (left side of the transposon) reads
            if transposon_end['left'] == transposon_end['right']:
                fo = '{path}/{sample}_trimmed.fastq'.format(
                    path=settings.path,
                    sample=barcode_dict[barcode])
            else:
                fo = '{path}/{sample}_TnL_trimmed.fastq'.format(
                    path=settings.path,
                    sample=barcode_dict[barcode])
            with open(fo, 'a') as fo:
                for read in demultiplex_dict[barcode]:
                    if read.trim and read.tn_side == 'left':
                        fo.write('@{n}\n{s}\n+\n{q}\n'.format(
                            n=read.name,
                            s=read.sequence[slice(read.trim[0] + 4, read.trim[1] + 4)],
                            q=read.quality[slice(read.trim[0] + 4, read.trim[1] + 4)]))
            # Write all of the TnR (right side of the transposon) reads
            if transposon_end['left'] != transposon_end['right']:
                fo = '{path}/{sample}_TnR_trimmed.fastq'.format(
                    path=settings.path,
                    sample=barcode_dict[barcode])
                with open(fo, 'a') as fo:
                    for read in demultiplex_dict[barcode]:
                        if read.trim and read.tn_side == 'right':
                            fo.write('@{n}\n{s}\n+\n{q}\n'.format(
                                n=read.name,
                                s=read.sequence[slice(read.trim[0] + 4, read.trim[1] + 4)],
                                q=read.quality[slice(read.trim[0] + 4, read.trim[1] + 4)]))



def main():
    '''Start here.'''

if __name__ == '__main__':
    main()

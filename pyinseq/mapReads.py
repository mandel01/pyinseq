#!/usr/bin/env python
import pyinseq.config # config file
import re
import subprocess


def bowtie_build(organism):
    '''Build a bowtie index given a fasta nucleotide file.'''
    fna = organism + '.fna'
    subprocess.check_call([pyinseq.config.bowtieBuild, '-q', fna, organism])


def bowtie_map(genome_index_path, reads, bowtie_output):
    '''Map comma-separated reads to a bowtie index.'''
    # String version of the shell command
    bashCommand = '{0} -m 1 --best --strata -a --fullref --suppress 1,6,7 -n 1 -l 17 {1} -c {2} {3} -p 2' \
        .format(pyinseq.config.bowtie, genome_index_path, reads, bowtie_output)
    # Convert bash command to run properly - no spaces; instead list of entries
    # that will appear in the shell as space-separated
    proc = subprocess.Popen(bashCommand.split(' '), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    bowtie_msg_out = proc.communicate()[0]
    # decode from bytes to unicode
    bowtie_msg_out = bowtie_msg_out.decode('utf8')
    print(bowtie_msg_out)
    return bowtie_msg_out


def parse_bowtie(bowtie_message):
    '''Parse bowtie stdout summary messages into a dictionary.'''
    bowtie_msg_dict = {}
    for line in bowtie_message.split('\n'):
        # extract counts from bowtie printing
        m = re.search('^(\#.+:) (\d+)', line)
        try:
            message, count = m.group(1), int(m.group(2))
            bowtie_msg_dict[message] = count
        except:
            pass
    return bowtie_msg_dict


# ===== Start here ===== #

def main():
    pass

if __name__ == '__main__':
    main()

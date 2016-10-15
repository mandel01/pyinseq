#!/usr/bin/env python
import csv
from collections import OrderedDict
#from io import StringIO
import pyinseq
import numpy as np
import pandas as pd

# pyinseq.utils

def test_filename_replace_spaces():
    assert pyinseq.utils.convert_to_filename('file name') == 'file_name'

def test_filename_leading_whitespace():
    assert pyinseq.utils.convert_to_filename('  filename') == 'filename'

def test_filename_trailing_whitespace():
    assert pyinseq.utils.convert_to_filename('filename  ') == 'filename'

def test_filename_funny_characters():
    assert pyinseq.utils.convert_to_filename('filename*&*&*$$') == 'filename'
    assert pyinseq.utils.convert_to_filename('file*&*&*$$name') == 'filename'


# pyinseq.pyinseq

def test_class_Settings_init():
    x = pyinseq.Settings('example')
    assert x.path == 'results/example/'
    assert x.genome_path == 'results/example/genome_lookup/'
    assert x.raw_path == 'results/example/raw_data/'
    assert x.samples_yaml == 'results/example/samples.yml'
    assert x.summary_yaml == 'results/example/summary.yml'
    assert x.keepall == False
    assert x.barcode_length == 4

def test_set_disruption():
    assert pyinseq.set_disruption(1.0) == 1.0
    assert pyinseq.set_disruption(0.8) == 0.8
    assert pyinseq.set_disruption(0.3) == 0.3
    assert pyinseq.set_disruption(0.0) == 0.0 # should probably return an error here though
    assert pyinseq.set_disruption(-1.0) == 1.0
    assert pyinseq.set_disruption(2.0) == 1.0

def test_tab_delimited_samples_to_dict_trailing_newline():
    #s = StringIO('sample_1\tAAAA\tsample_2\tTTTT\n')
    s = 'pyinseq/tests/testdata/sample01_01.txt'
    assert pyinseq.tab_delimited_samples_to_dict(s) == \
        OrderedDict([('sample_1', {'barcode': 'AAAA'}), ('sample_2', {'barcode': 'TTTT'})])

def test_tab_delimited_samples_to_dict_no_trailing_newline():
    #s = StringIO('sample_1\tAAAA\tsample_2\tTTTT\n')
    s = 'pyinseq/tests/testdata/sample01_02.txt'
    assert pyinseq.tab_delimited_samples_to_dict(s) == \
        OrderedDict([('sample_1', {'barcode': 'AAAA'}), ('sample_2', {'barcode': 'TTTT'})])

def test_tab_delimited_samples_to_dict_no_trailing_newline():
    #s = StringIO('sample_1\tAAAA\tsample_2\tTTTT\n')
    s = 'pyinseq/tests/testdata/sample01_02.txt'
    assert pyinseq.tab_delimited_samples_to_dict(s) == \
        OrderedDict([('sample_1', {'barcode': 'AAAA'}), ('sample_2', {'barcode': 'TTTT'})])

'''def test_map_counted_insertions_to_genes():
    df_sample = pd.DataFrame({'chromosome': [1, 1, 1, 1],
                              'insertion_nucleotide': [140, 180, 250, 500],
                              'count': [14, 600, 1000, 700]})
    df_genome_ftt = pd.DataFrame({'chromosome': [1, 1, 1, 2],
                                  'start': [100, 300, 500, 100],
                                  'end': [200, 600, 900, 250]})
    assert pyinseq.map_counted_insertions_to_genes(df_sample, df_genome_ftt).equals(
        pd.DataFrame({'sample_name': {(1, 100, 200): 614,
                                      (1, 300, 600): 700,
                                      (1, 500, 900): 700,
                                      (2, 100, 250): 0}}))'''

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

def test_insertion_falls_in_gene():
    # different disruption values; note difference in + vs - orientation
    assert pyinseq.insertion_falls_in_gene(1000, 3999, '+', 1600, 1.0)
    assert pyinseq.insertion_falls_in_gene(1000, 3999, '+', 1600, 0.5)
    assert not pyinseq.insertion_falls_in_gene(1000, 3999, '+', 1600, 0.1)
    assert pyinseq.insertion_falls_in_gene(1000, 3999, '-', 1600, 1.0)
    assert not pyinseq.insertion_falls_in_gene(1000, 3999, '-', 1600, 0.5)
    assert not pyinseq.insertion_falls_in_gene(1000, 3999, '-', 1600, 0.1)
    # margins of the disruption value
    assert pyinseq.insertion_falls_in_gene(1000, 3999, '+', 1599, 0.2)
    assert not pyinseq.insertion_falls_in_gene(1000, 3999, '+', 1600, 0.2)
    # margins of gene
    assert pyinseq.insertion_falls_in_gene(1000, 3999, '+', 1000, 1.0)
    assert pyinseq.insertion_falls_in_gene(1000, 3999, '-', 1000, 1.0)
    assert pyinseq.insertion_falls_in_gene(1000, 3999, '+', 3999, 1.0)
    assert pyinseq.insertion_falls_in_gene(1000, 3999, '-', 3999, 1.0)
    # outside of gene
    assert not pyinseq.insertion_falls_in_gene(1000, 3999, '+', 999, 1.0)
    assert not pyinseq.insertion_falls_in_gene(1000, 3999, '+', 4000, 1.0)
    assert not pyinseq.insertion_falls_in_gene(1000, 3999, '-', 999, 1.0)
    assert not pyinseq.insertion_falls_in_gene(1000, 3999, '-', 4000, 1.0)

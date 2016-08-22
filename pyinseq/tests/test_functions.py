#!/usr/bin/env python
import csv
from collections import OrderedDict
# from io import StringIO
import pyinseq


# pyinseq.utils

def test_filename_replace_spaces():
    assert pyinseq.utils.convert_to_filename('file name') == 'file_name'

def test_filename_leading_whitespace():
    assert pyinseq.utils.convert_to_filename('  filename') == 'filename'

def test_filename_trailing_whitespace():
    assert pyinseq.utils.convert_to_filename('filename  ') == 'filename'

# pyinseq.pyinseq

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

def test_set_disruption():
    assert pyinseq.set_disruption(1.0) == 1.0
    assert pyinseq.set_disruption(0.8) == 0.8
    assert pyinseq.set_disruption(0.3) == 0.3
    assert pyinseq.set_disruption(0.0) == 0.0 # should probably return an error here though
    assert pyinseq.set_disruption(-1.0) == 1.0
    assert pyinseq.set_disruption(2.0) == 1.0

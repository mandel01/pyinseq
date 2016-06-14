#!/usr/bin/env python
import csv
import pytest
from collections import OrderedDict
from io import StringIO
from pyinseq import temp, utils

# pyinseq.utils

def test_filename_replace_spaces():
    assert utils.convert_to_filename('file name') == 'file_name'

def test_filename_leading_whitespace():
    assert utils.convert_to_filename('  filename') == 'filename'

def test_filename_trailing_whitespace():
    assert utils.convert_to_filename('filename  ') == 'filename'

# pyinseq.pyinseq

## CORRECT THIS FROM pyinseq.temp once all is set

def test_tab_delimited_samples_to_dict_trailing_newline():
    #s = StringIO('sample_1\tAAAA\tsample_2\tTTTT\n')
    s = 'pyinseq/test/testdata/sample01_01.txt'
    assert temp.tab_delimited_samples_to_dict(s) == \
        OrderedDict([('sample_1', {'barcode': 'AAAA'}), ('sample_2', {'barcode': 'TTTT'})])

def test_tab_delimited_samples_to_dict_no_trailing_newline():
    #s = StringIO('sample_1\tAAAA\tsample_2\tTTTT\n')
    s = 'pyinseq/test/testdata/sample01_02.txt'
    assert temp.tab_delimited_samples_to_dict(s) == \
        OrderedDict([('sample_1', {'barcode': 'AAAA'}), ('sample_2', {'barcode': 'TTTT'})])


'''
# pyinseq.pyinseq

def test_set_disruption():
    assert pyinseq.set_disruption(1.0) == 1.0
'''

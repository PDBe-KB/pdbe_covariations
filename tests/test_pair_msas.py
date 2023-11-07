import os
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch
import tempfile

import multiprocessing
from multiprocessing import Pool

import pytest

from pdbe_covariations.utils.exceptions import CovariationsException
from pdbe_covariations import covariations, msa_utils, pair_msas
from pdbe_covariations.utils import path_utils, execute_command_utils



class TestPairMsas(TestCase):
        
    def test_parse_am3(self):
        # Test that the helper function returns data when a3m input exists                                           
        seq,lab = pair_msas.parse_a3m("./tests/data/mock_paired1.a3m")
        result_seq = ['MLRLLLRHH', '------------MRL------------VVLG-------EQLVD']
        result_lab = ['F5HCP3','UniRef100_D5KLH4']

        result = (result_seq,result_lab)
        
        self.assertEqual(result, (seq, lab))

    def test_parse_am3_not_exist(self):
        with self.assertRaises(FileNotFoundError):
            pair_msas.parse_a3m("./tests/data/invalid.xml")
    def test_uni2idx(self):
        ids = ['A0A3X9BQ16',
               'A0A5C2LZ81']

        hash = pair_msas.uni2idx(ids)
        result = [12946395570, 13563044281]
        
        self.assertEqual((result[0],result[1]), (hash[0],hash[1]))
        



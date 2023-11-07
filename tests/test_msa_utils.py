import os
from pathlib import Path
from unittest.mock import patch
import tempfile

import Bio
from Bio import SeqIO

import multiprocessing
from multiprocessing import Pool

import pytest
from pdbe_covariations.utils.exceptions import CovariationsException
from pdbe_covariations import covariations, msa_utils
from pdbe_covariations.utils import path_utils, execute_command_utils



def test_run_hhblits(args):
    # mock out_file content                                                                                                                          
    file_id = "foo"

    db=os.path.join('tests','data','pdbe_demo_db','pdbe_demo_db')
    out_file = os.path.join(args.out, path_utils.msa_file(file_id))
    with open(out_file, "w") as fp:
        fp.write("one\ntwo")

    input=str(args.input[0])
    with patch.object(execute_command_utils,"execute_command"):
        res = msa_utils.run_hhblits(input, file_id, args.out,db, args.threads)
        assert out_file == res


def test_run_hhfilter(args):
    # mock out_file content                                                                                                                          
    file_id = "foo"
    file_name = path_utils.msa_filtered_file(
        file_id, covariations.IDENTITY, covariations.COVERAGE
    )
    out_file = os.path.join(args.out, file_name)

    with open(out_file, "w") as fp:
        fp.write("one\ntwo")

    with patch.object(msa_utils, "execute_command"):
        res = msa_utils.run_hhfilter(file_id, args.out)
        assert out_file == res


def test_run_hmmbuild(args):
    # mock out_file content                                                                                  
    file_id = "F5HCP3"
    
    with tempfile.TemporaryDirectory() as temp_dir:
        a3m_file = os.path.join(temp_dir,'F5HCP3_id90cov75.a3m')
        with open(a3m_file,"w") as a3m:
            a3m.write(">F5HCP3\n MLRLLL")
        
        out_file = os.path.join(temp_dir, path_utils.msa_hmm_file(file_id))
        with open(out_file, "w") as fp:
            fp.write("one\ntwo")

        with patch.object(execute_command_utils,"execute_command"):
            res = msa_utils.run_hmmbuild(file_id, temp_dir)
        
        assert out_file == res



def test_run_hmmsearch(args):
    unp_id = 'F5HCP3'

    with tempfile.TemporaryDirectory() as out_dir:
        threads=args.threads
        db_hmmer = os.path.join(out_dir,'uniprot_trembl.fasta')
        with open(db_hmmer,"w") as db:
            db.write('>F5HCP3\n MLRLLL')
    
        hmm_file = os.path.join(out_dir,'{}_cov.hmm'.format(unp_id))
        with open(hmm_file,"w") as hmm:
            hmm.write("HMMER3/f [3.3.2 | Nov 2020]\nNAME  F5HCP3_id90cov75\nLENG  214")

        out_file = os.path.join(out_dir, path_utils.msa_sto_file(unp_id))
        with open(out_file, "w") as fp:
            fp.write("one\ntwo")

        with patch.object(execute_command_utils,"execute_command"):
            res = msa_utils.run_hmmsearch(unp_id,db_hmmer,out_dir,threads)

        assert out_file == res

def test_run_convert_sto_a3m(args):
    unp_id = 'F5HCP3'
    hhlib_path =''
    with tempfile.TemporaryDirectory() as out_dir:        
        sto_file = os.path.join(out_dir,'{}_cov.sto'.format(unp_id))
        with open(sto_file,"w") as sto:
            sto.write("# STOCKHOLM 1.0\n #=GF ID F5HCP3_id90cov75 \n #=GF AU hmmsearch (HMMER 3.3.2)")
        out_file = os.path.join(out_dir, path_utils.msa_sto_to_a3m_file(unp_id))
        with open(out_file, "w") as fp:
            fp.write("one\ntwo")

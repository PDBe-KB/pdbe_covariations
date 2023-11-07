import os
from pathlib import Path
from unittest import TestCase
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



class TestCovariations:
        
    def test_execute_command_parallel_with_log(self,tmpdir):
        tmp_file = os.path.join(tmpdir, "logfile.txt")
        inputs=[("ls", tmp_file),("ls", tmp_file)]
        #covariations.execute_command_parallel(inputs)
        execute_command_utils.execute_command_parallel(inputs)
    
        with open(tmp_file, "r") as fp:
            lines = fp.read().splitlines()
            assert len(lines) > 0


    def test_execute_command_with_log(self,tmpdir):
        tmp_file = os.path.join(tmpdir, "logfile.txt")
        input=["ls", tmp_file]
        #covariations.execute_command(input)
        execute_command_utils.execute_command(input)
        with open(tmp_file, "r") as fp:
            lines = fp.read().splitlines()
            assert len(lines) > 0


    def test_execute_command_no_log(self):
        input=["ls", None]
        covariations.execute_command(input)
        execute_command_utils.execute_command(input)
        
        print()


    def test_with_empty_args(self):
        """
        User passes no args, should produce a usage statement and then
        raise SystemExit. Usage statement will appear
        """
        parser = covariations.create_parser()
        with pytest.raises(SystemExit):
            parser.parse_args()


    def test_args_ok(self,fasta_file,uniclust_db,tmpdir):
        """
        User pases OK arguments
        """
        parser = covariations.create_parser()
        unp_id = 'ABCD'
        args = parser.parse_args(
            ["-i", str(fasta_file),"--unp_ids", unp_id, "-d", uniclust_db, "-t 10", "-o", str(tmpdir)]
        )

        assert args.input 
        assert args.unp_ids 
        assert args.out == tmpdir
        assert args.threads == 10
        assert args.db
        assert args.mode == "all"


    def test_args_nok_missing_input(self,uniclust_db, tmpdir):
        """
        User pases not OK arguments
        """
        parser = covariations.create_parser()

        with pytest.raises(SystemExit):
            parser.parse_args(["-d", uniclust_db, "-t 10", "-o", str(tmpdir)])


    def test_args_nok_wrong_treads(self,fasta_file, uniclust_db, tmpdir):
        """
        User pases not OK arguments
        """
        parser = covariations.create_parser()

        with pytest.raises(SystemExit):
            parser.parse_args(
                ["-i", str(fasta_file), "-d", uniclust_db, "-t 0", "-o", str(tmpdir)]
            )

    #def test_run_gremlin(args):
    #    with patch.object(covariations, "execute_command"):
    #        covariations.run_gremlin("foobar", args.out, args.threads)


    def test_get_covariation_pairs_one_id(self,args):
        test_data_dir = os.path.join("tests", "data")
        probability_file = Path(test_data_dir, "F5HCP3_prob.txt")
        score_file = Path(test_data_dir, "F5HCP3_score.txt")
        input_fasta = Path(test_data_dir, "F5HCP3.fasta")
        unp_ids = ["F5HCP3"]
    
        file_paths = (probability_file, score_file)
        fasta_sequence = list(SeqIO.parse(input_fasta, "fasta"))
        name, sequence = fasta_sequence[0].id, str(fasta_sequence[0].seq)
        lenseq =len(sequence)
    
        with patch.object(covariations, "run_gremlin", return_value=file_paths):
            pdbe_covariations = covariations.get_covariation_pairs(lenseq,unp_ids,input_fasta,sequence,args.out, args.threads)
            print (pdbe_covariations)
            assert pdbe_covariations
            assert len(pdbe_covariations) > 100
        
            pivot = pdbe_covariations[0]
            assert pivot == ["F5HCP3",24,"LEU","F5HCP3",39,"PRO",-0.0045671,0.539227]

    def test_get_covariation_pairs_two_ids(self,args):
        test_data_dir = os.path.join("tests", "data")
        probability_file = Path(test_data_dir, "F5HCP3_Q38DE2_prob.txt")
        score_file = Path(test_data_dir, "F5HCP3_Q38DE2_score.txt")
        input_msa1 = Path(test_data_dir, "F5HCP3.fasta")
        input_file = Path(test_data_dir, "paired.a3m")
        input_fasta = Path(test_data_dir, "paired.fasta")
        unp_ids = ["F5HCP3","Q38DE2"]

        file_paths = (probability_file, score_file)
        fasta_seq_msa1 = list(SeqIO.parse(input_msa1, "fasta"))
        name_msa1, seq_msa1 = fasta_seq_msa1[0].id, str(fasta_seq_msa1[0].seq)
        lenseq =len(seq_msa1)

        fasta_sequence = list(SeqIO.parse(input_fasta, "fasta"))
        name, sequence = fasta_sequence[0].id, str(fasta_sequence[0].seq)
    
        with patch.object(covariations, "run_gremlin", return_value=file_paths):
            pdbe_covariations = covariations.get_covariation_pairs(lenseq,unp_ids,input_file,sequence,args.out, args.threads)
            print (pdbe_covariations)
            assert pdbe_covariations
            assert len(pdbe_covariations) > 100

            pivot = pdbe_covariations[0]
            assert pivot == ["F5HCP3",1,"MET","Q38DE2",235,"ALA",1.19551e-06,0.837672]

    def test_run_covariations_three_ids(self,args):
        test_data_dir = os.path.join("tests", "data")
        db_hmmer= None
        db = None
        hhlib_path = None

        with tempfile.TemporaryDirectory() as temp_dir:
            fasta_file_1 = os.path.join(temp_dir,'F5HCP3,fasta')
            with open(fasta_file_1,"w") as fasta1:
                fasta1.write(">F5HCP3\n MLRLLL")
            fasta_file_2 = os.path.join(temp_dir,'FOO,fasta')
            with open(fasta_file_2,"w") as fasta2:
                fasta2.write(">FOO\n MLRLLL")
            fasta_file_3 = os.path.join(temp_dir,'FOOO,fasta')
            with open(fasta_file_3,"w") as fasta3:
                fasta3.write(">FOOO\n MLRLLL")

            unp_ids = ["F5HCP3","FOO","FOOO"]
            input_file_list = [fasta_file_1,fasta_file_2,fasta_file_3]

            pytest.raises(TypeError,covariations.run_covariations,unp_ids)

    def test_run_covariations_diff_ids_files(self,args):
        db_hmmer= None
        db = None
        hhlib_path = None

        with tempfile.TemporaryDirectory() as temp_dir:
            fasta_file_1 = os.path.join(temp_dir,'F5HCP3,fasta')
            with open(fasta_file_1,"w") as fasta1:
                fasta1.write(">F5HCP3\n MLRLLL")

            unp_ids = ["F5HCP3","FOO"]
            input_file_list = [fasta_file_1]

            pytest.raises(TypeError,covariations.run_covariations,(unp_ids,input_file_list))


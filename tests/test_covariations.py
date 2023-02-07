import os
from pathlib import Path
from unittest.mock import patch

import Bio
from Bio import SeqIO

import multiprocessing
import pytest
from pdbe_covariations.utils.exceptions import CovariationsException
from pdbe_covariations import covariations
from pdbe_covariations.utils import path_utils

threads = min(min(8, multiprocessing.cpu_count()), 8)

def test_execute_command_with_log(tmpdir):
    tmp_file = os.path.join(tmpdir, "logfile.txt")

    covariations.execute_command("ls", tmp_file)

    with open(tmp_file, "r") as fp:
        lines = fp.read().splitlines()
        assert len(lines) > 0


def test_execute_command_no_log():
    covariations.execute_command("ls")

    print()


def test_check_binaries_ok():
    with patch("subprocess.run"):
        covariations.__check_binaries()


#def test_binaries_not_available():
#    with pytest.raises(CovariationsException):
#        covariations.__check_binaries()


def test_with_empty_args():
    """
    User passes no args, should produce a usage statement and then
    raise SystemExit. Usage statement will appear
    """
    parser = covariations.create_parser()
    with pytest.raises(SystemExit):
        parser.parse_args()


def test_args_ok(fasta_file, uniclust_db, tmpdir):
    """
    User pases OK arguments
    """
    parser = covariations.create_parser()

    args = parser.parse_args(
        ["-i", str(fasta_file), "-d", uniclust_db, "-t 10", "-o", str(tmpdir)]
    )

    assert str(args.input) == str(fasta_file)
    assert args.out == tmpdir
    assert args.threads == 10
    assert args.db
    assert args.mode == "all"


def test_args_nok_missing_input(uniclust_db, tmpdir):
    """
    User pases not OK arguments
    """
    parser = covariations.create_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(["-d", uniclust_db, "-t 10", "-o", str(tmpdir)])


def test_args_nok_wrong_treads(fasta_file, uniclust_db, tmpdir):
    """
    User pases not OK arguments
    """
    parser = covariations.create_parser()

    with pytest.raises(SystemExit):
        parser.parse_args(
            ["-i", str(fasta_file), "-d", uniclust_db, "-t 0", "-o", str(tmpdir)]
        )


def test_run_hhblits(args):
    # mock out_file content
    file_id = "foo"
    out_file = os.path.join(args.out, path_utils.msa_file(file_id))

    with open(out_file, "w") as fp:
        fp.write("one\ntwo")

    with patch.object(covariations, "execute_command"):
        res = covariations.run_hhblits(args.input, file_id, args.out, args.db, args.threads)
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

    with patch.object(covariations, "execute_command"):
        res = covariations.run_hhfilter(file_id, args.out)
        assert out_file == res


def test_run_gremlin(args):
    with patch.object(covariations, "execute_command"):
        covariations.run_gremlin("foobar", args.out, args.threads)


def test_get_covariation_pairs(args):
    test_data_dir = os.path.join("tests", "data")
    probability_file = Path(test_data_dir, "F5HCP3_prob.txt")
    score_file = Path(test_data_dir, "F5HCP3_score.txt")
    input_fasta = Path(test_data_dir, "F5HCP3.fasta")
    
    file_paths = (probability_file, score_file)
    fasta_sequence = list(SeqIO.parse(input_fasta, "fasta"))
    name, sequence = fasta_sequence[0].id, str(fasta_sequence[0].seq)
    
    with patch.object(covariations, "run_gremlin", return_value=file_paths):
        pdbe_covariations = covariations.get_covariation_pairs("F5HCP3", sequence,args.out, args.threads)
        print (pdbe_covariations)
        assert pdbe_covariations
        assert len(pdbe_covariations) > 100

        pivot = pdbe_covariations[0]
        assert pivot == ["F5HCP3",24,"LEU","F5HCP3",39,"PRO",-0.0045671,0.539227]

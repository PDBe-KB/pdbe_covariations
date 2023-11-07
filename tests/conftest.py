import os
from pathlib import Path

import pytest
from pdbe_covariations import covariations


@pytest.fixture(scope="module")
def fasta_file(tmpdir_factory):
    test_file = tmpdir_factory.mktemp("sequence").join("file.fa")
    test_file.write(">1tqn\nABCD")

    return str(test_file)

@pytest.fixture(scope="module")
def uniclust_db(tmpdir_factory):
    db = tmpdir_factory.mktemp("db").join("UniRef30_2020_02_a3m")
    db.write("foobar")

    return str(db)

@pytest.fixture(scope="module")
def fasta_trembl_db(tmpdir_factory):
    db_hmmer = tmpdir_factory.mktemp("db").join("uniprot_trembl.fasta")
    db_hmmer.write("foobar")

    return str(db_hmmer)

@pytest.fixture(scope="module")
def args(fasta_file, uniclust_db,fasta_trembl_db,tmpdir_factory):
    parser = covariations.create_parser()
    tmpdir = str(tmpdir_factory.mktemp("db"))
    unp_id = "foo"
    args = parser.parse_args(
        ["-i", fasta_file,"--unp_ids", unp_id, "-d", uniclust_db, "-t 1", "-o", str(tmpdir),"--db_hmmer",fasta_trembl_db]
    )

    return args


@pytest.fixture(scope="module")
def test_data_dir():
    return os.path.join(Path(__file__).parent, "data")

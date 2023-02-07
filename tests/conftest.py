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
def args(fasta_file, uniclust_db, tmpdir_factory):
    parser = covariations.create_parser()
    tmpdir = str(tmpdir_factory.mktemp("db"))
    args = parser.parse_args(
        ["-i", fasta_file, "-d", uniclust_db, "-t 10", "-o", str(tmpdir)]
    )

    return args


@pytest.fixture(scope="module")
def test_data_dir():
    return os.path.join(Path(__file__).parent, "data")

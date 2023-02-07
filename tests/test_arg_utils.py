import os
from argparse import ArgumentTypeError

import pytest
from cov_pairs.utils import arg_utils


def test_path_exists_ok(test_data_dir):
    ok_file = os.path.join(test_data_dir, "O27725_cov.csv")

    assert arg_utils.path_exist(ok_file)


def test_path_exists_nok(test_data_dir):
    no_file = os.path.join(test_data_dir, "foo")
    with pytest.raises(ArgumentTypeError):
        arg_utils.path_exist(no_file)


def test_check_positive_float_ok():
    assert 4.2 == arg_utils.positive_float("4.2")


@pytest.mark.parametrize("val", ["-1.1", "bagr"])
def test_check_positive_float_nok(val):
    with pytest.raises(ArgumentTypeError):
        arg_utils.positive_float(val)


def test_check_positive_int_ok():
    assert 1 == arg_utils.positive_int("1")


@pytest.mark.parametrize("val", ["1.1", "bagr", "-8"])
def test_check_positive_int_nok(val):
    with pytest.raises(ArgumentTypeError):
        arg_utils.positive_int(val)


def test_uniclust_db_exist(tmpdir):
    db_name = "UniRef30_2020_02_a3m"
    db_file = f"{db_name}.ffdata"
    db_path = tmpdir.mkdir("db")

    check_path = db_path.join(db_name)
    db_file_path = db_path.join(db_file)

    with open(db_file_path, "w") as fp:
        fp.write("content")

    assert arg_utils.check_db(str(check_path))

    with pytest.raises(ArgumentTypeError):
        new_path = db_path.join("UniRef30_2020_06_a3m")
        arg_utils.check_db(new_path)

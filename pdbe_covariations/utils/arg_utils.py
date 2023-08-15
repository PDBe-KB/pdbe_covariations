#!/usr/bin/env python
# software from PDBe: Protein Data Bank in Europe; https://pdbe.org
#
# Copyright 2020 EMBL - European Bioinformatics Institute
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on
# an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied. See the License for the
# specific language governing permissions and limitations
# under the License.

__doc__ = "Helper methods to test and validate application arguments"

import os
from pathlib import Path
from argparse import ArgumentTypeError


def path_exist(file_path):
    """Check valid path

    Args:
        file_path (str): Path to file or directory to be checked

    Raises:
        argparse.ArgumentTypeError: If it is invalid
    """
    fp = Path(file_path)
    if not fp.exists():
        raise ArgumentTypeError(f"Path {file_path} is not a valid file path.")

    return fp


def check_db(db_path):
    """Check valid uniclust db.

    Args:
        db_path (str): Path to the Uniclust db e.g.:
            /scratch/db/UniRef30_2020_02_a3m

    Raises:
        argparse.ArgumentTypeError: If this is not a valid path to the
        Uniclust db
    """
    db_files = os.listdir(Path(db_path).parent)
    db_name = os.path.basename(db_path)

    if not any(db_name == x.split(".")[0] for x in db_files):
        raise ArgumentTypeError(
            "DB needs to be a valid path to the decompressed uniclust db."
        )

    return db_path


def positive_int(nmb):
    """Check if passed number is a positive integer.

    Args:
        nmb (str): String representation of an integer

    Raises:
        argparse.ArgumentTypeError: If argument is not a positive integer

    Returns:
        int: Parsed value.
    """
    try:
        val = int(nmb)
    except ValueError:
        raise ArgumentTypeError("Argument needs to be a positive int value.")

    if val <= 0:
        raise ArgumentTypeError(f"{val} is an invalid positive int value")

    return val


def positive_float(nmb):
    """Check if passed number is a positive integer. This

    Args:
        nmb (str): String representation of a float.

    Raises:
        ArgumentTypeError: If value is not a positive float

    Returns:
        float: Parsed value
    """
    try:
        val = float(nmb)
    except ValueError:
        raise ArgumentTypeError("Argument needs to be a positive float value.")

    if val <= 0.0:
        raise ArgumentTypeError(f"{val} is an invalid positive float value")

    return val

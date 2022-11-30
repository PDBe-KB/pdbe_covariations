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

__doc__ = """
PDBe covariations © 2020 Protein Data Bank in Europe

Calculate covariation pairs from protein sequence.
"""

import argparse
import datetime
import logging
import multiprocessing
import os
import subprocess
import time

import bioint
import numpy
import pandas
import Bio
from Bio import SeqIO
from bioint.utils.exceptions import CovariationsException
from bioint.utils import arg_utils, path_utils

# region mini config
GREMLIN = "gremlin3"
HHBLITS = "hhblits"
HHFILTER = "hhfilter"

IDENTITY = "90"
COVERAGE = "75"
THRESHOLD = 0.5  # Write out only those pairs with probability greated than 0.5

# endregion mini config

# region check args
def __check_binaries():
    """Check if all relevant binaries can be found in path

    Raises:
        Exception: If any of the binaries is not available
    """
    binaries = [HHBLITS, HHFILTER]

    for b in binaries:
        try:
            subprocess.run([b, "-h"], check=True, stdout=subprocess.DEVNULL)
        except FileNotFoundError:
            raise CovariationsException(
                f"Binary {b} was not found in the $PATH environment variable"
            )


def process_args(args, write_out_parameters=True):
    """Preprocess application arguments for their further use.

    Args:
        args (argparse.Namespace): Parsed application arguments

    Raises:
        AttributeError: If any of the binaries that should be used for
        the pipeline are not present as part of the $PATH variable.
    """
    __check_binaries()
    os.makedirs(args.out, exist_ok=True)

    if write_out_parameters:
        logging.info(f"Running covariations pipeline v. {bioint.__version__}")
        logging.info("Settings:")
        for k, v in vars(args).items():
            logging.info(f"  {k:25s}{v}")


# endregion check args


# region run external programs
def execute_command(cmd, log_file=None):
    """Execute command and measure execution time

    Args:
        cmd (list of str): Execution arguments
    """
    command = " ".join(cmd)
    logging.debug(f"Running command: {command}")
    start = time.perf_counter()

    log_file = open(log_file, "w") if log_file else subprocess.DEVNULL
    subprocess.run(cmd, check=True, stderr=log_file, stdout=log_file)

    runtime = time.perf_counter() - start
    runtime_str = datetime.timedelta(seconds=runtime)

    logging.debug(f"Finished in: {runtime_str}.")


def run_hhblits(input_file, unp_id, out, db, threads):
    """Run hhblits software to generate MSA.

    Args:
        input_file (str): Path to the input file with FASTA sequence.
        unp_id (str): File name that is processed (usually UNP id.)
        out (str): Path to the out directory
        db (str): Path to the uniclust db.
        threads (int): No of threads to be used

    Raises:
        Exception: If something goes wrong

    Returns:
        str: Path to the file with the MSA
    """
    out_file = out / path_utils.msa_file(unp_id)
    log_file = out / path_utils.msa_log_file(unp_id)

    command = [
        HHBLITS,
        "-cpu",
        str(threads),
        "-mact",
        "0.35",
        "-maxfilt",
        "100000000",
        "-neffmax",
        "20",
        "-nodiff",
        "-realign_max",
        "10000000",
        "-n",
        "4",
        "-d",
        db,
        "-i",
        str(input_file),
        "-oa3m",
        str(out_file),
        "-o",
        "/dev/null",  # otherwise a file with best aligning sequences is created
    ]

    try:
        execute_command(command, log_file)
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"Generating MSA for input sequence {unp_id} failed with error: {str(e)}"
        )

    if not os.path.isfile(out_file):
        raise IOError("MSA file was not created.")

    with open(out_file, "r") as fp:
        logging.info(f"There are {len(fp.read().splitlines())} sequences in the MSA.")

    return out_file


def run_hhfilter(unp_id, out):
    """Runs hhfilter to filter out hits from MSA given a set of parameters.

    Args:
        unp_id (str): File name that is processed (usually UNP id).
        out (str): Path to the out directory.

    Raises:
        Exception: If hhfilter software fails

    Returns:
        str: Path to the filtered MSA
    """
    input_file = out / path_utils.msa_file(unp_id)
    out_file = out / path_utils.msa_filtered_file(unp_id, IDENTITY, COVERAGE)
    log_file = out / path_utils.msa_filtered_log_file(unp_id)

    command = [
        HHFILTER,
        "-id",
        IDENTITY,
        "-cov",
        COVERAGE,
        "-i",
        str(input_file),
        "-o",
        str(out_file),
    ]

    try:
        execute_command(command, log_file)
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"Filtering MSA for input sequence {unp_id} failed with error: {str(e)}"
        )

    with open(out_file, "r") as fp:
        logging.info(
            f"There are {len(fp.read().splitlines())} sequences in the filtered MSA."
        )

    if not os.path.isfile(out_file):
        raise IOError("Filtered MSA not created.")

    return out_file


def get_covariation_pairs(unp_id,sequence,out,threads):
    """Calculate residue pairs with covariation score and probability
    given the filtered MSA.

    Args:
        unp_id (str): File name that is processed (usually UNP id).
        out (str): Path to the out directory.
        threads (int): No of threads to be used.

    Returns:
        `list of list of str`: Covariation pairs along with their score
        and probability.
    """
    probability_path, score_path = run_gremlin(unp_id, out, threads)

    score = numpy.loadtxt(score_path)
    probability = numpy.loadtxt(probability_path)

    covariation_pairs = []
    sequence_separation = 5  # to exclude short-range predictions
    sequence_length = score.shape[0]  # sequence length

    for i in range(sequence_length):
        for j in range(i + sequence_separation, sequence_length):
            covariation_pairs.append([i + 1, j + 1,sequence[i],sequence[j], score[i, j], probability[i, j]])

    covariation_pairs.sort(key=lambda x: x[3], reverse=True)

    return covariation_pairs


def run_gremlin(unp_id, out, threads):
    """Run Gremlin software to calculate covariation pairs

    Args:
        unp_id (str): File name that is processed (usually UNP id).
        out (str): Path to the out directory.
        threads (int): No of threads to be used.

    Raises:
        Exception: if the process crashes.

    Returns:
        `tuple of str, str`: Path to the probability matrix and score matrix.
    """
    filtered_msa = out / path_utils.msa_filtered_file(unp_id, IDENTITY, COVERAGE)

    prob_path = out / path_utils.probability_file(unp_id)
    score_path = out / path_utils.score_file(unp_id)

    prob_log = out / path_utils.probability_log_file(unp_id)
    score_log = out / path_utils.score_log_file(unp_id)

    score_cmd = [
        GREMLIN,
        "-t",
        str(threads),
        "-i",
        str(filtered_msa),
        "-o",
        str(score_path),
    ]
    prob_cmd = [
        GREMLIN,
        "-t",
        str(threads),
        "-i",
        str(filtered_msa),
        "-o",
        str(prob_path),
        "-R",
        "PROB8",
    ]

    inputs = [(score_cmd, score_log), (prob_cmd, prob_log)]

    for i in inputs:
        try:
            execute_command(i[0], i[1])
        except subprocess.CalledProcessError:
            raise Exception(
                f"Error occured while running GREMLIN for sequence {unp_id}."
            )

    return prob_path, score_path


# endregion run external programs


def get_msa(input_file, unp_id, out, db, threads):
    run_hhblits(input_file, unp_id, out, db, threads)
    run_hhfilter(unp_id, out)


def get_covariation_info(unp_id, sequence, out, threads):
    covariation_pairs = get_covariation_pairs(unp_id,sequence,out,threads)

    df = pandas.DataFrame(
        covariation_pairs, columns=["Residue A", "Residue B", "Score", "Probability"]
    )
    out_file = out / f"{unp_id}_cov.csv"
    df[df["Probability"] >= THRESHOLD].to_csv(out_file)

    logging.info(f"Covariations were written to: {out_file}")


def run_covariations(input_file, out, db, threads, mode):
    """Run covariations

    Args:
        input_file (str): Path to input structure
        out (Path): Path to the out directory
        db (str): Path to the uniclust db
        threads (int): Number of threads to be used for calculation
        mode (str): Mode to be used
    """
    unp_id = os.path.basename(input_file).split(".")[0]
    fasta_sequence = list(SeqIO.parse(input_file, "fasta"))
    name, sequence = fasta_sequence[0].id, str(fasta_sequence[0].seq)
    
    
    logging.info(f"Getting covariations for: {unp_id}")

    if mode == "msa":
        get_msa(input_file, unp_id, out, db, threads)

    elif mode == "cov":
        get_covariation_info(unp_id,sequence, out, threads)
    else:
        get_msa(input_file, unp_id, out, db, threads)
        get_covariation_info(unp_id, sequence, out, threads)

    logging.info("Finished. Time for beer now?")


def create_parser():
    """Set up a parser to get command line options.

    Returns:
         argparse.ArgumentParser parser
    """
    def_no_cpu = min(min(8, multiprocessing.cpu_count()), 8)
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument(
        "-i",
        "--input",
        type=arg_utils.path_exist,
        help="Input sequence in FASTA format.",
        required=True,
    )
    parser.add_argument(
        "-d",
        "--db",
        type=arg_utils.check_db,
        help="Path to the database of clustered sequences.",
        required=True,
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=arg_utils.positive_int,
        default=def_no_cpu,
        help="Number of threads to be used for this calculation.",
        required=False,
    )

    parser.add_argument(
        "-o",
        "--out",
        type=arg_utils.path_exist,
        help="Output directory.",
        required=True,
    )

    parser.add_argument(
        "--debug",
        action="store_true",
        help="Turn on debug information.",
        required=False,
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "-m",
        "--msa",
        action="store_const",
        const="msa",
        dest="mode",
        help="Run MSA only.",
    )
    group.add_argument(
        "-c",
        "--cov",
        action="store_const",
        const="cov",
        dest="mode",
        help="Run covariations calculation from pre-existing MSA.",
    )
    group.add_argument(
        "-a",
        "--all",
        action="store_const",
        const="all",
        dest="mode",
        help="Run full computation.",
    )
    group.set_defaults(mode="all")

    return parser


def main():
    """Application entry point"""
    parser = create_parser()
    args = parser.parse_args()

    lvl = logging.DEBUG if __name__ == "__main__" and args.debug else logging.INFO
    frm = "[%(asctime)-15s]  %(message)s"
    logging.basicConfig(level=lvl, format=frm, datefmt="%a, %d %b %Y %H:%M:%S")

    process_args(args)

    run_covariations(args.input, args.out, args.db, args.threads, args.mode)


if __name__ == "__main__":
    main()

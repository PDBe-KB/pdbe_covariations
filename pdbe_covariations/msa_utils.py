import argparse
import datetime
import logging
import multiprocessing
import os
import subprocess
import time

import pdbe_covariations
import numpy

from pdbe_covariations.utils.exceptions import CovariationsException
from pdbe_covariations.utils import arg_utils, path_utils
from pdbe_covariations.utils.execute_command_utils import execute_command
# region mini config
HHBLITS = "hhblits"
HHFILTER = "hhfilter"
HMMBUILD = "hmmbuild"
HMMSEARCH = "hmmsearch"

IDENTITY = "90"
COVERAGE = "75"
E_VALUE = "1e-10"
THRESHOLD = 0.5  # Write out only those pairs with probability greated than 0.5                    

#HHblits running commands

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
    out_file = os.path.join(out,path_utils.msa_file(unp_id))
    log_file = os.path.join(out,path_utils.msa_log_file(unp_id))


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
        cmd_input = [command, log_file]
        execute_command(cmd_input)
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
    input_file = os.path.join(out,path_utils.msa_file(unp_id))
    out_file = os.path.join(out,path_utils.msa_filtered_file(unp_id, IDENTITY, COVERAGE))
    log_file = os.path.join(out,path_utils.msa_filtered_log_file(unp_id))


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
        cmd_input = [command, log_file]
        execute_command(cmd_input)
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

#HMMER running commands
def run_hmmbuild(unp_id,out):
    """Run hmmer tool to convert a3m MSA format to hmm.

    Args:
        input_file (str): MSA in a3m format.
        out (str): HMM profile

    Raises:
        Exception: If something goes wrong

    Returns:
        str: Path to the file with the HMM profile
    """
    input_file = os.path.join(out,path_utils.msa_filtered_file(unp_id, IDENTITY, COVERAGE))
    output_file = os.path.join(out,path_utils.msa_hmm_file(unp_id))
    log_file = os.path.join(out,path_utils.msa_hmm_log_file(unp_id))

    command = [
        HMMBUILD,
        "--hand",
        "--amino",
        "--informat=a2m",
        str(output_file),
        str(input_file),
    ]
    try:
        cmd_input = [command, log_file]
        execute_command(cmd_input)
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"Generating HMM profile for input a3m MSA {unp_id} failed with error: {str(e)}"
        )

    if not os.path.isfile(output_file):
        raise IOError("HMM profile was not created.")

    return output_file

def run_hmmsearch(unp_id,db_hmmer,out,threads):
    """Calculate MSA searching against uniprot_trembl.fasta with hmmer.

    Args:
        input_file (str): MSA profile in hmm format.
        out (str): MSA in sto format 

    Raises:
        Exception: If something goes wrong

    Returns:
        str: Path to the file with the MSA (sto format)
    """
    input_file = os.path.join(out,path_utils.msa_hmm_file(unp_id))
    output_file = os.path.join(out,path_utils.msa_sto_file(unp_id))
    log_file = os.path.join(out,path_utils.msa_sto_log_file(unp_id))

    command = [
        HMMSEARCH,
        "-E",
        E_VALUE,
        "--domE", 
        E_VALUE,
        "--notextw",
        "--cpu",
        str(threads),
        "-A",
        str(output_file),
        str(input_file),
        str(db_hmmer),
    ]
    try:
        cmd_input = [command, log_file]
        execute_command(cmd_input)
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"Generating MSA with HMMER/uniprot_trembl.fasta for {unp_id} failed with error: {str(e)}"
        )

    if not os.path.isfile(output_file):
        raise IOError("MSA with HMMER/uniprot_trembl.fasta was not created.")

    return output_file

#HHsuite convert MSAs formats

def convert_sto_a3m(unp_id,hhlib_path,out):
    """Covert MSA in 'sto' format to 'a3m' format using HSuite scripts. 

    Args:
        input_file (str): MSA profile in sto format.
        out (str): MSA in a3m format 

    Raises:
        Exception: If something goes wrong

    Returns:
        str: Path to the file with the MSA (a3m format)x
    """
    input_file = os.path.join(out,path_utils.msa_sto_file(unp_id))
    output_file = os.path.join(out,path_utils.msa_sto_to_a3m_file(unp_id))
    log_file = os.path.join(out,path_utils.msa_convert_format_log_file(unp_id))

    command = [
        hhlib_path,
        "sto",
        "a3m",
        str(input_file), 
        str(output_file),
        "-l",
        "10000",
    ]
    try:
        cmd_input = [command, log_file]
        execute_command(cmd_input)
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"Converting sto to a3m format {unp_id} failed with error: {str(e)}"
        )

    #if not os.path.isfile(output_file):
    #    raise IOError("a3m file converted from sto file  was not created.")

    return output_file

def convert_a3m_fasta(unp_id,input_file,hhlib_path,out):
    """Covert MSA in 'a3m' format to 'fasta' format using HSuite scripts.                                                                                                                              Args:                                                                                                                                                                                                   input_file (str): MSA profile in a3m format.                                                                                                                                                       output_file (str): MSA in fasta format
    Raises:
         Exception: If something goes wrong
    Returns:
         str: Path to the file with the MSA (fasta format)
    """
    output_file = os.path.join(out,'paired.fasta')
    log_file = os.path.join(out,path_utils.msa_convert_format_log_file(unp_id))

    command = [
        hhlib_path,
        "a3m",
        "fas",
        str(input_file),
        str(output_file),
        "-l",
        "10000",
    ]
    try:
        cmd_input = [command, log_file]
        execute_command(cmd_input)
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"Converting sto to a3m format {unp_id} failed with error: {str(e)}"
        )

    #if not os.path.isfile(output_file):                                                                                                                                                            
    #    raise IOError("a3m file converted from sto file  was not created.")                                                                                                                        

    return output_file

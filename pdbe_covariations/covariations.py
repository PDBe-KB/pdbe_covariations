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
from subprocess import Popen
import time
from multiprocessing import Pool

import pdbe_covariations
import numpy
import pandas

import Bio
from Bio import SeqIO
from Bio.SeqUtils import seq3

import tempfile


from pdbe_covariations.utils.exceptions import CovariationsException
from pdbe_covariations.utils import arg_utils, path_utils
from pdbe_covariations.utils.execute_command_utils import execute_command
from pdbe_covariations.msa_utils import run_hhblits, run_hhfilter, run_hmmbuild, run_hmmsearch, convert_sto_a3m, convert_a3m_fasta
from pdbe_covariations.pair_msas import run_pair_alignment

# region mini config
GREMLIN = "gremlin3"

IDENTITY = "90"
COVERAGE = "75"
THRESHOLD = 0.5  # Write out only those pairs with probability greated than 0.5



def process_args(args, write_out_parameters=True):
    
    """Preprocess application arguments for their further use.

    Args:
        args (argparse.Namespace): Parsed application arguments

    Raises:
        AttributeError: If any of the binaries that should be used for
        the pipeline are not present as part of the $PATH variable.
    """
    #__check_binaries()
    os.makedirs(args.out, exist_ok=True)

    if write_out_parameters:
        #logging.info(f"Running covariations pipeline v. {pdbe_covariations.__version__}")
        #logging.info("Settings:")
        for k, v in vars(args).items():
            logging.info(f"  {k:25s}{v}")

    
def get_covariation_pairs(len_seq1,unp_ids,input_file,sequence,out,threads):
    """Calculate residue pairs with covariation score and probability
    given the filtered MSA.

    Args:
        len_seq1 : Sequence length. If two sequences are paired, first sequence length.
        unp_id (str): File name that is processed (usually UNP id).
        out (str): Path to the out directory.
        threads (int): No of threads to be used.

    Returns:
        `list of list of str`: Covariation pairs along with their score
        and probability.
    """
    if len(unp_ids)==1:
        unp_id = str(unp_ids[0])
        unp_id_a = unp_id_b = unp_id
    if len(unp_ids)==2:
        unp_id_1 = str(unp_ids[0])
        unp_id_2 = str(unp_ids[1])
        unp_id="{}_{}".format(unp_id_1,unp_id_2)
        
    probability_path, score_path = run_gremlin(unp_id,input_file, out, threads)

    score = numpy.loadtxt(score_path)
    probability = numpy.loadtxt(probability_path)

    covariation_pairs = []
    sequence_separation = 5  # to exclude short-range predictions
    sequence_length = score.shape[0]  # sequence length

        
    for i in range(sequence_length):
        for j in range(i + sequence_separation, sequence_length):
            if(probability[i, j]>=THRESHOLD):
                res1=seq3(sequence[i]).upper()
                res2=seq3(sequence[j]).upper()
                if len(unp_ids)==2:
                   if i+1 > len_seq1 : unp_id_a = unp_id_2
                   else: unp_id_a = unp_id_1
                   if j+1 > len_seq1 : unp_id_b	= unp_id_2
                   else: unp_id_b = unp_id_1
                covariation_pairs.append([unp_id_a, i + 1,res1,unp_id_b,j+1,res2, score[i, j], probability[i, j]])

    covariation_pairs.sort(key=lambda x: x[3], reverse=True)

    return covariation_pairs


def run_gremlin(unp_id,input_file, out, threads):
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
    #filtered_msa = os.path.join(out,path_utils.msa_filtered_file(unp_id, IDENTITY, COVERAGE))

    prob_path = os.path.join(out,path_utils.probability_file(unp_id))
    score_path = os.path.join(out,path_utils.score_file(unp_id))

    prob_log = os.path.join(out,path_utils.probability_log_file(unp_id))
    score_log = os.path.join(out,path_utils.score_log_file(unp_id))

    score_cmd = [
        GREMLIN,
        "-t",
        str(threads),
        "-i",
        str(input_file),
        "-o",
        str(score_path),
    ]
    prob_cmd = [
        GREMLIN,
        "-t",
        str(threads),
        "-i",
        str(input_file),
        "-o",
        str(prob_path),
        "-R",
        "PROB8",
    ]

    inputs = [(score_cmd, score_log), (prob_cmd, prob_log)]

    with Pool() as pool:

        pool.map(execute_command, inputs)


    return prob_path, score_path


# endregion run external programs

def get_msas_hmmer(unp_id_1,unp_id_2, out, db_hmmer, threads):

    inputs_hmmbuild = [(unp_id_1,out),(unp_id_2,out)]
    
    with Pool() as pool:
        pool.starmap(run_hmmbuild,inputs_hmmbuild)

    inputs_hmmsearch = [(unp_id_1,db_hmmer,out,threads),(unp_id_2,db_hmmer,out,threads)]

    with Pool() as pool:
        pool.starmap(run_hmmsearch,inputs_hmmsearch)
        

def get_msa(input_file, unp_id, out, db, threads):
    run_hhblits(input_file, unp_id, out, db, threads)
    run_hhfilter(unp_id, out)

def convert_to_a3m(unp_id_1,unp_id_2,hhlib_path,out):
    inputs = [(unp_id_1,hhlib_path,out),(unp_id_2,hhlib_path,out)]
    with Pool() as pool:
        pool.starmap(convert_sto_a3m,inputs)

def get_covariation_info(len_seq1,unp_ids,input_file, sequence, out, threads):
    covariation_pairs = get_covariation_pairs(len_seq1,unp_ids,input_file,sequence,out,threads)
    
    if len(unp_ids)==2:
        unp_id="{}_{}".format(unp_ids[0],unp_ids[1])
    else:
        unp_id = str(unp_ids[0])
    
    df = pandas.DataFrame(
         
         covariation_pairs, columns=["uniprot_accession_a","uniprot_residue_index_a","uniprot_residue_label_a", "uniprot_accession_b","uniprot_residue_index_b","uniprot_residue_label_b", "covariation_score", "covariation_probability"]
    )
    out_file = os.path.join(out,"{}_cov.csv".format(unp_id))
    
    df[df["covariation_probability"] >= THRESHOLD].to_csv(out_file,index=False)

    logging.info(f"Covariations were written to: {out_file}")

def run_covariations(input_file_list, out, db,db_hmmer, threads, mode,unp_ids,hhlib_path,force=False):
    """Run covariations

    Args:
        input_file_list (str): list of Paths to input sequences in fasta format
        out (Path): Path to the out directory
        db (str): Path to the uniclust db
        db_hmmer : Path to the uniprot_trembl.fasta db  
        threads (int): Number of threads to be used for calculation
        mode (str): Mode to be used
        unp_ids : List of UniProt accession ids
        hhlib_path : path to HHSuite/script
        force : force calculation of MSAs 
    """
    if len(unp_ids)==1 and len(input_file_list)==1:

        input_file = str(input_file_list[0])
        fasta_sequence = list(SeqIO.parse(input_file, "fasta"))
        name, sequence = fasta_sequence[0].id, str(fasta_sequence[0].seq)
        len_seq1=len(sequence)
        
        filtered_msa_file = os.path.join(out,path_utils.msa_filtered_file(unp_id, IDENTITY, COVERAGE))
        msa_file =  os.path.join(out,path_utils.msa_file(unp_id))
    
        logging.info(f"Getting covariations for homomeric complex: {unp_id}")

        
        if force :
            get_msa(input_file, unp_id, out, db, threads)
            get_covariation_info(len_seq1,unp_ids,filtered_msa_file,sequence, out, threads)    
        elif mode == "msa" :
            get_msa(input_file, unp_id, out, db, threads)
        elif mode == "cov":
            get_covariation_info(len_seq1,unp_ids,filtered_msa_file,sequence, out, threads)
        else:

            if os.path.exists(msa_file) and os.path.exists(filtered_msa_file) :
                get_covariation_info(len_seq1,unp_ids,filtered_msa_file, sequence, out, threads)
            else:
                get_msa(input_file, unp_id, out, db, threads)
                get_covariation_info(len_seq1,unp_ids,filtered_msa_file,sequence, out, threads)

        logging.info("Finished. Time for beer now?")
        
    if len(unp_ids)==2 and len(input_file_list)==2 :
        
        unp_id_1 = str(unp_ids[0])
        unp_id_2 = str(unp_ids[1])
        input_file_1=str(input_file_list[0])
        input_file_2=str(input_file_list[1])
        fasta_sequence_1 = list(SeqIO.parse(input_file_1, "fasta"))
        fasta_sequence_2 = list(SeqIO.parse(input_file_2, "fasta"))
        name_1, sequence_1 = fasta_sequence_1[0].id, str(fasta_sequence_1[0].seq)
        name_2, sequence_2 = fasta_sequence_2[0].id, str(fasta_sequence_2[0].seq)
        
        filtered_msa_file_1 = os.path.join(out,path_utils.msa_filtered_file(unp_id_1, IDENTITY, COVERAGE))
        filtered_msa_file_2 = os.path.join(out,path_utils.msa_filtered_file(unp_id_2, IDENTITY, COVERAGE))
        msa_file_1 =  os.path.join(out,path_utils.msa_file(unp_id_1))
        msa_file_2 =  os.path.join(out,path_utils.msa_file(unp_id_2))
        msa_sto_file_1=os.path.join(out,path_utils.msa_sto_file(unp_id_1))
        msa_sto_file_2=os.path.join(out,path_utils.msa_sto_file(unp_id_2))
        msa_sto_to_a3m_file_1=os.path.join(out,path_utils.msa_sto_to_a3m_file(unp_id_1))
        msa_sto_to_a3m_file_2=os.path.join(out,path_utils.msa_sto_to_a3m_file(unp_id_2))
        pair_aligment_file= os.path.join(out,"paired.a3m")

        unp_id="{}_{}".format(unp_id_1,unp_id_2)
        
        if force :
            #Run MSAs calculations with HHblits
            inputs_msas = [(input_file_1, unp_id_1, out, db, threads),
                           (input_file_2, unp_id_2, out, db, threads)]
            with Pool() as pool:
                pool.starmap(get_msa,inputs_msas)

            #Calculate MSAs with HMMER from input MSA obtained with HHblits
            
            get_msas_hmmer(unp_id_1,unp_id_2, out, db_hmmer, threads)
            convert_to_a3m(unp_id_1,unp_id_2,hhlib_path,out)

            #Align MSAs and calculate length of sequences aligned
            
            len_msa1, len_msa2 = run_pair_alignment(msa_sto_to_a3m_file_1,msa_sto_to_a3m_file_2)

            #Calculate covariation pairs
            
            with tempfile.TemporaryDirectory() as tmp_dir:
                paired_fasta = convert_a3m_fasta(unp_id,pair_aligment_file,hhlib_path,tmp_dir)
                fasta_sequence = list(SeqIO.parse(paired_fasta, "fasta"))
                name, sequence = fasta_sequence[0].id, str(fasta_sequence[0].seq)
                
                logging.info(f"Getting covariations for heteromeric pair: {unp_id_1,unp_id_2}")
                get_covariation_info(len_msa1,unp_ids,pair_aligment_file,
                                     sequence, out, threads)
        else:
            #If MSAs exist, skip MSA calculations and compute covariation pairs
            
            if not os.path.exists(msa_sto_file_1) or not os.path.exists(msa_sto_file_2):
                
                if os.path.exists(filtered_msa_file_1) and os.path.exists(filtered_msa_file_2) :
                    
                    get_msas_hmmer(unp_id_1,unp_id_2, out, db_hmmer, threads)
                
                else:
                    inputs_msas=[(input_file_1, unp_id_1, out, db, threads),
                                 (input_file_2, unp_id_2, out, db, threads)]
                    
                    with Pool() as pool:
                        pool.starmap(get_msa,inputs_msas)

                    get_msas_hmmer(unp_id_1,unp_id_2, out, db_hmmer, threads)
                    
            convert_to_a3m(unp_id_1,unp_id_2,hhlib_path,out)
            len_msa1, len_msa2 = run_pair_alignment(msa_sto_to_a3m_file_1,msa_sto_to_a3m_file_2)
            
            with tempfile.TemporaryDirectory() as tmp_dir:
                paired_fasta = convert_a3m_fasta(unp_id,pair_aligment_file,hhlib_path,tmp_dir)
                fasta_sequence = list(SeqIO.parse(paired_fasta, "fasta"))
                name, sequence = fasta_sequence[0].id, str(fasta_sequence[0].seq)

                logging.info(f"Getting covariations for heteromeric pair: {unp_id_1,unp_id_2}")
                get_covariation_info(len_msa1,unp_ids,pair_aligment_file,sequence, out, threads)
        
    if len(unp_ids)>2:
        logging.error("Sorry, at the moment covariation analysis handles only pairs of uniprot ids")
        
    if len(unp_ids)!=len(input_file_list):
        logging.error("Missing input files or uniprot ids, the number of input sequences should match the number of uniprot ids")

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
        nargs="+",
        help="List of input sequence(s) in FASTA format.",
        required=True,
    )
    parser.add_argument(
        "--unp_ids",
        nargs="+",
        help="List of accession ids (maximun 2 uniprot accession ids for current method)",
        required=True,
    )
    parser.add_argument(
        "-d",
        "--db",
        #type=arg_utils.check_db,
        help="Path to the uniclust database of clustered sequences.",
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
        "--db_hmmer",
        help="Path to the database uniprot_trembl.fasta of clustered sequences.",
        #default=argparse.SUPPRESS
        required=False,
    )
    parser.add_argument(
        "--hhlib_path",
        help="Provide path to HHSuite/script to convert .sto to .a3m format",
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

    parser.add_argument(
        "--force",
        action="store_true",
        help="Always run msa calculation with hhblits/hhfilter",
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
    parser.add_argument(
        "-a",
        "--all",
        action="store_const",
        const="all",
        dest="mode",
        help="Run full computation",
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

    db_hmmer= None
    hhlib = None
    
    if len(args.unp_ids) > 1:
        if not args.db_hmmer or not args.hhlib_path :
            logging.error("add flag --db_hmmer and --hhlib_path ,required for covariation pipeline for heteromeric pairs")
        else:
            db_hmmer=args.db_hmmer
            hhlib=args.hhlib_path
            
    run_covariations(args.input, args.out, args.db,db_hmmer, args.threads, args.mode,args.unp_ids,hhlib, args.force)


if __name__ == "__main__":
    main()

# Basic information

This code performs the calculation of covariation pairs from protein sequence. 

The code will: 
1) Perform multiple sequence aligment using HHSuite
2) Calculate covariation pairs using Gremlin3

# Installation

git clone https://github.com/PDBe-KB/covariation_pairs
cd covariation-pairs
python setup.py install

# Dependencies

The module runs the packages listed below as a subprocesses and requires a-priori compilation of these packages:

HHsuite   https://github.com/soedinglab/hh-suite
Gremlin3  https://github.com/gjoni/gremlin3

The modules requires to provide the path to the database used for MSA in HH-suite:

Uniclust - https://uniclust.mmseqs.com/

The process requires environment variables for binaries hhblits and gremlin3 (HHsuite and Gremlin3), which can be set in the path:
PATH="$PATH:/your_path/hh-suite/build/bin:$PATH:/your_path/gremlin3/bin"

# Usage

After installing the package, the 'covariation' module can be run in terminal as:

```
covariation -i fasta_input/MSA_input -d clustered_sequences_database -t no_threads -o output_path
```

Required:
```
-i / --input : Path to the input file with FASTA sequence 
-d /--db     : Path to the database of clustered sequences
-t           :  No. of threads for calculation
-o           : output directory 
```

Optional:

```
--debug  :  Turn on debug information
-m / --msa : run MSA only
-c / --cov : Run covariations calculation from pre-existing MSA
-a / --all : Run full computation (default)
```

The process is as follows:

1. The process first takes as an input a sequence in FASTA format for the uniprot accession 
2. Next the process runs hhblits to perform multiple sequence analysis (MSA) and generates a file:
   name_file.a3m where name_file is the same as the FASTA file
3. Next step, the process runs hhfilter to filter out hits from MSA and generates a new file :
   filtered_file.a3m(unp_id, IDENTITY, COVERAGE)
6. An output 

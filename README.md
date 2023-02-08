Protein residue covariation pipeline
=

## Basic information

The code in this repository runs a pipeline that calculates covariation pairs from protein sequences. 

The main steps of this pipeline are: 
1) Create and filter multiple sequence aligment (MSA) using [HHSuite](https://github.com/soedinglab/hh-suite)
2) Calculate covariation pairs using [Gremlin3](https://github.com/gjoni/gremlin3)

## Installation

```
git clone https://github.com/PDBe-KB/covariation_pairs

cd covariation-pairs

python setup.py install
```

## Dependencies

The process runs the packages listed below as subprocesses and requires a-priori compilation of:

HHsuite   https://github.com/soedinglab/hh-suite

Gremlin3  https://github.com/gjoni/gremlin3

The process requires a path to a database with clustered sequences, used for MSA in HH-suite. We use the Uniclust database by default:

Uniclust - https://uniclust.mmseqs.com/

The process requires environment variables for binaries hhblits, hhfilter and gremlin3 (HHsuite and Gremlin3), which can be set in the path:

```
PATH="$PATH:/your_path/hh-suite/build/bin:$PATH:/your_path/gremlin3/bin"

```
Other dependencies can be installed with:

```
pip install -r requirements.txt
```
See  [requirements.txt](https://github.com/PDBe-KB/covariation_pairs/blob/main/requirements.txt)

For development: 

**pre-commit usage**

```
pip install pre-commit
pre-commit
pre-commit install
```

## Usage

After installing the package, the `covariation` module can be run in terminal as:

```
covariation -i fasta_input/MSA_input -d clustered_sequences_database -t no_threads -o output_path
```

Required:
```
-i / --input : Path to the input file with FASTA sequence or file with MSA  
-d /--db     : Path to the database of clustered sequences
-o           : Output directory 
```

Optional:

```
--debug     : Turn on debug information
-m / --msa  : Run MSA only
-t          : No. of threads for calculation
-c / --cov  : Run covariations calculation from pre-existing MSA
-a / --all  : Run full computation (default)
```

## Overview of the process

1. The process first reads an input file which contains a sequence in FASTA format for a UniProt accession. The input file can also be a pre-existing MSA:
   - UNP_acc.fasta  (the name of the file UNP_ACC must be the UniProt accession number)
   - UNP_acc.a3m (pre-existing MSA file, the file must be named as the UniProt accession number UNP_acc)
   
2. Next, if -c flag is not used,  the process runs `hhblits` to create an MSA and generates a file:
   - UNP_acc.a3m (where the name of the file UNP_acc is the UniProt accession number)
3. Next step, the process runs `hhfilter` to filter out hits from the MSA and generates a new file:
   - UNP_acc_IDENTITY_COVERAGE.a3m (file name: UniProt id, identity, coverage)
6. Then, the process runs `gremlin3` to calculate covariation pairs and outputs two files:
   - UNP_ACC_prob.txt (where UNP_ACC refers to the UniProt id)
   - UNP_ACC_score.txt (where UNP_ACC refers to the UniProt id)
7. Finally, the process reads through scores and probabilities files and outputs a CSV file with the covariation pairs (with probability larger than 0.5):
   - UNP_ACC_cov.csv (UNP_ACC: UniProt id)
   
## Expected output CSV file

The output CSV file looks as follows:
```
uniprot_accession_a,uniprot_residue_index_a,uniprot_residue_label_a,uniprot_accession_b,uniprot_residue_index_b,uniprot_residue_label_b,covariation_score,covariation_probability
F5HCP3,24,LEU,F5HCP3,39,PRO,-0.0045671,0.539227
F5HCP3,24,LEU,F5HCP3,41,TRP,-0.0045671,0.695802
F5HCP3,24,LEU,F5HCP3,46,TYR,-0.0045671,0.569487
F5HCP3,24,LEU,F5HCP3,52,ALA,-0.0045682,0.516928
F5HCP3,24,LEU,F5HCP3,56,TYR,-0.00456711,0.569487
F5HCP3,24,LEU,F5HCP3,57,CYS,-0.00456711,0.671099
```

## Versioning

We use [SemVer](https://semver.org) for versioning.

## Authors
* [Grisell Diaz Leines](https://github.com/grisell) - Developer
* [Lukas Pravda](https://github.com/grisell) - Developer
* [Mihaly Varadi](https://github.com/mvaradi) - Review and management 

See all contributors [here](https://github.com/PDBe-KB/pisa-analysis/graphs/contributors).

## License

See  [LICENSE](https://github.com/PDBe-KB/pisa-analysis/blob/main/LICENSE)

#!/bin/bash
#
#PBS -l select=1:ncpus=10:mem=50gb:scratch_local=150gb
#PBS -l walltime=20:00:00
#PBS -N 05_parse_tables

# Set the database directory (files ending "genomic.fna")
DATABASE_DIR="/full/path/to/database/"

# Set the directory containing results from infernal (files for all models must be in one directory and ending "genomic.csv" and "genomic-alignment")
INFERNAL_DIR="/full/path/to/database/all"

# Set taxonomy.csv 
FILE_TAXONOMY="/full/path/to/database/taxonomy.csv"

# Set the directory for outputs -> the directory is going to be created in the next step
RESULT_DIR="/full/path/to/results_run2/"

# Set the path to the parse_table.py script
PARSE_TABLES="/full/path/to/parse_tables.py"

# Create the directory for outputs
mkdir -p $RESULT_DIR

# Load required modules
module add python36-modules-gcc
source /cvmfs/software.metacentrum.cz/modulefiles/5.1.0/loadmodules
module load seqtk

# Run parse_tables script
python3 $PARSE_TABLES -i $INFERNAL_DIR -t $FILE_TAXONOMY -r 200 -d $DATABASE_DIR -o $RESULT_DIR







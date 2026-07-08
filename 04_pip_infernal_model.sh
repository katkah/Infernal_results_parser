#!/bin/bash
#
# Runs INFERNAL cmsearch across all covariance models × all genomes.
#
# The cmsearch invocation and output-file naming convention are taken from the
# GERONIMO pipeline (rule `infernal_search`) by Agata M. Kilar, reused here under
# the MIT License:
#   Kilar A.M., Fajkus P., Fajkus J. (2023) GERONIMO: a tool for systematic
#   retrieval of structural RNAs in a broad evolutionary context.
#   GigaScience 12, giad080.  https://github.com/amkilar/GERONIMO
#   Copyright (c) Agata M. Kilar — MIT License.

#Script to run infernal
# inputs: 
	# calibrated models
	# database (reference) in fasta format

DIR_MODELS=/full/path/to/models
DIR_DATABASE=/full/path/to/database

#preparation of envirnment for infernal
echo "Preparing environment"
module add conda-modules-py37
conda activate infernal 
chmod 744 $DIR_DATABASE
cd $DIR_MODELS/

#LOOP OVER ALL MODELS PER EACH .FNA IN ./Database
THREADS=$PBS_NUM_PPN

# change cov_* for the specific model name to run single analysis (for example cov_model_lncRNA1)
for model in cov_*
do
 echo "${model} entering the loop"
 SPECIES=${model#"cov_model_"} #drop prefix
 echo "making directory for ${SPECIES}"
 mkdir $DIR_DATABASE/$SPECIES

 for ref in $DIR_DATABASE/*.fna
 do
  echo "${ref} entering the loop"
  REF=$(basename "${ref%.fna}") #drop suffix and path
  RESULT="result_${SPECIES}_vs_${REF}"
  ALIN="result_${SPECIES}_vs_${REF}-alignment"
  TAB="result_${SPECIES}_vs_${REF}.csv"
  
  ################ SEARCH ################
  #search a sequence database with cmsearch
 
  cmsearch --cpu $THREADS --notextw -A $DIR_DATABASE/$SPECIES/$ALIN -o $DIR_DATABASE/$SPECIES/$RESULT --tblout $DIR_DATABASE/$SPECIES/$TAB $DIR_MODELS/$model $ref 

 done

done


################ CLEAN-UP ################
echo "Cleaning-up".

conda deactivate


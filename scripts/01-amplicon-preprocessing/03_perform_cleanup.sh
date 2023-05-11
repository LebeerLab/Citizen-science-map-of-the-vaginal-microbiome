#!/usr/bin/env bash 

# The goal of this script is to perform cleanup of intermediate output files for 
# the current pipeline. 

dio=../../results

# determine pipeline name 
pipeline=$(basename $(pwd))

# perform cleanup if main output file is present and not empty
if [ -s $dio/$pipeline/*tsv ] 
then
  echo MAIN RESULT PRESENT AND NOT EMPTY - PERFORMING CLEANUP
  rm -r $dio/$pipeline/samples_filtered
  rm -r $dio/$pipeline/intermediate
else 
  echo MAIN RESULT ABSENT OR EMPTY - ABORTING CLEANUP
  exit 1
fi

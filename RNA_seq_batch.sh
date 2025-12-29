#!/usr/bin/bash

set -e

# call this script like this:
# $ ./RNA_seq_batch.sh ./path/to/sample/dirs/ SAMPLE_1,SAMPLE_2,SAMPLE_3
# the sample argument should be a comma-separated list with no spaces between names

# the directory in which the sample directories can be found, including the trailing /
directoryPrefix=$1

# split the sample list into an array
IFS=',' read -r -a samples <<< "$2"

# run the readFiltering script for each given sample
for sample in ${samples[@]}
do
  ./RNA_seq_readFiltering.sh "${directoryPrefix}" $sample
done

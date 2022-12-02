#!/bin/bash

#SBATCH --time=02:00:00

### Usage
# sbatch combine_sv_data.sh

# date name of your matrices folder
input_date="2022_11_11"

# define dirs
base_dir=$(dirname $(dirname $(dirname $PWD)))

# input dirs
input_dir=$base_dir/results/11_downsampling/02_extract_mutational_signatures_from_vcfs/${input_date}/results/matrices/sv

# create output dirs if they dont exist
output_file=$base_dir/results/11_downsampling/02_extract_mutational_signatures_from_vcfs/${input_date}/results/matrices/sv_contexts.txt

# init counter
counter=0

# loop over each individual matrix and combine them into one
for i in $input_dir/*txt; do
	let counter=counter+1

	if [[ $counter -eq 1 ]]; then
		head -n1 $i > $output_file
	fi

	tail -n+2 $i >> $output_file
done

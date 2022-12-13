#!/bin/bash

### Usage
# sbatch snv_create_hmf_coverage_info_table.sh

# get current date as variable
current_date=$(date +"%Y_%m_%d")

# define dirs
base_dir=$(dirname $(dirname $PWD))

output_dir=$base_dir/results/07_create_coverage_tables/${current_date}/results; mkdir -p $output_dir

# run script
Rscript $base_dir/code/07_create_coverage_tables/snv_create_hartwig_coverage_info_table.R $output_dir

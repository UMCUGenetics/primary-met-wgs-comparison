#!/bin/bash

#SBATCH --time=04:00:00

### Usage
# sbatch get_fusion_contingency_matrix.sh -l 2021_10_25 -e no

# get command line input that gets passed to the Rscript
while getopts ":l:e:" opt; do
  case ${opt} in
    l )
    linx_fusion_table_results_date=$OPTARG
      ;;
    e )
    exclude_hypermutators=$OPTARG
      ;;
    \? ) echo "Usage: cmd [-d] <linx_fusion_table_results_date> [-e] <yes/no>"; exit 1
      ;;
		: ) echo "Invalid option: $OPTARG requires an argument" 1>&2
	    ;;
  esac
done

# define dirs
output_dir=/path/to/output; mkdir -p $output_dir

# run Rscript
Rscript /path/to/Rscript/get_fusion_contingency_matrix.R \
$linx_fusion_table_results_date $exclude_hypermutators $output_dir

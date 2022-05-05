#!/bin/bash

#SBATCH --time=04:00:00

### Usage
# sbatch get_contingency_matrix.sh -l 2021_10_25 -d 2021_10_26 -e no -f yes -c clonal

# get command line input that gets passed to the Rscript
while getopts ":l:d:e:f:c:" opt; do
  case ${opt} in
    l )
    linx_table_results_date=$OPTARG
      ;;
    d )
    linx_driver_clonality_results_date=$OPTARG
      ;;
    e )
    exclude_hypermutators=$OPTARG
      ;;
    f )
    filter_clonality=$OPTARG
      ;;
    c )
    clonality=$OPTARG
      ;;
    \? ) echo "Usage: cmd [-d] <linx_table_results_date> [-e] <yes/no>"; exit 1
      ;;
		: ) echo "Invalid option: $OPTARG requires an argument" 1>&2
	    ;;
  esac
done

# define dirs
output_dir=/path/to/output; mkdir -p $output_dir

# run Rscript
Rscript /path/to/Rscript/get_contingency_matrix.R \
$linx_table_results_date $linx_driver_clonality_results_date $exclude_hypermutators $filter_clonality $clonality $output_dir

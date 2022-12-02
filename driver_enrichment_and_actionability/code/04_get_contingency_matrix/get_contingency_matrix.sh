#!/bin/bash

#SBATCH --time=04:00:00

### Usage
# WARNING: -c, -s -n, -p & -m are mutually exclusive! Only one at a time can be 'yes'. Or all can be 'no'.
# sbatch get_contingency_matrix.sh -l 2021_12_03 -d 2021_12_03 -f no -c clonal -s yes -n no -p no -m no

# get command line input that gets passed to the Rscript
while getopts ":l:d:f:c:s:n:p:m:" opt; do
  case ${opt} in
    l )
    linx_table_results_date=$OPTARG
      ;;
    d )
    linx_driver_clonality_results_date=$OPTARG
      ;;
    f )
    filter_clonality=$OPTARG
      ;;
    c )
    clonality=$OPTARG
      ;;
    s )
    by_subtype=$OPTARG
      ;;
    n )
    by_subtype_met_location=$OPTARG
      ;;
    p )
    by_progression=$OPTARG
      ;;
    m )
    by_met_location=$OPTARG
      ;;
    \? ) echo "Usage: cmd [-l] <linx_table_results_date> [-d] <linx_clonality_results_date> [-f] <yes/no> [-c] <clonal/subclonal> [-s] <yes/no> [-n] <yes/no> [-p] <yes/no> [-m] <yes/no>"; exit 1
      ;;
		: ) echo "Invalid option: $OPTARG requires an argument" 1>&2
	    ;;
  esac
done

# get current date as variable
current_date=$(date +"%Y_%m_%d")

# define dirs
base_dir=$(dirname $(dirname $PWD))

# define dirs
output_dir=$base_dir/results/04_get_contingency_matrix/${current_date}/results; mkdir -p $output_dir

# run Rscript
Rscript $base_dir/code/04_get_contingency_matrix/get_contingency_matrix.R \
$linx_table_results_date $linx_driver_clonality_results_date \
$filter_clonality $clonality $by_subtype $by_subtype_met_location $by_progression $by_met_location \
$output_dir

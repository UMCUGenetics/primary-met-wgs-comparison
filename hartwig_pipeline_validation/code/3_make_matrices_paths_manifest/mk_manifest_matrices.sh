#!/bin/bash

### Usage
# sbatch mk_manifest_matrices.sh -f <main_context_matrices_folder>
# e.g. sbatch mk_manifest_matrices.sh -f 2021_12_03_matrices
# <matrices_folder> is the folder name where you wrote your matrices to.

# get the length_cutoff as commandline input
while getopts ":f:" opt; do
  case ${opt} in
    f )
		matrices_folder=$OPTARG
      ;;
    \? ) echo "Usage: cmd [-f]"
      ;;
		: ) echo "Invalid option: $OPTARG requires an argument" 1>&2
	    ;;
  esac
done

# define dirs
base_dir=$(dirname $(dirname $PWD))

# run script
Rscript $base_dir/code/3_make_matrices_paths_manifest/mk_manifest_matrices.R $matrices_folder

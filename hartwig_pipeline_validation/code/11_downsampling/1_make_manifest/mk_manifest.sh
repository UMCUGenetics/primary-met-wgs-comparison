#!/bin/bash

### Usage
# sbatch mk_manifest.sh -p <postprocessed_vcf yes/no>
# e.g. sbatch mk_manifest.sh -p yes

### Order of execution
# 1. run sbatch mk_manifest.sh -p no (creates manifest for downsampled samples)
# 2. run sbatch postprocess_vcfs.sh
# 3. run sbatch mk_manifest.sh -p yes (creates manifest for downsampled & postprocessed samples)

# get commandline input
while getopts ":p:" opt; do
  case ${opt} in
    p )
    postprocessed_vcf=$OPTARG
      ;;
    \? ) echo "Usage: cmd [-p]"
      ;;
		: ) echo "Invalid option: $OPTARG requires an argument" 1>&2
	    ;;
  esac
done

# define dirs
base_dir=$(dirname $(dirname $(dirname $PWD)))

# run script
Rscript $base_dir/code/11_downsampling/1_make_manifest/mk_manifest.R $postprocess_vcf

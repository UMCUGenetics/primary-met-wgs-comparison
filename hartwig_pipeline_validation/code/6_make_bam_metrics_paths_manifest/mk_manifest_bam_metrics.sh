#!/bin/bash

### Usage
# sbatch mk_manifest_bam_metrics.sh

# define dirs
base_dir=$(dirname $(dirname $PWD))

# run script
Rscript $base_dir/code/6_make_bam_metrics_paths_manifest/mk_manifest_bam_metrics.R

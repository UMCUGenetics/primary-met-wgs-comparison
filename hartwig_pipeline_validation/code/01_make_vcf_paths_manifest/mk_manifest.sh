#!/bin/bash

### Usage
# sbatch mk_manifest.sh

# define dirs
base_dir=$(dirname $(dirname $PWD))

Rscript $base_dir/code/01_make_vcf_paths_manifest/mk_manifest.R

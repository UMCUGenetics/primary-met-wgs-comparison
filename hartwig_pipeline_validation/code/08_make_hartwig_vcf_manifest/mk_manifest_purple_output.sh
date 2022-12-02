#!/bin/bash

#SBATCH --time=01:00:00

### Usage
# sbatch mk_manifest_purple_output.sh -d <data_request_version> -m <metadata_folder_name> -p <purple_version>
# e.g. sbatch 08_make_hartwig_vcf_manifest/mk_manifest_purple_output.sh -d DR-104-update4 -m metadata_2021 -p purplesoft3.3 -r yes

# get commandline input
while getopts ":d:m:p:r:" opt; do
  case ${opt} in
    d )
		data_request_version=$OPTARG
      ;;
    m )
    metadata_folder=$OPTARG
      ;;
    p )
    purple_version=$OPTARG
      ;;
    r )
    postprocessed_vcf=$OPTARG
      ;;
    \? ) echo "Usage: cmd [-d] [-m] [-p] [-r]"
      ;;
		: ) echo "Invalid option: $OPTARG requires an argument" 1>&2
	    ;;
  esac
done

# define dirs
base_dir=$(dirname $(dirname $PWD))

Rscript $base_dir/code/08_make_hartwig_vcf_manifest/mk_manifest_purple_output.R \
$data_request_version $metadata_folder $purple_version $postprocessed_vcf

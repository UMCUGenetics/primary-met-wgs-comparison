#!/bin/bash

### Usage
# sbatch mk_linx_driver_manifest.sh -d DR-104-update4 -l linxsoft1.17

# get the linx version as commandline input to pass down to Rscript
while getopts ":d:l:" opt; do
  case ${opt} in
    d )
    data_request_version=$OPTARG
      ;;
    l )
    linx_version=$OPTARG
      ;;
    \? ) echo "Usage: cmd [-d] [-l]"
      ;;
    : ) echo "Invalid option: $OPTARG requires an argument" 1>&2
	    ;;
  esac
done

Rscript /path/to/Rscript/mk_linx_driver_manifest.R \
$data_request_version $linx_version

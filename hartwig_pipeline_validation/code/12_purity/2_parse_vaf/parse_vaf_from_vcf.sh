#!/bin/bash

#SBATCH --time=02:00:00

### Description
# this script parses the raw VAF values from VCF files of 100% pure tumor samples.
# requires SnpSift version 5.1d (build 2022-04-19 15:50), by Pablo Cingolani

### Input
# this script needs a manifest file specified on the command line (generated in step 12.1)

### Output
# A Table containing the sample_id, chromosome, position, gene name, subclnal likelihood and VAF of each PASS variant in that sample.

### Usage
# sbatch parse_driver_clonality_from_vcf.sh -m <manifest_file>
# e.g. sbatch parse_vaf_from_vcf.sh -m purple_output_paths_postprocessed_pure.txt.gz

# get the commandline input
while getopts ":m:" opt; do
  case ${opt} in
    m )
    manifest_file=$OPTARG
      ;;
    \? ) echo "Usage: cmd [-m]"
      ;;
    : ) echo "Invalid option: $OPTARG requires an argument" 1>&2
	    ;;
  esac
done

# get current date as variable
current_date=$(date +"%Y_%m_%d")

# define dirs
base_dir=$(dirname $(dirname $(dirname $PWD)))

manifest_path=$base_dir/data/processed/metadata/manifest/$manifest_file
job_dir=$base_dir/results/12_purity/2_parse_vaf/${current_date}/scripts/jobs/; mkdir -p $job_dir; cd $job_dir

output_dir=$base_dir/results/12_purity/2_parse_vaf/${current_date}/results/; mkdir -p $output_dir

# get header line of manifest file
manifest_column_names=$(zcat $manifest_path | head -n 1)

echo "the column names are: $manifest_column_names"

# init counter
counter=0

# read the manifest line by line, create one job script for one line and submit that job if it is not .done yet
zcat $manifest_path | tail -n +2 | while read $manifest_column_names; do
	let counter=counter+1

	job_file=$job_dir/pvfv_${sample}.job

	echo "#!/bin/bash
	#SBATCH --job-name=$(basename $job_file)
	#SBATCH --time=00:20:00
	#SBATCH --mem=10G

	zcat ${purple_vcf} | java -jar /path/to/snpEff/SnpSift.jar filter \"( FILTER = 'PASS' )\" | \
  java -jar /path/to/snpEff/SnpSift.jar extractFields - CHROM POS ANN[0].GENE SUBCL GEN[1].AF | \
  awk -F'\t' 'BEGIN {OFS = FS} NR==1{print \"sample_id\",\"chrom\",\"pos\",\"gene\",\"subclonal_likelihood\",\"vaf\";next}{print \"${sample}\",\$1,\$2,\$3,\$4,\$5}' - > $output_dir/${sample}_vaf.temp && touch ${job_file}.done
	" > $job_file

	if [[ ! -f ${job_file}.done ]]; then
		sbatch -o $(basename $job_file).out --job-name=parse_vaf_from_vcf $job_file
	fi

	sleep 0.5s

	# if [[ $counter -eq 2 ]]; then break; fi
done

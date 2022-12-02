#!/bin/bash

#SBATCH --time=02:00:00

### Usage
# sbatch extract_sv_sigs_length_cutoff.sh -m <manifest_file> -l <sv_length_cutoff>

# get the length_cutoff as commandline input
while getopts ":m:l:" opt; do
  case ${opt} in
    m )
    manifest_file=$OPTARG
      ;;
    l )
		length_cutoff=$OPTARG
      ;;
    \? ) echo "Usage: cmd [-m] [-l]"
      ;;
		: ) echo "Invalid option: $OPTARG requires an argument" 1>&2
	    ;;
  esac
done

# get current date as variable
current_date=$(date +"%Y_%m_%d")

# define dirs
base_dir=$(dirname $(dirname $PWD))

manifest_path=$base_dir/data/processed/metadata/manifest/$manifest_file
job_dir=$base_dir/results/2_extract_mutational_signatures_from_vcfs/${current_date}/scripts/sv_jobs/; mkdir -p $job_dir

# create output dirs if they dont exist
output_dir=$base_dir/results/2_extract_mutational_signatures_from_vcfs/${current_date}/results/matrices; mkdir -p $output_dir; cd $output_dir
mkdir sv_hmf_len_cutoff sv_pcawg_len_cutoff

# Rscript path
R_SCRIPT=$sv_dir/code/2_extract_mutational_signatures_from_vcfs/extract_sv_sigs_length_cutoff.R

# init counter
counter=0

# get header line of manifest file
manifest_column_names=$(zcat $manifest_path | head -n 1)

# read the manifest line by line, create one job script for one line and submit that job if it is not .done yet
zcat $manifest_path | tail -n +2 | while read $manifest_column_names; do
	let counter=counter+1

	job_file=$job_dir/svcutoff_${sample}.job

	echo "#!/bin/bash
	#SBATCH --job-name=$(basename $job_file)
	#SBATCH --output=$(basename $job_file).o
	#SBATCH --time=00:20:00
	#SBATCH --mem=10G

	Rscript $R_SCRIPT $sample $sv_hmf $sv_pcawg $length_cutoff $output_dir && touch ${job_file}.done
	" > $job_file

	if [[ ! -f ${job_file}.done ]]; then
		sbatch $job_file
	fi

	# if [[ $counter -eq 2 ]]; then break; fi
done

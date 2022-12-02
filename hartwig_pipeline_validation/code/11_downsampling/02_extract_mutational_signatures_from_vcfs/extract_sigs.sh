#!/bin/bash

#SBATCH --time=02:00:00

### Usage
# sbatch extract_sigs.sh

# get current date as variable
current_date=$(date +"%Y_%m_%d")

# define dirs
base_dir=$(dirname $(dirname $(dirname $PWD)))

manifest_path=$base_dir/data/processed/metadata/manifest/downsampled_vcf_paths.txt.gz
job_dir=$base_dir/results/11_downsampling/02_extract_mutational_signatures_from_vcfs/${current_date}/scripts/jobs/; mkdir -p $job_dir

# create output dirs if they dont exist
output_dir=$base_dir/results/11_downsampling/02_extract_mutational_signatures_from_vcfs/${current_date}/results/matrices; mkdir -p $output_dir

cd $job_dir

# Rscript path
R_SCRIPT=$base_dir/code/11_downsampling/02_extract_mutational_signatures_from_vcfs/extract_sigs.R

# initate counter
counter=0

# get header line of manifest file
manifest_column_names=$(zcat $manifest_path | head -n 1)

echo "the column names are: $manifest_column_names"

# read the manifest line by line, create one job script for one line and submit that job if it is not .done yet
zcat $manifest_path | tail -n +2 | while read $manifest_column_names; do
	let counter=counter+1

	job_file=$job_dir/extract_smnv_${sample}.job

	echo "#!/bin/bash
	#SBATCH --job-name=$(basename $job_file)
	#SBATCH --time=00:30:00
	#SBATCH --mem=8G

	Rscript $R_SCRIPT $sample $som_vcf $output_dir && touch ${job_file}.done
	" > $job_file

	if [[ ! -f ${job_file}.done ]]; then
		echo -n "[$counter] Submitting: $sample; "
		sbatch $job_file --job-name=$(basename $job_file)
	else
		echo "[$counter] Skipping: $sample"
	fi

	#if [[ $counter -eq 1 ]]; then break; fi
done

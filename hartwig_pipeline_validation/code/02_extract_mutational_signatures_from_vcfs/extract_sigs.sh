#!/bin/bash

#SBATCH --time=02:00:00

### Usage
# sbatch extract_sigs.sh -m <manifest_file> -l <sv_langth_cutoff> (e.g. sbatch extract_sigs.sh -l 500)
# e.g. sbatch extract_sigs.sh -m vcf_paths.txt.gz -l 500
# for <sv_langth_cutoff> 500 is recommended, because that is the minimum SV length that the PCAWG pipeline calls

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
job_dir=$base_dir/results/02_extract_mutational_signatures_from_vcfs/${current_date}/scripts/smnv_jobs/; mkdir -p $job_dir

# create output dirs if they dont exist
output_dir=$base_dir/results/02_extract_mutational_signatures_from_vcfs/${current_date}/results/matrices; mkdir -p $output_dir; cd $output_dir
mkdir snv_hartwig snv_pcawg dbs_hartwig dbs_pcawg indel_hartwig indel_pcawg sv_hartwig sv_pcawg

cd $job_dir

# Rscript path
R_SCRIPT=$base_dir/code/02_extract_mutational_signatures_from_vcfs/extract_sigs.R

# initate counter
counter=0

# get header line of manifest file
manifest_column_names=$(zcat $manifest_path | head -n 1)

echo "the column names are: $manifest_column_names"

# read the manifest line by line, create one job script for one line and submit that job if it is not .done yet
zcat $manifest_path | tail -n +2 | while read $manifest_column_names; do
	let counter=counter+1

	job_file=$job_dir/extract_${sample}.job

	echo "#!/bin/bash
	#SBATCH --job-name=$(basename $job_file)
	#SBATCH --time=00:20:00
	#SBATCH --mem=10G

	Rscript $R_SCRIPT $sample $som_hartwig $som_pcawg_snv $som_pcawg_indel $sv_hartwig $sv_pcawg $output_dir && touch ${job_file}.done
	" > $job_file

	if [[ ! -f ${job_file}.done ]]; then
		sbatch -o $(basename $job_file).out $job_file
	fi

	# if [[ $counter -eq 2 ]]; then break; fi
done

# execute the SV context profile extraction with different length cutoff
sbatch $base_dir/code/02_extract_mutational_signatures_from_vcfs/extract_sv_sigs_length_cutoff.sh -l $length_cutoff

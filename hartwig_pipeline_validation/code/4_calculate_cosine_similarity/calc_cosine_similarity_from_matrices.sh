#!/bin/bash

#SBATCH --time=02:00:00

### Usage
# sbatch calc_cosine_similarity_from_matrices.sh -m <manifest_file>

# get commandline input
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
base_dir=$(dirname $(dirname $PWD))

manifest_path=$base_dir/data/processed/metadata/manifest/$manifest_file
job_dir=$base_dir/results/4_calculate_cosine_similarity/${current_date}/scripts/jobs/; mkdir -p $job_dir; cd $job_dir

output_dir=$base_dir/results/4_calculate_cosine_similarity/${current_date}/results/; mkdir -p $output_dir

# path to R_SCRIPT
R_SCRIPT=$base_dir/code/4_calculate_cosine_similarity/calc_cosine_similarity_from_matrices.R

# create all the output files with headers
echo -e "sample_id\tsnv_cosine_similarity" > $output_dir/snv_cosine_similarity.tsv
echo -e "sample_id\tdbs_cosine_similarity" > $output_dir/dbs_cosine_similarity.tsv
echo -e "sample_id\tindel_cosine_similarity" > $output_dir/indel_cosine_similarity.tsv
echo -e "sample_id\tsv_cosine_similarity" > $output_dir/sv_cosine_similarity.tsv
echo -e "sample_id\tsv_cosine_similarity_len_cutoff" > $output_dir/sv_cosine_similarity_len_cutoff.tsv

# get header line of manifest file
manifest_column_names=$(zcat $manifest_path | head -n 1)

echo "the column names are: $manifest_column_names"

# initialize counter
counter=0

# read the manifest line by line, create one job script for one line and submit that job if it is not .done yet
zcat $manifest_path | tail -n +2 | while read $manifest_column_names; do
	let counter=counter+1

	job_file=$job_dir/cos_${sample}.job

	echo "#!/bin/bash
	#SBATCH --job-name=$(basename $job_file)
	#SBATCH --time=00:20:00
	#SBATCH --mem=10G

	Rscript $R_SCRIPT $sample $snv_hartwig $dbs_hartwig $indel_hartwig $sv_hartwig $sv_hartwig_len_cutoff $snv_pcawg $dbs_pcawg $indel_pcawg $sv_pcawg $sv_pcawg_len_cutoff $output_dir && touch ${job_file}.done
	" > $job_file

	if [[ ! -f ${job_file}.done ]]; then
		sbatch -o $(basename $job_file).out $job_file
	fi

	# if [[ $counter -eq 2 ]]; then break; fi
	sleep 0.5s

done

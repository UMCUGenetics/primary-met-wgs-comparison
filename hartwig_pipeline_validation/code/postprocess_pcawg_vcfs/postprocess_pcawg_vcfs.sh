#!/bin/bash

#SBATCH --time=04:00:00

### Usage
# sbatch postprocess_pcawg_vcfs.sh -m <manifest_file
# e.g. sbatch postprocess_pcawg_vcfs.sh -m purple_output_paths.txt.gz

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

manifest_path=$base_dir/data/processed/metadata/manifest/$manifest_path
job_dir=$base_dir/results/postprocess_pcawg_vcfs/scripts/jobs/; mkdir -p $job_dir; cd $job_dir

# path to the Rscript
R_SCRIPT=$base_dir/snv_mnv/analysis/postprocess_pcawg_vcfs/postprocess_pcawg_vcfs.R

# get header line of manifest file
manifest_column_names=$(zcat $manifest_path | head -n 1)

echo "the column names are: $manifest_column_names"
echo "using this manifest file: $(basename $manifest_path)"

# init counter
counter=0

# read the manifest line by line, create one job script for one line and submit that job if it is not .done yet
zcat $manifest_path | grep -E '^DO' | tail -n +1 | while read $manifest_column_names; do
	let counter=counter+1

	# define output dir (same as original VCF)
	output_dir=$(dirname $purple_vcf)

	# name the processed VCF
	zcat $purple_vcf | grep -B1000 '#CHROM' > $output_dir/${sample}T.purple.somatic.postprocessed.vcf.header
	vcf_body="${sample}T.purple.somatic.postprocessed.vcf.body"

	job_file=$job_dir/postprocess_${sample}.job

	echo "#!/bin/bash
	#SBATCH --job-name=$(basename $job_file)
	#SBATCH --time=01:00:00
	#SBATCH --mem=10G

	Rscript $R_SCRIPT $sample $vcf_body $purple_vcf $output_dir && touch ${job_file}.done
	" > $job_file

	if [[ ! -f ${job_file}.done ]]; then
		sbatch -o $(basename $job_file).out --time=01:00:00 --mem=16G --job-name=postprocess_pcawg_vcfs_${sample} $job_file
	fi

	# combine header and body
	sbatch -o $(basename $job_file).combine.out --dependency=singleton --job-name=postprocess_pcawg_vcfs_${sample}\
	--wrap="gzip -c $output_dir/${sample}T.purple.somatic.postprocessed.vcf.header $output_dir/${sample}T.purple.somatic.postprocessed.vcf.body > $output_dir/${sample}T.purple.somatic.postprocessed.vcf.gz
	rm $output_dir/${sample}T.purple.somatic.postprocessed.vcf.header $output_dir/${sample}T.purple.somatic.postprocessed.vcf.body"

	sleep 0.5s

	# if [[ $counter -eq 1 ]]; then break; fi
done

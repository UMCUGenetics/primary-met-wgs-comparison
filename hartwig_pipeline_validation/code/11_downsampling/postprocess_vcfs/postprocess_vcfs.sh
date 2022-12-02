#!/bin/bash

#SBATCH --time=04:00:00

### Usage
# sbatch postprocess_vcfs.sh

# get current date as variable
current_date=$(date +"%Y_%m_%d")

# define dirs
base_dir=$(dirname $(dirname $(dirname $PWD)))

manifest_path=$base_dir/data/processed/metadata/manifest/downsampled_vcf_paths.txt.gz
job_dir=$base_dir/results/11_downsampling/postprocess_vcfs/scripts/jobs/; mkdir -p $job_dir; cd $job_dir

# path to the Rscript
R_SCRIPT=$base_dir/code/11_downsampling/postprocess_vcfs/postprocess_vcfs.R

# get header line of manifest file
manifest_column_names=$(zcat $manifest_path | head -n 1)

echo "the column names are: $manifest_column_names"
echo "using this manifest file: $(basename $manifest_path)"

# init counter
counter=0

# read the manifest line by line, create one job script for one line and submit that job if it is not .done yet
zcat $manifest_path | tail -n +2 | while read $manifest_column_names; do
	let counter=counter+1

	# define output dir
	output_dir=$(dirname $som_vcf)

	# name the processed VCF
	zcat $som_vcf | grep -B1000 '#CHROM' > $output_dir/${sample}T.purple.somatic.postprocessed.vcf.header
	vcf_body="${sample}T.purple.somatic.postprocessed.vcf.body"

	job_file=$job_dir/postprocess_${sample}.job

	echo "#!/bin/bash
	#SBATCH --job-name=$(basename $job_file)
	#SBATCH --time=01:00:00
	#SBATCH --mem=10G

	Rscript $R_SCRIPT $sample $vcf_body $som_vcf $output_dir && touch ${job_file}.done
	" > $job_file

	if [[ ! -f ${job_file}.done ]]; then
		sbatch -o $(basename $job_file).out --time=01:00:00 --mem=16G --job-name=postprocess_vcfs_${sample} $job_file

		# combine header and body
		sbatch -o $(basename $job_file).combine.out --dependency=singleton --job-name=postprocess_vcfs_${sample}\
		--wrap="rm $output_dir/${sample}T.purple.somatic.postprocessed.vcf.gz
		gzip -c $output_dir/${sample}T.purple.somatic.postprocessed.vcf.header $output_dir/${sample}T.purple.somatic.postprocessed.vcf.body > $output_dir/${sample}T.purple.somatic.postprocessed.vcf.gz
		rm $output_dir/${sample}T.purple.somatic.postprocessed.vcf.header $output_dir/${sample}T.purple.somatic.postprocessed.vcf.body"
	fi

	sleep 0.5s

	# if [[ $counter -eq 2 ]]; then break; fi
done

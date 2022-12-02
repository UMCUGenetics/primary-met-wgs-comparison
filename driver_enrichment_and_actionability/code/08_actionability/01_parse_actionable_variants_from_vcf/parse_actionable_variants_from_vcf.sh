#!/bin/bash

#SBATCH --time=02:00:00

### Usage
# sbatch parse_actionability_from_vcf.sh -m <manifest_file>
# e.g. sbatch parse_actionability_from_vcf.sh -m purple_output_paths_postprocessed.txt.gz

### input
# this script needs a manifest file created in step 1

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
job_dir=$base_dir/results/08_actionability/${current_date}/scripts/jobs/; mkdir -p $job_dir; cd $job_dir

output_dir=$base_dir/results/08_actionability/${current_date}/results/; mkdir -p $output_dir

# path to the Rscript
R_SCRIPT=$base_dir/code/08_actionability/parse_actionable_variants_from_vcf.R

# get header line of manifest file
manifest_column_names=$(zcat $manifest_path | head -n 1)

echo "the column names are: $manifest_column_names"

# init counter
counter=0

# read the manifest line by line, create one job script for one line and submit that job if it is not .done yet
# sample purple_vcf mutationaltimer purity
zcat $manifest_path | tail -n +2 | while read $manifest_column_names; do
	let counter=counter+1

	job_file=$job_dir/pac_${sample}.job

	echo "#!/bin/bash

	#SBATCH --job-name=$(basename $job_file)
	#SBATCH --time=00:40:00
	#SBATCH --mem=32G

	Rscript $R_SCRIPT $sample $purple_vcf $output_dir && touch ${job_file}.done
	" > $job_file

	if [[ ! -f ${job_file}.done ]]; then
		sbatch -o $(basename $job_file).out --job-name=parse_actionability --mem=64G --time=01:00:00 $job_file
	fi

	sleep 0.5s

	# if [[ $counter -eq 2 ]]; then break; fi
done

### combine all the separate table into one big table
cd $output_dir

# create the file
echo -e "sample_id\tgene\tvar_type\tintron_exon_rank\tprotein_residue_change" \
> $output_dir/actionable_variants_prefilter.tsv

# create combine_tables.sh, a Bash script that combines all the short *.temp tables into a big table
echo '#!/bin/bash

#SBATCH --job-name=combine_tables
#SBATCH --output=combine_tables.o
#SBATCH --time=02:00:00
#SBATCH --mem=8G

echo "combining tables..."

for file in *.temp; do
	cat $file >> actionable_variants_prefilter.tsv
	rm $file
done

gzip actionable_variants_prefilter.tsv

echo "...done"
' > combine_tables.sh

sbatch --dependency=singleton --job-name=parse_actionability combine_tables.sh

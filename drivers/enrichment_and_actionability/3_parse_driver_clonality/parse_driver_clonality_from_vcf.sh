#!/bin/bash

#SBATCH --time=02:00:00

### Usage
# sbatch parse_driver_clonality_from_vcf.sh -l <linx_driver_catalog_date>
# e.g. sbatch parse_driver_clonality_from_vcf.sh -m 2022_02_16_purple_output_paths_postprocessed.txt.gz -l 2021_12_03

### input
# this script needs a manifest file created by /P0022_HMF_validation/snv_mnv/analysis/9_make_hmf_vcf_manifest/mk_manifest_purple_output.sh

# get the commandline input
while getopts ":l:" opt; do
  case ${opt} in
    m )
    manifest_file=$OPTARG
      ;;
    l )
    linx_catalog_date=$OPTARG
      ;;
    \? ) echo "Usage: cmd [-l]"
      ;;
    : ) echo "Invalid option: $OPTARG requires an argument" 1>&2
	    ;;
  esac
done

# get current date as variable
current_date=$(date +"%Y_%m_%d")

# define dirs
base_dir=/base/path

manifest_path=$base_dir/path/to/$manifest_file
job_dir=$base_dir/path/to/job/scripts/${current_date}_scripts/jobs/; mkdir -p $job_dir; cd $job_dir

output_dir=$base_dir/path/to/output/${current_date}/; mkdir -p $output_dir

# path to the Rscript
R_SCRIPT=$base_dir/path/to/Rscript/parse_driver_clonality_from_vcf.R

# get header line of manifest file
manifest_column_names=$(zcat $manifest_path | head -n 1)

echo "the column names are: $manifest_column_names"

# init counter
counter=0

# read the manifest line by line, create one job script for one line and submit that job if it is not .done yet
# sample purple_vcf clonality
zcat $manifest_path | tail -n +2 | while read $manifest_column_names; do
	let counter=counter+1

	job_file=$job_dir/pdc_${sample}.job

	echo "#!/bin/bash
	#SBATCH --job-name=$(basename $job_file)
	#SBATCH --time=00:20:00
	#SBATCH --mem=10G

	Rscript $R_SCRIPT $sample $purple_vcf $output_dir $linx_catalog_date && touch ${job_file}.done
	" > $job_file

	if [[ ! -f ${job_file}.done ]]; then
		sbatch -o $(basename $job_file).out --job-name=parse_driver_clonality $job_file
	fi

	sleep 0.5s

	# if [[ $counter -eq 2 ]]; then break; fi
done

### combine all the separate table into one big table
cd $output_dir

# create the file + header for the combined SV length cutoff
echo -e "sample_id\tgene\tsubclonal_likelihood\tclonality\ttier" \
> $output_dir/driver_clonality.tsv

# create combine_tables.sh, a Bash script that combines all the short *.temp tables into a big table
echo '#!/bin/bash
#SBATCH --job-name=combine_tables
#SBATCH --output=combine_tables.o
#SBATCH --time=02:00:00
#SBATCH --mem=8G

echo "combining tables..."

for file in *.temp; do
	cat $file >> driver_clonality.tsv
	rm $file
done

echo "...done"
' > combine_tables.sh

sbatch --dependency=singleton --job-name=parse_driver_clonality combine_tables.sh

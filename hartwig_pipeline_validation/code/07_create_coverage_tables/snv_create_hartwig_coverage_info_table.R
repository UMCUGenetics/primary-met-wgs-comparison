############### 90X vs 60X/30X coverage comparison ######################################
############### Create the hartwig processed hartwig coverage table ###########################

# author: remy (sascha)
# date: 04/06/2021
# last updated: 04/06/2021

### Description
# This script fetches the tumor metadata information about hartwig generated samples
# specifically the sample_id and cohort column, then adds a 3rd and 4th column
# called 'coverage' and 'coverage_label' representing the average read coverage of that sample.
# All samples are marked as 90X since there is currently no coverage info available for hartwig samples.

### Input
# The metadata.tsv file
# PCAWG specimen table

### Output
# A table that contains the Hartwig sample_id, the mean_coverage of 90 and a category column 'coverage_label' of '90X',
# because all Hartwig samples are sequenced at 90X.

### Usage
# sbatch snv_create_hartwig_coverage_info_table.sh

# libraries
source(paste0(here::here(), '/code/r_objects/libs.R'))

#========= Path prefixes =========#
base_dir <- list(
  path=paste0(here::here())
)

for(i in base_dir){
  if(dir.exists(i)){
    base_dir <- i
    break
  }
}

# pass commandline directory arguments to 'args' vector
args = commandArgs(trailingOnly=TRUE)

# parse commandline input
output_folder <- args[1]

# ------------------------- METADATA
metadata <- read_tsv(paste0(base_dir, '/data/processed/metadata/metadata.tsv'))

# ------------------------- PCAWG specimen table
specimen <- read_tsv(paste0(base_dir, '/data/external/pcawg_specimen_histology_August2016_v9_complete.tsv'))

# select only cohort and sample_id col
metadata <- metadata %>%
  select(sample_id, cohort)

# keep only unique combinations
metadata <- metadata %>%
  distinct(sample_id, cohort)

# filter for hartwig data only
metadata_hartwig <- metadata %>%
  filter(cohort == 'Hartwig')

# add 'coverage' and 'coverage_label' column
metadata_hartwig <- metadata_hartwig %>%
  mutate(mean_coverage = 90,
         coverage_label = '90X')

# write to disk
write.table(metadata_hartwig, 
            file = paste0(output_folder, 'hartwig_coverage_info.tsv'),
            sep = '\t', quote = FALSE, append = FALSE, row.names = FALSE,
            col.names = TRUE)
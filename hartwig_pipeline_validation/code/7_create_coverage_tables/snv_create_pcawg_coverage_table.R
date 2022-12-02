############### 90X vs 60X/30X coverage comparison ######################################
############### Create the HMF processed PCAWG coverage table ###########################

# author: remy (sascha)
# date: 04/06/2021
# last updated: 22/04/2022

### Description
# This script gets the mean coverage from the bam metrics files of the HMF processed PCAWG
# samples, stores those values in a table with one row being one sample, then creates a second
# coverage_label column based in the mean coverage value
# mean coverage < 49 ~ '38X'
# mean coverage >= 49 ~ '60X'

### Input
# The manifest file you generated in step 7.

### Output
# A temporary table that holds the sample_id and the mean_coverage of that sample 
# as well as a factor column denoting whether the sample is 30X or 60X sequenced.

### Usage
# sbatch snv_create_pcawg_coverage_table.sh

# global options
options(stringsAsFactors=F)

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

# get argument input
args <- commandArgs(trailingOnly=TRUE)
sample_id <- args[1]
snv_pcawg_wgsmetrics <- args[2]
output_dir <- args[3]

###################################### PCAWG

# read in the Hartwig pipeline processed PCAWG sample, rename column names
pcawg_mean_coverage <- read_tsv(snv_pcawg_wgsmetrics,
                                  col_names = TRUE) %>%
  rename_with(tolower)

# add two columns: the sample name and the coverage_label
pcawg_mean_coverage <- pcawg_mean_coverage %>%
  mutate(sample_id = sample_id,
         cohort = 'PCAWG',
         coverage_label = case_when(mean_coverage < 49 ~ '30X',
                                    mean_coverage >= 49 ~ '60X')) %>%
  select(sample_id, cohort, everything())

# write this table to disk
write.table(pcawg_mean_coverage, 
            file = paste0(output_dir, '/', sample_id,  '_pcawg_coverage.temp'),
            sep = '\t', quote = FALSE, append = FALSE, row.names = FALSE,
            col.names = FALSE)
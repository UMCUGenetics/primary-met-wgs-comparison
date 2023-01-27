############### Clean the extracted variants ###############
# author: remy (sascha)
# date: 24/01/2023
# last updated: 24/01/2023

### Description
# This script loads the extracted variant table from a VCF file, cleans it and saves the cleaned up table to disk

### Input
# a (dirty) table of extracted variants from a VCF file.

### Output
# A (cleaned up) table of extracted variants with sample_id added as an additional column

# global options
options(stringsAsFactors = FALSE)

# command line input
## nextflow
args <- commandArgs(trailingOnly = TRUE)
dirty_variants_table_date <- args[1]
sample_id <- args[2]
sage_filtered <- as.logical(args[3])
output_dir <- args[4]
dirty_variants_table_input <- file("stdin")

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

# libraries
source(paste0(here::here(), '/code/r_objects/libs.R'))

# get today's date
current_date <- format(Sys.Date(), '%Y_%m_%d')

# --------------------------------- DIRTY VARIANTS TABLE

dirty_variants_table <- read_tsv(
  file = dirty_variants_table_input,
  na = 'NA',
  col_names = TRUE
  ) %>%
  rename_with(str_to_title)

# add sampleId to the table
clean_variants_table <- dirty_variants_table %>%
  mutate(SampleId = sample_id, .before = Chrom)

# write to standard out
cat(format_tsv(
  clean_variants_table
))
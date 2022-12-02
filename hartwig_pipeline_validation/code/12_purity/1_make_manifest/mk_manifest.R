############### Make purity manifest ############### 
# author: remy (sascha)
# date: 01/11/2022
# last updated: 01/11/2022

### Description
# makes a manifest file of samples with 100% purity

### Input
# - The manifest file of ALL samples you generated in step 8.
# - metadata.tsv

### Output
# A filtered manifest file that only contains 100% pure tumor samples.

# global options
options(stringsAsFactors = F)

# libs
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

# ---------------------- MANIFEST

input_file <- gzfile(paste0(base_dir, '/purple_output_paths_postprocessed.txt.gz'))
manifest <- read_tsv(
  file = input_file
)

# ---------------------- PURITY

pure_samples <- read_tsv(paste0(base_dir, '/data/processed/metadata/metadata.tsv'),
                     col_types = 'icccccccccccccccllllcccdiiildclccllllllliiiiiiiiiii') %>%
  # select relevant columns
  select(sample_id, cohort, cancer_type, tumor_purity) %>%
  filter(near(tumor_purity, 1))

# filter the full manifest file to only include samples that have a purity of 1
manifest <- manifest %>%
  semi_join(pure_samples, by = c('sample' = 'sample_id'))

# write gzipped table to disk
write.table(
  manifest,
  file = gzfile(paste0(base_dir, '/purple_output_paths_postprocessed_pure.txt.gz')),
  sep='\t',quote=F,row.names=F
)
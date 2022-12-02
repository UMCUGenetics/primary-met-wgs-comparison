############### Create the VCF manifest file ############
# author: remy (sascha)
# date: 22/03/2022
# last updated: 22/03/2022

### Description
# This script creates a manifest file of the downsampled hartwig samples (30X, 38X and 60X coverage) together with the original 100X VCFs.

### Input
# The metadata subset you created in step 15.0.

### Output
# A manifest file that contains the VCF files paths to the downsampled VCFs: downsampled_vcf_paths.txt.gz

### Usage
# sbatch mk_manifest.sh

# global options
options(stringsAsFactors=F)

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

# commandline input
args <- commandArgs(trailingOnly=TRUE)
postprocessed_vcf <- args[1]

#========= Main =========#
# load downsampled_samples, rename column names to lower case, filter the table for PCAWG cohort only
downsampled_samples <- read_tsv(paste0(base_dir, '/data/processed/downsampling/cov_downsampling_samples.tsv')) %>%
  rename_with(tolower) %>%
   mutate(cov_30x_id = paste0(sample_id, 'B'),
          cov_38x_id = paste0(sample_id, 'C'),
          cov_60x_id = paste0(sample_id, 'A')) %>%
   select(-sample_id) %>%
   pivot_longer(., cols = c(cov_30x_id:cov_60x_id), names_to = 'coverage', values_to = 'sample_id') %>%
   relocate(sample_id, .before = sample_id)

original_samples <- read_tsv(paste0(base_dir, '/data/processed/downsampling/cov_downsampling_samples.tsv')) %>%
   rename_with(tolower)

# hartwig processed pcawg vcfs base path
if (postprocessed_vcf == 'yes') {
  downsampled_data_dir <- paste0(base_dir, 'results/11_downsampling/postprocess_vcfs/results')
  original_data_dir <- paste0(base_dir, '/data/raw')
  
  # initiate list to stare individual VCF file paths
  downsampled_file_paths <- list()
  original_file_paths <- list()
  
  ## PURPLE SOMATIC VCF
  downsampled_file_paths$som_vcf <- with(downsampled_samples,{ paste0(downsampled_data_dir, '/', sample_id, 'T.purple.somatic.postprocessed.vcf.gz') })
  original_file_paths$som_vcf <- with(original_samples,{ paste0(original_data_dir, '/data_request/somatics/', sample_id,'/purple/', sample_id, '.purple.somatic.vcf.gz') }) 
  
  ## PURPLE SV VCF
  downsampled_file_paths$sv_vis_data <- with(downsampled_samples,{ paste0(downsampled_data_dir, '/', sample_id,'/', 'linx_downsampled/', sample_id, 'T.linx.vis_sv_data.tsv') })
  original_file_paths$sv_vis_data <- with(original_samples,{ paste0(original_data_dir, '/data_request/somatics/', sample_id,'/linx/', sample_id, '.linx.vis_sv_data.tsv') })
} else {
  downsampled_data_dir <- paste0(base_dir, '/data/raw')
  original_data_dir <- paste0(base_dir, '/data/raw')
  
  # initiate list to stare individual VCF file paths
  downsampled_file_paths <- list()
  original_file_paths <- list()
  
  ## PURPLE SOMATIC VCF
  downsampled_file_paths$som_vcf <- with(downsampled_samples,{ paste0(downsampled_data_dir, '/', sample_id,'/', 'purple_downsampled/', sample_id, 'T.purple.somatic.vcf.gz') })
  original_file_paths$som_vcf <- with(original_samples,{ paste0(original_data_dir, '/data_request/somatics/', sample_id,'/purple/', sample_id, '.purple.somatic.vcf.gz') }) 
  
  ## PURPLE SV VCF
  downsampled_file_paths$sv_vis_data <- with(downsampled_samples,{ paste0(downsampled_data_dir, '/', sample_id,'/', 'linx_downsampled/', sample_id, 'T.linx.vis_sv_data.tsv') })
  original_file_paths$sv_vis_data <- with(original_samples,{ paste0(original_data_dir, '/data_request/somatics/', sample_id,'/linx/', sample_id, '.linx.vis_sv_data.tsv') })
}

# transform list to data frame, give each row the corresponding sample_id
downsampled_file_paths <- as.data.frame(downsampled_file_paths)
rownames(downsampled_file_paths) <- downsampled_samples$sample_id

original_file_paths <- as.data.frame(original_file_paths)
rownames(original_file_paths) <- original_samples$sample_id

# combine downsampled and original manifest df
file_paths <- rbind(downsampled_file_paths, original_file_paths)

#========= Get existing files =========#
counter <- 0
pb <- txtProgressBar(max=nrow(downsampled_samples))
exist_files <- as.data.frame(lapply(file_paths, function(i){
   counter <<- counter+1
   setTxtProgressBar(pb, counter)
   file.exists(i)
}))
rownames(exist_files) <- c(downsampled_samples$sample_id,original_samples$sample_id)

selectExistingFile <- function(keys){
   idx <- apply(exist_files[,keys],1,which)

   file_paths_ss <- file_paths[,keys]
   counter <- 0
   out <- unname(apply(file_paths_ss,1,function(i){
      counter <<- counter+1
      i[ idx[counter] ]
   }))

   return(out)
}

# initate second list to store those file paths that actually exist on the HPC
file_paths_2 <- list()

# fill the list
file_paths_2$som_vcf <- file_paths$som_vcf
file_paths_2$sv_vis_data <- file_paths$sv_vis_data

# transform to data frame
file_paths_2 <- as.data.frame(file_paths_2)

#========= Export =========#
vcf_paths <- as.data.frame(lapply(file_paths_2, function(i){ gsub('/local/path','',i) }))
vcf_paths <- cbind(sample=c(downsampled_samples$sample_id, original_samples$sample_id), vcf_paths)

# arrange columns
vcf_paths <- vcf_paths %>%
   # filter(!str_detect(sample, pattern = regex('TI?I?$'))) %>% # optional filter for postprocessed samples
   arrange(sample)

output_file <- ifelse(
  postprocessed_vcf == 'yes',
  paste0(base_dir, '/data/processed/metadata/manifest/downsampled_postprocessed_vcf_paths.txt.gz'),
  paste0(base_dir, '/data/processed/metadata/manifest/downsampled_vcf_paths.txt.gz')
)

message(paste0('Generating file... ', output_file))

# write to disk
write.table(
   vcf_paths, gzfile(output_file),
   sep='\t', quote = F, row.names = F
)
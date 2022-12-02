############### Create the bam metrics manifest file ############
# author: remy (sascha)
# date: 19/05/2021
# last updated: 13/07/2021

### Description
# This script creates a manifest file of the purple bam wgsmetrics file that stores values like mean_coverage.
# The mean coverage is important for the coverage comparison plots.
# ! Bam metrics are only available for PCAWG samples, not for Hartwig (HMF) samples !

### Input
# metadata.tsv file

### Output
# A manifest file that contains the file paths to the BAM WGS metrics files that are part of the BAM output files 
# of the Hartwig pipeline.

### Usage
# sbatch mk_manifest_bam_metrics.sh

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

#========= Main =========#
# read metadata in, rename columns to lower then filter for the PCAWG cohort only
metadata <- read_tsv(paste0(base_dir, '/data/processed/metadata/metadata.tsv'),
                     col_types = 'icccccccccccccccllllcccdiiildclccllllllliiiiiiiiiii') %>%
  rename_with(tolower) %>%
  filter(cohort == "PCAWG")

# bam wgsmetrics files base path
bam_metrics_data_dir <- '/path/to/bam/metrics/'
if(dir.exists('/local/path')){
  bam_metrics_data_dir <- paste0('/local/path', bam_metrics_data_dir)
}

# intiate list where the paths to the individual wgsmetrics files will be stored
file_paths <- list()

# add the pcawg wgsmetrics file paths to the list using the metadata file
file_paths$snv_pcawg_wgsmetrics <- with(metadata,{ paste0(bam_metrics_data_dir, '/', sample_id, '/bam_metrics/', sample_id, 'T.wgsmetrics') })

# transform list to data frame and give it the sample names as rownames
file_paths <- as.data.frame(file_paths)
rownames(file_paths) <- metadata$sample_id

#========= Get existing files =========#
counter <- 0
pb <- txtProgressBar(max=nrow(metadata))
exist_files <- as.data.frame(lapply(file_paths, function(i){
   counter <<- counter+1
   setTxtProgressBar(pb, counter)
   file.exists(i)
}))
rownames(exist_files) <- metadata$sample_id

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

# initiate the second list which stores only files that exist
file_paths_2 <- list()

# add the pcawg wgsmetrics existing file paths to the list
file_paths_2$snv_pcawg_wgsmetrics <- file_paths$snv_pcawg_wgsmetrics

#========= Export =========#
bam_metrics_paths <- as.data.frame(lapply(file_paths_2, function(i){ gsub('/local/path','',i) }))
bam_metrics_paths <- cbind(sample=metadata$sample_id, bam_metrics_paths)

# write the manifest as gzipped file to the metadata/ folder
write.table(
  bam_metrics_paths, gzfile(paste0(base_dir, '/data/processed/metadata/manifest/bam_metrics_paths.txt.gz')),
   sep='\t',quote=F,row.names=F
)
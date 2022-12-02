############### Create the VCF manifest file ############
# author: remy (sascha)
# date: 07/05/2021
# last updated: 01/10/2021

### Description
# This script creates a manifest file of the Hartwig and PCAWG somatic VCF file paths to be used for
# subsequent analysis of the same.

### Input
# metadata.tsv file

### Output
# a manifest file that contains the file paths to the context matrices you generated in step 2.

### Usage
# sbatch mk_manifest_matrices.sh -f <matrices_folder>
# e.g. sbatch mk_manifest_matrices.sh -f 2021_12_03_matrices
# <matrices_folder> is the folder name where you wrote your matrices to.

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

# set todays date as the current_date
current_date <- format(Sys.Date(), '%Y_%m_%d')

# commandline argument parsing
args <- commandArgs(trailingOnly=TRUE)
matrices_folder <- args[1]

#========= Main =========#
# read in the matedata file to use it as an anchor for the generation of matrices file paths
metadata <- read_tsv(paste0(base_dir, "/data/processed/metadata/metadata.tsv")) %>%
  rename_with(tolower) %>%
  filter(cohort == "PCAWG")

# matrices base dir
matrices_data_dir <- paste0(base_dir, '/results/02_extract_mutational_signatures_from_vcfs/results')
if(dir.exists('/local/path')){
  matrices_data_dir <- paste0('/local/path', matrices_data_dir)
}

# combine path with the selected folder
matrices_data_dir <- paste0(matrices_data_dir, '/', matrices_folder, '/matrices')

# intiate list where the matrices paths will be stored
file_paths <- list()

# add hartwig matrices paths to the list
file_paths$snv_hartwig <- with(metadata,{ paste0(matrices_data_dir, '/snv_hartwig/', sample_id, '_snv_hartwig.txt') })
file_paths$dbs_hartwig <- with(metadata,{ paste0(matrices_data_dir, '/dbs_hartwig/', sample_id, '_dbs_hartwig.txt') })
file_paths$indel_hartwig <- with(metadata,{ paste0(matrices_data_dir, '/indel_hartwig/', sample_id, '_indel_hartwig.txt') })
file_paths$sv_hartwig <- with(metadata,{ paste0(matrices_data_dir, '/sv_hartwig/', sample_id, '_sv_hartwig.txt') })
file_paths$sv_hartwig_len_cutoff <- with(metadata,{ paste0(matrices_data_dir, '/sv_hartwig_len_cutoff/', sample_id, '_sv_hartwig_len_cutoff.txt') })
# add pcawg matrices paths to the list
file_paths$snv_pcawg <- with(metadata,{ paste0(matrices_data_dir, '/snv_pcawg/', sample_id, '_snv_pcawg.txt') })
file_paths$dbs_pcawg <- with(metadata,{ paste0(matrices_data_dir, '/dbs_pcawg/', sample_id, '_dbs_pcawg.txt') })
file_paths$indel_pcawg <- with(metadata,{ paste0(matrices_data_dir, '/indel_pcawg/', sample_id, '_indel_pcawg.txt') })
file_paths$sv_pcawg <- with(metadata,{ paste0(matrices_data_dir, '/sv_pcawg/', sample_id, '_sv_pcawg.txt') })
file_paths$sv_pcawg_len_cutoff <- with(metadata,{ paste0(matrices_data_dir, '/sv_pcawg_len_cutoff/', sample_id, '_sv_pcawg_len_cutoff.txt') })

# transform list to data frame and add rownames accordingly
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

# initiate second matrices file paths list that contain existing files
file_paths_2 <- list()

# add hartwig matrices paths to the list
file_paths_2$snv_hartwig <- file_paths$snv_hartwig
file_paths_2$dbs_hartwig <- file_paths$dbs_hartwig
file_paths_2$indel_hartwig <- file_paths$indel_hartwig
file_paths_2$sv_hartwig <- file_paths$sv_hartwig
file_paths_2$sv_hartwig_len_cutoff <- file_paths$sv_hartwig_len_cutoff
# add pcawg matrices paths to the list
file_paths_2$snv_pcawg <- file_paths$snv_pcawg
file_paths_2$dbs_pcawg <- file_paths$dbs_pcawg
file_paths_2$indel_pcawg <- file_paths$indel_pcawg
file_paths_2$sv_pcawg <- file_paths$sv_pcawg
file_paths_2$sv_pcawg_len_cutoff <- file_paths$sv_pcawg_len_cutoff
file_paths_2 <- as.data.frame(file_paths_2)

#========= Export =========#
matrices_paths <- as.data.frame(lapply(file_paths_2, function(i){ gsub('/local/path','',i) }))
matrices_paths <- cbind(sample=metadata$sample_id, matrices_paths)

# write to disk
write.table(
   matrices_paths, gzfile(paste0(base_dir,'/data/processed/metadata/manifest/', matrices_folder, '_matrices_paths.txt.gz')),
   sep='\t',quote=F,row.names=F
)
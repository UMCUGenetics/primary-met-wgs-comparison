############### Create the linx driver catalogue manifest file ############
# author: remy (sascha)
# date: 19/07/2021
# last updated: 12/04/2022

### Description
# This script creates a manifest file of the linx driver output file paths.

# global options
options(stringsAsFactors=F)

# libs
library(readr)
library(dplyr)
library(stringr)

#========= Path prefixes =========#
base_dir <- list(
  hpc='/base/path/to/your/files',
  mnt='/base/path/to/your/files',
  umc='/base/path/to/your/files'
)

for(i in base_dir){
   if(dir.exists(i)){
      base_dir <- i
      break
   }
}

hmf_metadata_dir <- list(
  hpc='/base/path/to/hmf/metadata',
  mnt='/base/path/to/hmf/metadata',
  umc='/base/path/to/hmf/metadata'
)

for(i in hmf_metadata_dir){
  if(dir.exists(i)){
    hmf_metadata_dir <- i
    break
  }
}

# set todays date as the output_date
output_date <- format(Sys.Date(), '%Y_%m_%d')

# command line input
args <- commandArgs(trailingOnly = TRUE)
data_request_version <- args[1]
linx_version <- args[2]

# local
# data_request_version <- 'DR-104-update4'
# linx_version <- 'linxsoft1.16'

#========= Main =========#
# read the metadata in, make all the column names lower case
metadata <- read_tsv(paste0(base_dir, '/subpath/to/metadata.tsv'), quote = '#') %>%
  rename_with(tolower)

# read the hmf metadata in (variable setName is important)
hmf_metadata <- read_tsv(paste0(hmf_metadata_dir, '/metadata.tsv')) %>%
  select(hmfSampleId, setName)

# join setName from the hmf_metadata to the metadata table
metadata <- metadata %>%
  left_join(hmf_metadata, by = c('sample_id_2' = 'hmfSampleId')) %>%
  relocate(setName, .after = sample_id_2)

# split by cohort
metadata_hmf <- metadata %>%
  filter(cohort == 'HMF')
metadata_pcawg <- metadata %>%
  filter(cohort == 'PCAWG')

hmf_data_dir <- '/path/to/hartwig/driver/catalogs/'
if(dir.exists('/mount/path')){
  hmf_data_dir <- paste0('/mount/path',
                         hmf_data_dir)
} else if (dir.exists('/different/mount/path')) {
  hmf_data_dir <- paste0('/different/mount/path',
                         hmf_data_dir)
}

pcawg_hmf_data_dir <- '/path/to/pcawg/driver/catalogs/'
if(dir.exists('/mount/path')){
  pcawg_hmf_data_dir <- paste0('/mount/path',
                               pcawg_hmf_data_dir)
} else if (dir.exists('/different/mount/path')) {
  pcawg_hmf_data_dir <- paste0('/different/mount/path',
                               pcawg_hmf_data_dir)
}

# initiate list to store the file paths in
file_paths <- list()

# add hmf linx driver catalog file path to the list
file_paths$linx_driver_hmf <- with(metadata_hmf,{ paste0(hmf_data_dir, data_request_version, '/somatics/', 
                                                         if (data_request_version == 'DR-104-update3') {
                                                           paste0(setName, '/', linx_version)
                                                         } else {
                                                           paste0(sample_id, '/linx')
                                                         },
                                                         '/', sample_id, '.linx.driver.catalog.tsv') })
# add hmf linx driver fusion file path to the list
file_paths$linx_fusion_hmf <- with(metadata_hmf,{ paste0(hmf_data_dir, data_request_version, '/somatics/',
                                                         if (data_request_version == 'DR-104-update3') {
                                                           paste0(setName, '/', linx_version)
                                                         } else {
                                                           paste0(sample_id, '/linx')
                                                         },
                                                         '/', sample_id, '.linx.fusion.tsv') })
# add hmf linx driver breakpoints file path to the list
file_paths$linx_breakpoints_hmf <- with(metadata_hmf,{ paste0(hmf_data_dir, data_request_version, '/somatics/',
                                                              if (data_request_version == 'DR-104-update3') {
                                                                paste0(setName, '/', linx_version)
                                                              } else {
                                                                paste0(sample_id, '/linx')
                                                              },
                                                              '/', sample_id, '.linx.vis_sv_data.tsv') })

# add pcawg linx driver catalog file path to the list
file_paths$linx_driver_pcawg <- with(metadata_pcawg,{ paste0(pcawg_hmf_data_dir, patient_id, '-from-jar', '/', linx_version, '/', sample_id, 'T.linx.driver.catalog.tsv') })
# add pcawg linx driver fusion file path to the list
file_paths$linx_fusion_pcawg <- with(metadata_pcawg,{ paste0(pcawg_hmf_data_dir, patient_id, '-from-jar', '/', linx_version, '/', sample_id, 'T.linx.fusion.tsv') })
# add pcawg linx driver breakpoints file path to the list
file_paths$linx_breakpoints_pcawg <- with(metadata_pcawg,{ paste0(pcawg_hmf_data_dir, patient_id, '-from-jar', '/', linx_version, '/', sample_id, 'T.linx.vis_sv_data.tsv') })

# transform list into a dataframe and change the column names
# driver catalog
linx_driver_hmf <- as.data.frame(file_paths$linx_driver_hmf) %>%
  rename('linx_driver_catalog' = `file_paths$linx_driver_hmf`)
linx_driver_pcawg <- as.data.frame(file_paths$linx_driver_pcawg) %>%
  rename('linx_driver_catalog' = `file_paths$linx_driver_pcawg`)

# fusions
linx_fusion_hmf <- as.data.frame(file_paths$linx_fusion_hmf) %>%
  rename('linx_fusion_catalog' = `file_paths$linx_fusion_hmf`)
linx_fusion_pcawg <- as.data.frame(file_paths$linx_fusion_pcawg) %>%
  rename('linx_fusion_catalog' = `file_paths$linx_fusion_pcawg`)

# breakpoints
linx_breakpoints_hmf <- as.data.frame(file_paths$linx_breakpoints_hmf) %>%
  rename('linx_breakpoints_catalog' = `file_paths$linx_breakpoints_hmf`)
linx_breakpoints_pcawg <- as.data.frame(file_paths$linx_breakpoints_pcawg) %>%
  rename('linx_breakpoints_catalog' = `file_paths$linx_breakpoints_pcawg`)

# bind columns, add to the list
file_paths <- rbind(linx_driver_hmf, linx_driver_pcawg)
file_paths[2] <- rbind(linx_fusion_hmf, linx_fusion_pcawg)
file_paths[3] <- rbind(linx_breakpoints_hmf, linx_breakpoints_pcawg)

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

# initate second list to store those file paths that actually exist on the HPC
file_paths_2 <- list()

# fill the list
file_paths_2$linx_driver_catalog <- file_paths$linx_driver_catalog
file_paths_2$linx_fusion_catalog <- file_paths$linx_fusion_catalog
file_paths_2$linx_breakpoints_catalog <- file_paths$linx_breakpoints_catalog
# transform to data frame
file_paths_2 <- as.data.frame(file_paths_2)

#========= Export =========#
if (dir.exists('/mount/path')) {
  export_file_paths <- as.data.frame(lapply(file_paths_2, function(i){ gsub('/mount/path','',i) })) 
} else {
  export_file_paths <- as.data.frame(lapply(file_paths_2, function(i){ gsub('/different/mount/path','',i) }))
}
export_file_paths <- cbind(sample=metadata$sample_id, export_file_paths)

write.table(
  export_file_paths, gzfile(paste0(base_dir,'/path/to/manifest/output/', output_date, '_linx_output_paths.txt.gz')),
   sep='\t',quote=F,row.names=F
)


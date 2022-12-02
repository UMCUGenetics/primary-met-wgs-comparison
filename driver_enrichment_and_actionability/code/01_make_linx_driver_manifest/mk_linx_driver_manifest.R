############### Create the linx driver catalogue manifest file ############
# author: remy (sascha)
# date: 19/07/2021
# last updated: 12/04/2022

### Description
# This script creates a manifest file of the linx driver output file paths.

### Input 
# metadata.tsv
# hartwig metadata file (you get this file with your data request)

### Output
# a manifest file that contains the sample_id and paths to LINX output files

# global options
options(stringsAsFactors=F)

# libs
here::set_here(path='../../')
source(paste0(here::here(), '/code/r_objects/libs.R'))

#========= Path prefixes =========#
base_dir <- list(
  path=paste0(here::here(), '/data')
)

for(i in base_dir){
   if(dir.exists(i)){
      base_dir <- i
      break
   }
}

hartwig_metadata_dir <- list(
  path=paste0(here::here(), '/data')
)

for(i in hartwig_metadata_dir){
  if(dir.exists(i)){
    hartwig_metadata_dir <- i
    break
  }
}

# set todays date as the current_date
current_date <- format(Sys.Date(), '%Y_%m_%d')

# command line input for newer linx versions
args <- commandArgs(trailingOnly = TRUE)
data_request_version <- args[1]
linx_version <- args[2]

#========= Main =========#
# read the metadata in, make all the column names lower case
metadata <- read_tsv(paste0(base_dir, '/processed/metadata/metadata.tsv'),
                     col_types = 'icccccccccccccccllllcccdiiildclccllllllliiiiiiiiiii') %>%
  rename_with(tolower)

# read the hmf metadata in (variable setName is important); you will get this metadata with your data request from Hartwig Medical Foundation
hartwig_metadata <- read_tsv(paste0(hartwig_metadata_dir, '/processed/metadata/hartwig_sample_metadata.tsv')) %>%
  select(hmfSampleId, setName)

# join setName from the hartwig_metadata to the metadata table
metadata <- metadata %>%
  left_join(hartwig_metadata, by = c('sample_id_2' = 'hmfSampleId')) %>%
  relocate(setName, .after = sample_id_2)

# split by cohort
metadata_hartwig <- metadata %>%
  filter(cohort == 'Hartwig')
metadata_pcawg <- metadata %>%
  filter(cohort == 'PCAWG')

#========= Path prefixes to your hartwig and PCAWG samples =========#
hartwig_data_dir <- '/hartwig/samples/path/' ### CHANGE base path to your Hartwig samples here
if(dir.exists('/local/path')){
  hartwig_data_dir <- paste0('/local/path',
                             hartwig_data_dir)
} else if (dir.exists('/other/local/path')) {
  hartwig_data_dir <- paste0('/other/local/path',
                             hartwig_data_dir)
}

pcawg_data_dir <- '/pcawg/sample/path/' ### CHANGE base path to your PCAWG samples here
if(dir.exists('/local/path')){
  pcawg_data_dir <- paste0('/local/path',
                                   pcawg_data_dir)
} else if (dir.exists('/other/local/path')) {
  pcawg_data_dir <- paste0('/other/local/path',
                                   pcawg_data_dir)
}

# initiate list to store the file paths in
file_paths <- list()

# add hartwig linx driver catalog file path to the list
file_paths$linx_driver_hartwig <- with(metadata_hartwig,{ paste0(hartwig_data_dir, data_request_version, '/somatics/', 
                                                                 sample_id, '/linx/', sample_id, '.linx.driver.catalog.tsv') })
# add hartwig linx driver fusion file path to the list
file_paths$linx_fusion_hartwig <- with(metadata_hartwig,{ paste0(hartwig_data_dir, data_request_version, '/somatics/', 
                                                                 sample_id, '/linx', sample_id, '.linx.fusion.tsv') })
# add hartwig linx driver breakpoints file path to the list
file_paths$linx_breakpoints_hartwig <- with(metadata_hartwig,{ paste0(hartwig_data_dir, data_request_version, '/somatics/', 
                                                                      sample_id, '/linx', sample_id, '.linx.vis_sv_data.tsv') })

# add pcawg linx driver catalog file path to the list
file_paths$linx_driver_pcawg <- with(metadata_pcawg,{ paste0(pcawg_data_dir, patient_id, '-from-jar', '/', linx_version, '/', sample_id, 'T.linx.driver.catalog.tsv') })
# add pcawg linx driver fusion file path to the list
file_paths$linx_fusion_pcawg <- with(metadata_pcawg,{ paste0(pcawg_data_dir, patient_id, '-from-jar', '/', linx_version, '/', sample_id, 'T.linx.fusion.tsv') })
# add pcawg linx driver breakpoints file path to the list
file_paths$linx_breakpoints_pcawg <- with(metadata_pcawg,{ paste0(pcawg_data_dir, patient_id, '-from-jar', '/', linx_version, '/', sample_id, 'T.linx.vis_sv_data.tsv') })

# transform list into a dataframe and change the column names
# driver catalog
linx_driver_hartwig <- as.data.frame(file_paths$linx_driver_hartwig) %>%
  rename('linx_driver_catalog' = `file_paths$linx_driver_hartwig`)
linx_driver_pcawg <- as.data.frame(file_paths$linx_driver_pcawg) %>%
  rename('linx_driver_catalog' = `file_paths$linx_driver_pcawg`)

# fusions
linx_fusion_hartwig <- as.data.frame(file_paths$linx_fusion_hartwig) %>%
  rename('linx_fusion_catalog' = `file_paths$linx_fusion_hartwig`)
linx_fusion_pcawg <- as.data.frame(file_paths$linx_fusion_pcawg) %>%
  rename('linx_fusion_catalog' = `file_paths$linx_fusion_pcawg`)

# breakpoints
linx_breakpoints_hartwig <- as.data.frame(file_paths$linx_breakpoints_hartwig) %>%
  rename('linx_breakpoints_catalog' = `file_paths$linx_breakpoints_hartwig`)
linx_breakpoints_pcawg <- as.data.frame(file_paths$linx_breakpoints_pcawg) %>%
  rename('linx_breakpoints_catalog' = `file_paths$linx_breakpoints_pcawg`)

# bind columns, add to the list
file_paths <- rbind(linx_driver_hartwig, linx_driver_pcawg)
file_paths[2] <- rbind(linx_fusion_hartwig, linx_fusion_pcawg)
file_paths[3] <- rbind(linx_breakpoints_hartwig, linx_breakpoints_pcawg)

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
if (dir.exists('/local/path')) {
  export_file_paths <- as.data.frame(lapply(file_paths_2, function(i){ gsub('/local/path','',i) })) 
} else {
  export_file_paths <- as.data.frame(lapply(file_paths_2, function(i){ gsub('/other/local/path','',i) }))
}
export_file_paths <- cbind(sample=metadata$sample_id, export_file_paths)

# create dir
dir.create(paste0(base_dir,'/processed/metadata/manifest'), recursive = TRUE)

# write to disk
write.table(
  export_file_paths, gzfile(paste0(base_dir,'/processed/metadata/manifest/', current_date, '_linx_output_paths.txt.gz')),
  sep='\t',quote=F,row.names=F
)
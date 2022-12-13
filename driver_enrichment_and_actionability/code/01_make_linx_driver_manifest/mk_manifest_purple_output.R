############### Create the purple output manifest file ############
# author: remy (sascha)
# date: 01/06/2021
# last updated: 22/10/2021

### Description
# This script creates a manifest file of the purple output file paths including VCF, purity and others.

### Input
# The metadata.tsv file & the Hartwig sample metadata file that you get when you get access to the Hartwig dataset (should include a setName column, e.g. 210831_HMFregACTN_FR30543039_FR30543053_ACTN01020002)

### Output
# A manifest file that contains the file paths to all the PURPLE VCF files of your dataset.

### Usage
# sbatch mk_manifest_purple_output.sh -d <data_request_version> -m <metadata_folder_name> -p <purple_version> -r <postprocessed_vcf: yes/no>
# e.g. sbatch 8_make_hartwig_vcf_manifest/mk_manifest_purple_output.sh -d DR-104-update4 -m metadata_2021 -p purplesoft3.3 -r yes
# <data_request_version> is the version of your data request from Hartwig Medical Foundation
# <metadata_folder_name> the folder where the Hartwig metadata file is stored
# <purple_version> the folder name of your PURPLE output files (must be the same name per sample)
# <postprocessed_vcf: yes/no> if you ran the postprocess_pcawg_vcf;sh script then this should be 'yes'

# global options
options(stringsAsFactors=F)

# libs
source(paste0(here::here(), '/code/r_objects/libs.R'))

# pass commandline directory arguments to 'args' vector
args = commandArgs(trailingOnly=TRUE)

# parse commandline input
data_request_version <- args[1]
metadata_folder <- args[2]
purple_version <- args[3]
postprocessed_vcf <- args[4]

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

hartwig_data_dir <- list(
  path=paste0(here::here(), '/data/raw/', data_request_version, '/', metadata_folder)
)

for(i in hartwig_data_dir){
  if(dir.exists(i)){
    hartwig_data_dir <- i
    break
  }
}

#========= Main =========#
# read the metadata in, make all the column names lower case
metadata <- read_tsv(paste0(base_dir, '/data/processed/metadata/metadata.tsv'),
                     col_types = 'icccccccccccccccllllcccdiiildclccllllllliiiiiiiiiii') %>%
  rename_with(tolower)

# read the hartwig metadata in (variable setName is important)
hartwig_metadata <- read_tsv(paste0(base_dir, '/data/processed/metadata/hartwig_metadata.tsv')) %>%
  select(hmfSampleId, setName)

# join setName from the hartwig_metadata to the metadata table
metadata <- metadata %>%
  left_join(hartwig_metadata, by = c("sample_id" = "hmfSampleId")) %>%
  relocate(setName, .after = sample_id)

# split by cohort
metadata_hartwig <- metadata %>%
  filter(cohort == "Hartwig")
metadata_pcawg <- metadata %>%
  filter(cohort == "PCAWG")

hartwig_data_dir <- hartwig_data_dir
if(dir.exists('/local/path')){
  hartwig_data_dir <- paste0('/local/path',
                         hartwig_data_dir)
}

pcawg_hartwig_data_dir <- '/path/to/pcawg/samples'
if(dir.exists('/local/path')){
  pcawg_hartwig_data_dir <- paste0('/local/path',
                               pcawg_hartwig_data_dir)
}

# initiate list to store the file paths in
file_paths <- list()

# add hartwig hartwig purple vcf file path to list
file_paths$purple_hartwig <- with(metadata_hartwig,{ paste0(hartwig_data_dir, '/somatics/', sample_id,
                                                    '/purple/', sample_id, '.purple.somatic.vcf.gz') })
# add pcawg hartwig purple vcf file path to list
if (postprocessed_vcf == 'yes') {
  file_paths$purple_pcawg_hartwig <- with(metadata_pcawg,{ paste0(pcawg_hartwig_data_dir, '/', sample_id, '/', purple_version, '/', sample_id, 'T.purple.somatic.postprocessed.vcf.gz') })
} else {
  file_paths$purple_pcawg_hartwig <- with(metadata_pcawg,{ paste0(pcawg_hartwig_data_dir, '/', sample_id, '/', purple_version, '/', sample_id, 'T.purple.somatic.vcf.gz') })
}

# add purity hartwig hartwig file path to the list
file_paths$purity_hartwig <- with(metadata_hartwig,{ paste0(hartwig_data_dir, '/somatics/', sample_id,
                                                    '/purple/', sample_id, '.purple.purity.tsv') })
# add purity pcawg hartwig file path to the list
file_paths$purity_pcawg <- with(metadata_pcawg,{ paste0(pcawg_hartwig_data_dir, '/', sample_id, '/', purple_version, '/', sample_id, 'T.purple.purity.tsv') })

# transform list into a dataframe and change the column names
purple_hartwig <- as.data.frame(file_paths$purple_hartwig) %>%
  rename("purple_vcf" = `file_paths$purple_hartwig`)

purple_pcawg <- as.data.frame(file_paths$purple_pcawg_hartwig) %>%
  rename("purple_vcf" = `file_paths$purple_pcawg_hartwig`)

purity_hartwig <- as.data.frame(file_paths$purity_hartwig) %>%
  rename("purity" = `file_paths$purity_hartwig`)

purity_pcawg <- as.data.frame(file_paths$purity_pcawg) %>%
  rename("purity" = `file_paths$purity_pcawg`)

# bind columns, add to the list
file_paths <- rbind(purple_hartwig, purple_pcawg)
file_paths[3] <- rbind(purity_hartwig, purity_pcawg)

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
   #keys=c('germ1','germ2')
   idx <- apply(exist_files[,keys],1,which)

   file_paths_ss <- file_paths[,keys]
   counter <- 0
   out <- unname(apply(file_paths_ss,1,function(i){
      counter <<- counter+1
      i[ idx[counter] ]
   }))

   return(out)
}

# initiate second list to store paths to files that exist
file_paths_2 <- list()

# Hartwig
file_paths_2$purple_vcf <- file_paths$purple_vcf
file_paths_2$purity <- file_paths$purity

#========= Export =========#
export_file_paths <- as.data.frame(lapply(file_paths_2, function(i){ gsub('/local/path','',i) }))
export_file_paths <- cbind(sample=metadata$sample_id, export_file_paths)

# write to disk
write.table(
  export_file_paths, 
  file = if (postprocessed_vcf == 'yes') {
    gzfile(paste0(base_dir,'/data/processed/metadata/manifest/purple_output_paths_postprocessed.txt.gz'))
  } else {
    gzfile(paste0(base_dir,'/data/processed/metadata/manifest/purple_output_paths.txt.gz'))
  },
   sep='\t',quote=F,row.names=F
)
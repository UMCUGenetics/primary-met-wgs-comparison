############### Create the VCF manifest file ############
# author: remy (sascha)
# date: 07/05/2021
# last updated: 16/02/2022

### Description
# This script creates a manifest file of the Hartwig (= HMF) and PCAWG somatic VCF file paths to be used for
# subsequent analysis.

### Input
# metadata.tsv table

### Output
# A manifest file: vcf_paths.txt.gz

### Usage
# sbatch mk_manifest.sh

# global options
options(stringsAsFactors=F)

# libs
here::set_here(path='../../')
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

#========= Main =========#
# load metadata, rename column names to lower case, filter the table for PCAWG cohort only
metadata <- read_tsv(paste0(base_dir, '/data/processed/metadata/metadata.tsv')) %>%
  rename_with(tolower) %>%
  filter(cohort == 'PCAWG')

# hartwig processed pcawg vcfs base path
hartwig_data_dir <- '/path/to/hartwig/VCFs/'
if(dir.exists('/local/path')){
   hartwig_data_dir <- paste0('/local/path', hartwig_data_dir)
}

# pcawg processed pcawg vcfs base path
pcawg_data_dir <- '/path/to/PCAWG/VCFs/'
if(dir.exists('/local/path')){
  pcawg_data_dir <- paste0('/local/path', pcawg_data_dir)
}

# initiate list to stare individual VCF file paths
file_paths <- list()

## Somatic SNV/MNV/Indel
file_paths$som_hartwig <- with(metadata,{ paste0(hartwig_data_dir, '/', sample_id, '/', sample_id, 'T.purple.somatic.postprocessed.vcf.gz') }) # postprocessed VCF files have removed duplicates and only include the highest base substitution at each position (e.g. if there is a SBS and DBS reported at the same position then only the DBS is kept)
file_paths$som_pcawg_snv <- with(metadata,{ paste0(pcawg_data_dir, '/', sample_id, '/', sample_id, '.consensus.20160830.somatic.snv_mnv_PASS.vcf.gz') })
file_paths$som_pcawg_indel <- with(metadata,{ paste0(pcawg_data_dir, '/', sample_id, '/', sample_id, '.consensus.20161006.somatic.indel_PASS.vcf.gz') })

## SV 
file_paths$sv_hartwig <- with(metadata,{ paste0(hartwig_data_dir, '/', sample_id, '/', sample_id, 'T.purple.sv.vcf.gz') })
file_paths$sv_pcawg <- with(metadata,{ paste0(pcawg_data_dir, '/', sample_id, '/', sample_id, '.pcawg_consensus_1.6.161116.somatic.sv.vcf.gz') })

# transform list to data frame, give each row the corresponding sample_id
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

# initate second list to store those file paths that actually exist on the HPC
file_paths_2 <- list()

# fill the list
file_paths_2$som_hartwig <- file_paths$som_hartwig
file_paths_2$som_pcawg_snv <- file_paths$som_pcawg_snv
file_paths_2$som_pcawg_indel <- file_paths$som_pcawg_indel
file_paths_2$sv_hartwig <- file_paths$sv_hartwig
file_paths_2$sv_pcawg <- file_paths$sv_pcawg
# transform to data frame
file_paths_2 <- as.data.frame(file_paths_2)

#========= Export =========#
vcf_paths <- as.data.frame(lapply(file_paths_2, function(i){ gsub('/local/path','',i) }))
vcf_paths <- cbind(sample=metadata$sample_id, vcf_paths)

message(paste0('Generating file... ', base_dir, '/data/processed/metadata/manifest/vcf_paths.txt.gz'))

# create output dir
dir.create(path = paste0(base_dir, '/data/processed/metadata/manifest'), recursive = TRUE)

# write manifest to disk
write.table(
   vcf_paths, gzfile(paste0(base_dir,'/data/processed/metadata/manifest/', current_date, '_vcf_paths.txt.gz')),
   sep='\t',quote=F,row.names=F
)

## Init --------------------------------

expanding_path <- '/home/ali313/Documents/studies/master/umc-project'

path_prefix <- c(
  hpc  ='/hpc/',
  local = paste0(expanding_path, '/hpc/')
)
path_prefix <- path_prefix[min(which(dir.exists(path_prefix)))]
path_prefix <- sub('/hpc/','',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')


## Loading basic packages --------------------------------

library(readxl)
library(stringr)
library(dplyr)
library(tidyverse)
library(stringr)

## Define functions --------------------------------

`%notin%` <- Negate(`%in%`)


## Assign values --------------------------------


patsam_ids <- list(PD11458 = c("a", "b", "c"), PD11459 = c("a", "b", "c"), PD11460 = c("a", "b", "c", "d") ,
                   PD11461 = c("a", "b", "c"), PD13596 = c("a", "b", "c"),
                   PD14780 = c("a", "b", "d", "e"), PD4243 = c("a", "b", "c"), PD4248 = c("a", "b", "c"),
                   PD4252 = c("a", "b", "c"), PD4820 = c("a", "b", "c"),
                   PD5956 = c("a", "b", "c"), PD6728 = c("a", "b", "c"), PD8948 = c("b", "d", "e"),
                   PD9193 = c("a", "b", "c"), PD9194 = c("a", "b", "c"),
                   PD9195 = c("a", "b", "c", "d"), PD9771 = c("a", "b", "c", "d", "e"))


patient_ids <- names(patsam_ids)



sample_ids <- c()


for (i in 1:length(patsam_ids)){
  sample_ids <- c(sample_ids, paste0(patient_ids[i], patsam_ids[[i]]))
}


# read in the metadata of Yates et al.

meta_lyates_raw <- as.data.frame(read_xlsx(path = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/lyates_cancer_cell_2014/metadata.xlsx"), sheet = 5,skip = 1))

# only keep selected samples
meta_lyates <- meta_lyates_raw[meta_lyates_raw$sample %in% sample_ids,]

# homogenize the tumor cohorts
meta_lyates$SIMPLIFIED_SAMPLE_CODE <- NA

meta_lyates$SIMPLIFIED_SAMPLE_CODE[meta_lyates$SAMPLE_CODE == "DISTANT_METASTASIS"] <- "METASTASIS"
meta_lyates$SIMPLIFIED_SAMPLE_CODE[meta_lyates$SAMPLE_CODE == "LOCAL_RELAPSE"] <- "METASTASIS"
meta_lyates$SIMPLIFIED_SAMPLE_CODE[meta_lyates$prefix == "PD11459" & meta_lyates$sample == "PD11459a"] <- "METASTASIS"
meta_lyates$SIMPLIFIED_SAMPLE_CODE[meta_lyates$SAMPLE_CODE == "PRIMARY"] <- "PRIMARY"

meta_lyates <- meta_lyates[!is.na(meta_lyates$SIMPLIFIED_SAMPLE_CODE),]
nrow(meta_lyates)
# only keep samples that have both primary and metastatic samples available
for (pat in unique(meta_lyates$prefix)){
  print(pat)
  if ("METASTASIS" %in% meta_lyates$SIMPLIFIED_SAMPLE_CODE[meta_lyates$prefix == pat] & "PRIMARY" %in% meta_lyates$SIMPLIFIED_SAMPLE_CODE[meta_lyates$prefix == pat]){
    
  } else {
    meta_lyates <- meta_lyates[meta_lyates$prefix != pat,]
  }
}

# remove extra (redundant) samples
for (pat in unique(meta_lyates$prefix)){
  print(pat)
  print(meta_lyates$SIMPLIFIED_SAMPLE_CODE[meta_lyates$prefix == pat])
}

meta_lyates <- meta_lyates[meta_lyates$sample %notin% c("PD9771a", "PD9771c"),]
meta_lyates <- meta_lyates[meta_lyates$sample %notin% "PD9195d",]



# update the patient/sample info based on the above blacklisting criteria

patsam_ids <- patsam_ids[names(patsam_ids) %in% unique(meta_lyates$prefix)]

patient_ids <- names(patsam_ids)

sample_ids <- c()


for (i in 1:length(patsam_ids)){
  sample_ids <- c(sample_ids, paste0(patient_ids[i], patsam_ids[[i]]))
}

sample_ids <- sample_ids[sample_ids %in% meta_lyates$sample]




## Init --------------------------------

expanding_path <- '/home/ali313/Documents/studies/master/umc-project'

path_prefix <- c(
  hpc  ='/hpc/',
  local = paste0(expanding_path, '/hpc/')
)
path_prefix <- path_prefix[min(which(dir.exists(path_prefix)))]
path_prefix <- sub('/hpc/','',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')

## Dependencies + defined values --------------------------------

source(paste0(base_dir, 'drivers/analysis/dna-rep-ann/final-update/code_for_github/source_script.R'))


## Sample metadata --------------------------------
sample_metadata <- read.delim(
  paste0(base_dir,'/passengers/processed/metadata/main/metadata_20221111_1251.txt.gz'),
  na.strings=c('NA','#N/A')
)


# prepare the metadata

sample_metadata <- sample_metadata[!sample_metadata$is_blacklisted,]
row.names(sample_metadata) <- 1:nrow(sample_metadata)

sample_metadata[,"cohort"] <- factor(sample_metadata[,"cohort"], levels = cohort_order)
sample_metadata[,"tissue_group"] <- factor(sample_metadata[,"tissue_group"], levels = tissue_group_order)
sample_metadata[,"cancer_type"] <- factor(sample_metadata[,"cancer_type"], levels = cancer_type_order)
sample_metadata[,"cancer_type_code"] <- factor(sample_metadata[,"cancer_type_code"], levels = cancer_type_code_order)
sample_metadata[,"gender"] <- factor(sample_metadata[,"gender"], levels = c("MALE", "FEMALE"))

sample_metadata$total_tmb <- rowSums(sample_metadata[,c("sbs_load", "dbs_load", "indel_load")])




## Summarizing the timing data of all mutations  --------------------------------


all_timing <- data.frame()
for (i in 1:nrow(sample_metadata)){


  sample_id <- sample_metadata$sample_id[i]


  cohort <- tolower(as.vector(sample_metadata$cohort[i]))
  cancer_type <- as.character(sample_metadata$cancer_type[i])

  print(i)
  print(sample_id)


  # assign the path to timing information
  if (sample_metadata$cohort[i] == "Hartwig") {
    file_path <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/timing/", sample_id, "/", sample_id, ".mutationaltiming.tsv.gz")
  } else if (sample_metadata$cohort[i] == "PCAWG") {
    file_path <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/timing/", sample_id, "/", sample_id, ".mutationaltiming.tsv.gz")
  }
  
  # read in the timing data
  combined_df <- tryCatch(read.csv(file = file_path, header = T, sep = "\t", stringsAsFactors = F), error=function(e) NULL)
  
  # count mutation number for each timing group
  tmb_clonal_late <- nrow(combined_df[combined_df$timing_class == "clonal [late]",])
  tmb_clonal_early <- nrow(combined_df[combined_df$timing_class == "clonal [early]",])
  tmb_clonal_na <- nrow(combined_df[combined_df$timing_class == "clonal [NA]",])
  tmb_mutationTimeR_subclonal <- nrow(combined_df[combined_df$timing_class == "subclonal",])

  if (!is.null(combined_df)){
    all_timing <- rbind(all_timing, c(sample_id, tmb_clonal_late, tmb_clonal_early, tmb_clonal_na, tmb_mutationTimeR_subclonal))
  } else {
    all_timing <- rbind(all_timing, c(sample_id, NA, NA, NA, NA))
  }


}

# assign the correct column names
colnames(all_timing) <- c("sample_id", "tmb_clonal_late", "tmb_clonal_early", "tmb_clonal_na", "tmb_mutationTimeR_subclonal")

# save the data
write.table(all_timing, file = gzfile(paste0(wd, "r-objects/rebuttal/all_timing.tsv.gz")), sep = "\t", row.names = F, quote = F)

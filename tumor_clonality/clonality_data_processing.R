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




## Summarize the clonality data  --------------------------------
no_rows <- nrow(sample_metadata)

purple_timing <- data.frame(sample_id = character(no_rows), clonal = numeric(no_rows), probably_clonal = numeric(no_rows), probably_subclonal = numeric(no_rows), subclonal = numeric(no_rows), clonal_tmb_80_cutoff = numeric(no_rows), subclonal_tmb_80_cutoff = numeric(no_rows))

for (i in 1:no_rows){
  
  sample_id <- sample_metadata$sample_id[i]
  
  print(i)
  print(sample_id)
  
  # assign the vcf paths
  if (sample_metadata$cohort[i] == "Hartwig") {
      path_to_vcf <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/HMF_data/DR-104-update5/somatics/", sample_id, "/purple/", sample_id, ".purple.somatic.vcf.gz")
  } else if (sample_metadata$cohort[i] == "PCAWG") {
      path_to_vcf <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/", sample_id, "-from-jar/purplesoft3.3/", sample_id, "T.purple.somatic.vcf.gz")
  }
  
  # read in the vcf file
  vcf <- tryCatch(variantsFromVcf(vcf.file = path_to_vcf,
                                  merge.consecutive = T,
                                  vcf.filter = "PASS",
                                  vcf.fields = c("CHROM", "POS", "REF", "ALT", "FILTER", "INFO")), error=function(e) NULL)
  
  if(!is.null(vcf)){
    if (nrow(vcf) > 0){
      
      # get the TNC and SUBCL information from the INFO field
      selelcted_info_fields <- getInfoValues(vcf$info, keys = c("TNC", "SUBCL"))
      
      selelcted_info_fields$SUBCL <- as.numeric(selelcted_info_fields$SUBCL)
      
      selelcted_info_fields$SUBCL[is.na(selelcted_info_fields$SUBCL)] <- 0
      
      # method 1: determining clonality probability using binned data
      percentiles <- quantile(selelcted_info_fields$SUBCL[selelcted_info_fields$SUBCL != 0], probs = c(0, 0.33, 0.66,1))
      
      tbl_subcl <- table(selelcted_info_fields$SUBCL)


      clonal <- sum(as.numeric(tbl_subcl["0"]))
      probably_clonal <- sum(as.numeric(tbl_subcl[-1][as.numeric(names(tbl_subcl))[-1] <= as.numeric(percentiles[2])]))
      probably_subclonal <- sum(as.numeric(tbl_subcl[-1][as.numeric(percentiles[2]) < as.numeric(names(tbl_subcl))[-1] &  as.numeric(names(tbl_subcl))[-1] <= as.numeric(percentiles[3])]))
      subclonal <- sum(as.numeric(tbl_subcl[-1][as.numeric(percentiles[3]) < as.numeric(names(tbl_subcl))[-1] &  as.numeric(names(tbl_subcl))[-1] <= as.numeric(percentiles[4])]))
      
      # method 2: determining clonality status using hard cutoff of 0.80 of subclonality score
      clonal_80 <- sum(selelcted_info_fields$SUBCL < 0.8)
      subclonal_80 <- sum(selelcted_info_fields$SUBCL >= 0.8)

      purple_timing[i,"sample_id"] <- sample_id
      purple_timing[i,2:7] <- c(clonal, probably_clonal, probably_subclonal, subclonal, clonal_80, subclonal_80)

    } else {
      purple_timing[i,"sample_id"] <- sample_id
      purple_timing[i,2:7] <- rep(NA, 6)
    }
  } else {
    purple_timing[i,"sample_id"] <- sample_id
    purple_timing[i,2:7] <- rep(NA, 6)
  }
}


# save the data
write.table(purple_timing, file = gzfile(paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/all-purple-timing.txt.gz")), quote = F, row.names = F, sep = "\t", col.names = F)



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




## Select clocklike mutations for each sample and add the timing info  --------------------------------

for (i in 1:nrow(sample_metadata)){


  sample_id <- sample_metadata$sample_id[i]


  cohort <- tolower(as.vector(sample_metadata$cohort[i]))
  cancer_type <- as.character(sample_metadata$cancer_type[i])


  # assign the vcf paths
  if (sample_metadata$cohort[i] == "Hartwig") {
    path_to_vcf <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/HMF_data/DR-104-update5/somatics/", sample_id, "/purple/", sample_id, ".purple.somatic.vcf.gz")
  } else if (sample_metadata$cohort[i] == "PCAWG") {
    path_to_vcf <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/", sample_id, "-from-jar/purplesoft3.3/", sample_id, "T.purple.somatic.vcf.gz")
  }
  


  # set the sample_id names correct
  if (cohort == "pcawg") {
    sample1 <- sample_id
    sample2 <- paste0(sample_id, "T")
  } else if (cohort == "hartwig"){
    if (substr(sample_id, nchar(sample_id),nchar(sample_id)) == "I" | substr(sample_id, nchar(sample_id),nchar(sample_id)) == "V"){
      sample_id2 <- str_split(sample_id, pattern = "TI")[[1]][1]
    } else {
      sample_id2 <- substr(sample_id, 1,(nchar(sample_id)-1))
    }
    sample1 <- sample_id2
    sample2 <- sample_id
  }
  
  # read in the vcf
  tryCatch(vcf <- variantsFromVcf(vcf.file = path_to_vcf,
                                  merge.consecutive = T,
                                  vcf.fields = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER","INFO", "FORMAT", paste0(sample1, "R"), sample2)), error=function(e) NULL)

  if(!is.null(vcf)){
    if (nrow(vcf) > 0){
      
      # filter for PASS mutation records
      vcf <- vcf[vcf$filter == "PASS",]
      
      # assign tri-nucleotide context for snvs
      vaf <- extractSigsSnv(df=vcf[,c(1:2,4:11,3)], output='df')
      vaf$chrom <- sapply(str_split(vaf$chrom, pattern = "r"), "[[", 2)

      vaf <- vaf[,c(1:2, 11, 3:10, 12:15)]
      
      # we consider clock like mutations as single base CpG > TpG mutations in NpCpG context. 
      # we excluded CpG > TpG in TpCpG which is also a characteristic peak in COSMIC SBS2 signature profile for Breast cohort
      # for skin melanoma CpG > TpG in [C/T]pCpG which overlaps with SBS7a was excluded
      if (cancer_type == "Breast cancer"){
        vaf <- vaf[grepl('\\w\\[C>T\\]G', vaf$context),]
        vaf <- vaf[substr(vaf$context, 1,1) !="T",]
      } else if (cancer_type == "Skin melanoma"){
        vaf <- vaf[grepl('\\w\\[C>T\\]G', vaf$context),]
        vaf <- vaf[substr(vaf$context, 1,1) %notin% c("C", "T"),]
      } else {
        vaf <- vaf[grepl('\\w\\[C>T\\]G', vaf$context),]
        vaf <- vaf[substr(vaf$context, 1,1) !="T",]
      }

      vaf <- vaf[,1:11]
      
      # format the data frame
      colnames(vaf) <- toupper(colnames(vaf))
      rownames(vaf) <- 1:nrow(vaf)

      selelcted_info_fields <- getInfoValues(vaf$INFO, keys = c("TNC", "SUBCL"))
      selelcted_info_fields$SUBCL[is.na(selelcted_info_fields$SUBCL)] <- 0
      selelcted_info_fields$SUBCL <- as.numeric(selelcted_info_fields$SUBCL)
      selelcted_info_fields$clonality <- NA
      selelcted_info_fields$clonality[which(selelcted_info_fields$SUBCL >= 0.8)] <- "subclonal"
      selelcted_info_fields$clonality[which(selelcted_info_fields$SUBCL < 0.8)] <- "clonal"

      vaf <- cbind(vaf, selelcted_info_fields)

    }
  }

  
  if (sample_metadata$cohort[i] == "Hartwig") {
      file_path <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/timing/", sample_id, "/", sample_id, ".mutationaltiming.tsv.gz")
  } else if (sample_metadata$cohort[i] == "PCAWG") {
      file_path <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/timing/", sample_id, "/", sample_id, ".mutationaltiming.tsv.gz")
  }

  # read in the result of mutationtimeR (contains timing info for mutations)
  combined_df <- tryCatch(read.csv(file = file_path, header = T, sep = "\t", stringsAsFactors = F), error=function(e) NULL)

  # add the timing info to the mutation data frame
  if (!is.null(combined_df)){
    vaf$coorPos <- paste0(vaf$CHROM, "_", vaf$REF, "_", vaf$POS, "_", vaf$ALT)
    combined_df$coorPos <- paste0(combined_df$CHROM, "_", combined_df$REF, "_", combined_df$POS, "_", combined_df$ALT)
    vaf <- merge(vaf, combined_df[,c("coorPos", "timing_class")], by = "coorPos")
    vaf <- vaf[,-1]
  } else {
    vaf$timing_class <- NA
  }
  
  # save the data
  write.table(vaf, file = gzfile(paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/clock-like-vafs-age-investigation/", sample_id, ".purple.somatic.clocklike-subset.vcf.gz")), sep = "\t", row.names = F, quote = F)

}


## Summarize the counts of clocklike mutations and their timing categories  --------------------------------

sbs1_timing <- data.frame()

for (i in 1:nrow(sample_metadata)){

  sample_id <- sample_metadata$sample_id[i]
  print(i)

  # read in the clocklike mutation timing data 
  clock_like_muts <- tryCatch(read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/clock-like-vafs-age-investigation/", sample_id, ".purple.somatic.clocklike-subset.vcf.gz"), sep = "\t", header = T, stringsAsFactors = F), error=function(e) NULL)
  if(is.null(clock_like_muts)){
    sbs1_count <- NA
    sbs1_subcl_subclonal <- NA
    sbs1_clonal_late <- NA
    sbs1_clonal_early <- NA
    sbs1_na_clonal <- NA
    sbs1_mutationTimeR_subclonal <- NA
  } else {
    sbs1_count <- nrow(clock_like_muts)
    sbs1_subcl_subclonal <- nrow(clock_like_muts[clock_like_muts$clonality == "subclonal",])
    if (all(is.na(clock_like_muts$timing_class))) {
      sbs1_clonal_late <- NA
      sbs1_clonal_early <- NA
      sbs1_na_clonal <- NA
      sbs1_mutationTimeR_subclonal <- NA
    } else {
      sbs1_clonal_late <- nrow(clock_like_muts[clock_like_muts$timing_class == "clonal [late]",])
      sbs1_clonal_early <- nrow(clock_like_muts[clock_like_muts$timing_class == "clonal [early]",])
      sbs1_na_clonal <- nrow(clock_like_muts[clock_like_muts$timing_class == "clonal [NA]",])
      sbs1_mutationTimeR_subclonal <- nrow(clock_like_muts[clock_like_muts$timing_class == "subclonal",])
    }
  }

  sbs1_timing <- rbind(sbs1_timing, c(sample_id, sbs1_count, sbs1_subcl_subclonal, sbs1_clonal_late, sbs1_clonal_early,
                                      sbs1_na_clonal, sbs1_mutationTimeR_subclonal))

}

# assign the correct column names
colnames(sbs1_timing) <- c("sample_id", "sbs1_count", "sbs1_subcl_subclonal", "sbs1_clonal_late", "sbs1_clonal_early",
                         "sbs1_na_clonal", "sbs1_mutationTimeR_subclonal")

# save the sumary data
write.table(sbs1_timing, file = gzfile(paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/rebuttal/sbs1_timing.tsv.gz")), quote = F, row.names = F, sep = "\t")



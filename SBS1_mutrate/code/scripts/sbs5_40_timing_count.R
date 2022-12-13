
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




## Summarizing the timing data of SBS5/40 mutations (based on mutation signature analysis) --------------------------------

sbs5_timing <- data.frame()
sbs40_timing <- data.frame()

for (i in 1:nrow(sample_metadata)){  
  
  
  sample_id <- sample_metadata$sample_id[i]
  
  
  cohort <- tolower(as.vector(sample_metadata$cohort[i]))
  cancer_type <- as.character(sample_metadata$cancer_type[i])
  
  # for sbs 5 and 40 mutations we rely on the outcome of signature analysis
  path_to_ms <- paste0(base_dir, "/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/12_fixed_seeds/sig_contrib/muts_assigned/", sample_id, ".txt.gz")
  
  # read in the mutation table with assigned signatures
  tryCatch(ms_df <- read.csv(file = path_to_ms, header = T, stringsAsFactors = F, sep = "\t"), error=function(e) NULL)
  
  ms_df <- ms_df[!is.na(ms_df$sig_name),]
  
  if (nrow(ms_df) > 0){
    ms_df$chrom <- sapply(str_split(ms_df$chrom, pattern = "r"), "[[", 2)
    
    # assign the path to timing information
    if (sample_metadata$cohort[i] == "Hartwig") {
        file_path <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/timing/", sample_id, "/", sample_id, ".mutationaltiming.tsv.gz")
      } else if (sample_metadata$cohort[i] == "PCAWG") {
          file_path <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/timing/", sample_id, "/", sample_id, ".mutationaltiming.tsv.gz")
      }
      
      # read in the timing data
      combined_df <- tryCatch(read.csv(file = file_path, header = T, sep = "\t", stringsAsFactors = F), error=function(e) NULL)
    
    
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
    vcf <- tryCatch(variantsFromVcf(vcf.file = path_to_vcf,
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
        vaf <- vaf[,1:11]
        
        colnames(vaf) <- toupper(colnames(vaf))
        rownames(vaf) <- 1:nrow(vaf)
        
        # get the TNC and SUBCL information from the INFO field and use the 0.8 cut off score to assign their clonality
        selelcted_info_fields <- getInfoValues(vaf$INFO, keys = c("TNC", "SUBCL"))
        selelcted_info_fields$SUBCL[is.na(selelcted_info_fields$SUBCL)] <- 0
        selelcted_info_fields$SUBCL <- as.numeric(selelcted_info_fields$SUBCL)
        selelcted_info_fields$clonality <- NA
        selelcted_info_fields$clonality[which(selelcted_info_fields$SUBCL >= 0.8)] <- "subclonal"
        selelcted_info_fields$clonality[which(selelcted_info_fields$SUBCL < 0.8)] <- "clonal"
        
        # combine vcf data and clonality data
        vaf <- cbind(vaf, selelcted_info_fields)
        
        
        vaf$coorPos <- paste0(vaf$CHROM, "_", vaf$REF, "_", vaf$POS, "_", vaf$ALT)
        
      }
    }
    
    
    
    ms_df$coorPos <- paste0(ms_df$chrom, "_", ms_df$ref, "_", ms_df$pos, "_", ms_df$alt)
    
    if (!is.null(combined_df)){
      combined_df$coorPos <- paste0(combined_df$CHROM, "_", combined_df$REF, "_", combined_df$POS, "_", combined_df$ALT)
      ms_df <- merge(ms_df, combined_df[,c("coorPos", "timing_class")], by = "coorPos")
    } else {
      ms_df$timing_class <- NA
    }
    
    # add clonality data to assigned signature data
    ms_df <- merge(ms_df, vaf[,c("coorPos", "clonality")], by = "coorPos")
    
    ms_df <- ms_df[,-1]
    
    # select for mutations assigned to SBS5 and SBS40
    ms_df_sbs5 <- ms_df[ms_df$sig_name == "SBS5",]
    ms_df_sbs40 <- ms_df[ms_df$sig_name == "SBS40",]
    
    # count numer of mutations for their timing of SBS5
    sbs5_count <- nrow(ms_df_sbs5)
    sbs5_subcl_subclonal <- nrow(ms_df_sbs5[ms_df_sbs5$clonality =="subclonal",])
    if (all(is.na(ms_df_sbs5$timing_class))){
      sbs5_clonal_late <- NA
      sbs5_clonal_early <- NA
      sbs5_clonal_na <- NA
      sbs5_timing_subclonal <- NA
      
    } else {
      sbs5_clonal_late <- nrow(ms_df_sbs5[ms_df_sbs5$timing_class =="clonal [late]",])
      sbs5_clonal_early <- nrow(ms_df_sbs5[ms_df_sbs5$timing_class =="clonal [early]",])
      sbs5_clonal_na <- nrow(ms_df_sbs5[ms_df_sbs5$timing_class =="clonal [NA]",])
      sbs5_timing_subclonal <- nrow(ms_df_sbs5[ms_df_sbs5$timing_class =="subclonal",])
    }
    
    sbs5_timing <- rbind(sbs5_timing, c(sample_id,sbs5_count,sbs5_subcl_subclonal,sbs5_clonal_late,sbs5_clonal_early,sbs5_clonal_na,sbs5_timing_subclonal))
    
    # count numer of mutations for their timing of SBS40
    sbs40_count <- nrow(ms_df_sbs40)
    sbs40_subcl_subclonal <- nrow(ms_df_sbs40[ms_df_sbs40$clonality =="subclonal",])
    if (all(is.na(ms_df_sbs5$timing_class))){
      sbs40_clonal_late <- NA
      sbs40_clonal_early <- NA
      sbs40_clonal_na <- NA
      sbs40_timing_subclonal <- NA
    } else {
      sbs40_clonal_late <- nrow(ms_df_sbs40[ms_df_sbs40$timing_class =="clonal [late]",])
      sbs40_clonal_early <- nrow(ms_df_sbs40[ms_df_sbs40$timing_class =="clonal [early]",])
      sbs40_clonal_na <- nrow(ms_df_sbs40[ms_df_sbs40$timing_class =="clonal [NA]",])
      sbs40_timing_subclonal <- nrow(ms_df_sbs40[ms_df_sbs40$timing_class =="subclonal",])
    }
    
    sbs40_timing <- rbind(sbs40_timing, c(sample_id,sbs40_count,sbs40_subcl_subclonal,sbs40_clonal_late,sbs40_clonal_early,sbs40_clonal_na,sbs40_timing_subclonal))
    
  } else {
    sbs5_timing <- rbind(sbs5_timing, c(sample_id,0,0,0,0,0,0))
    sbs40_timing <- rbind(sbs40_timing, c(sample_id,0,0,0,0,0,0))
    
    
  }
  
}


# assign the correct column names
colnames(sbs5_timing) <- c('sample_id',"sbs5_count","sbs5_subcl_subclonal","sbs5_clonal_late","sbs5_clonal_early","sbs5_clonal_na","sbs5_mutationTimeR_subclonal")
colnames(sbs40_timing) <- c('sample_id',"sbs40_count","sbs40_subcl_subclonal","sbs40_clonal_late","sbs40_clonal_early","sbs40_clonal_na","sbs40_mutationTimeR_subclonal")


# save the data
write.table(sbs5_timing, file = gzfile(paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/rebuttal/sbs5_timing.tsv.gz")), sep = "\t", row.names = F, quote = F)
write.table(sbs40_timing, file = gzfile(paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/rebuttal/sbs40_timing.tsv.gz")), sep = "\t", row.names = F, quote = F)



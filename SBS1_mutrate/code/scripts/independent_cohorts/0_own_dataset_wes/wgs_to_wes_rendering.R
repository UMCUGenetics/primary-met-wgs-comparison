
### In this script we restrict the VCF files of WGS PCAWG and Hartwig dataset to WES defined regions for comparison with external WES datasets



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

source(paste0(base_dir, 'drivers/analysis/dna-rep-ann/final-update/code_for_github/SBS1/independent_cohorts/0_own_dataset_wes/wgs_to_wes_rendering_source.R'))


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

sample_metadata$tmb_total <- rowSums(sample_metadata[,c("sbs_load", "dbs_load", "indel_load")])
sample_metadata$subclonal_tmb <- rowSums(sample_metadata[,c("sbs_load.subclonal", "dbs_load.subclonal", "indel_load.subclonal")])
sample_metadata$clonal_tmb <- sample_metadata$tmb_total - sample_metadata$subclonal_tmb


## Add age info  --------------------------------

# read in and process raw Hartwig age information


hartwig_meta <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external-files/hmf-clinical-metadata-20211112.tsv"), stringsAsFactors = F, header = T, sep = "\t")

hartwig_meta$biopsyDateProcessed <- NA
hartwig_meta$biopsyDateProcessed <- unlist(lapply(str_split(hartwig_meta$biopsyDate, pattern = "-"), "[[", 1))

hartwig_meta$age_at_biopsy <- NA

hartwig_meta[,c("birthYear","biopsyDateProcessed")] <- apply(hartwig_meta[,c("birthYear","biopsyDateProcessed")], as.numeric, MARGIN = 2)

hartwig_meta$age <- hartwig_meta$biopsyDateProcessed - hartwig_meta$birthYear


colnames(hartwig_meta)[2] <- "sample_id"

hartwig_meta <- hartwig_meta[,c(2,ncol(hartwig_meta))]




# read in and process raw PCWG age information


pcawg_meta <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/PCAWG/Metadata/donor.info.followup.all_projects.csv"), header = T, stringsAsFactors = F, sep = "\t")


colnames(pcawg_meta)[c(1,9)] <- c("sample_id", "age")
pcawg_meta <- pcawg_meta[,c(1,9)]




# merge the age info of both cohorts to the sample metadata object
age_meta_all <- rbind(hartwig_meta, pcawg_meta)
sample_metadata <- merge(sample_metadata, age_meta_all, by = "sample_id")




# define cancer types that we are interested in their WES mutation 
cancer_type_codes_of_interest <- c("BRCA", "COREAD", "OV", "KIRC", "PRAD", "THCA")


# select and format vcf files of samples that are in the cacner type of interests

for (i in 1:length(cancer_type_codes_of_interest)){
  cancer_type_code <- cancer_type_codes_of_interest[i]
  
  sub_sample_metadata <- sample_metadata[sample_metadata$cancer_type_code == cancer_type_code,]

  for (j in 1:nrow(sub_sample_metadata))
  
    sample_id <- sub_sample_metadata$sample_id[j]
  
    cohort <- tolower(as.vector(sub_sample_metadata$cohort[j]))
    cancer_type <- as.character(sub_sample_metadata$cancer_type[j])
  
    
    # assign the vcf paths
    if (sample_metadata$cohort[i] == "Hartwig") {
      path_to_vcf <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/HMF_data/DR-104-update5/somatics/", sample_id, "/purple/", sample_id, ".purple.somatic.vcf.gz")
    } else if (sample_metadata$cohort[i] == "PCAWG") {
      path_to_vcf <- paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/", sample_id, "-from-jar/purplesoft3.3/", sample_id, "T.purple.somatic.vcf.gz")
    }
    
    
  
    # read in the vcf
    tryCatch(vcf <- variantsFromVcf(vcf.file = path_to_vcf,
                           merge.consecutive = T,
                           vcf.fields = c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER","INFO", "FORMAT")), error=function(e) NULL)
  
      if (!is.null(vcf)){
      
      # format the mutation calls
      colnames(vcf)[1:2] <- c("chr", "pos1")
      vcf$pos2 <- vcf$pos1 + nchar(vcf$ref)
      vcf <- vcf[,c(1,2,ncol(vcf), 3:9)]
      vcf$chr <- sapply(strsplit(vcf$chr, split = "r"), "[[", 2)
      
      # save the data
      write.table(vcf, file = gzfile(paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/prostate_exome/raw_vcfs/", cancer_type_code, "/", j, ".tsv.gz")), sep = "\t", quote = F, row.names = F)
  
  }
}


  # -------------------------------------------------------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  ## Bash code to retrive WES data by intersecting WES ranges with mutation calls
  
  # for i in {1..nrow(sub_sample_metadata)}; do echo ${i} & ~/bedtools intersect -a ./raw_vcfs/${i}.tsv.gz -b S04380110_Regions_processed.bed.gz > ./exome_vcfs/${i}.tsv; done
  
  # -------------------------------------------------------------------------------------------------------------------------------------------------------------
  # -------------------------------------------------------------------------------------------------------------------------------------------------------------
  
  

# annotate tri-nucleotide context of WES mutations



# load required packages and objects
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
ref_genome <- getBSgenome(ref_genome)




# add the tri-nucleotide context of mutations
for (i in 1:length(cancer_type_codes_of_interest)){
  cancer_type_code <- cancer_type_codes_of_interest[i]
  
  sub_sample_metadata <- sample_metadata[sample_metadata$cancer_type_code == cancer_type_codes_of_interest,]
  
  for (j in 1:nrow(sub_sample_metadata)){
    print(paste0("===============================",k))
  
    path_to_wes <- paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/prostate_exome/exome_vcfs/", cancer_type_code, "/", k, ".tsv")
  
    vcf <- read.csv(file = path_to_wes, sep = "\t", stringsAsFactors = F, header = F)
    colnames(vcf) <- c("CHROM", "POS1", "POS2","ID", "REF", "ALT", "QUAL", "FILTER","INFO", "FORMAT")
    vcf <- vcf[vcf$FILTER == "PASS",]
    vcf <- vcf[nchar(vcf$REF) == 1 & nchar(vcf$ALT) == 1,]
  
    vcf$CHROM <- paste0("chr", vcf$CHROM)
  
  
    vcf$tri_context_ref <- as.character(Biostrings::getSeq(ref_genome, vcf$CHROM, vcf$POS1 - 1, vcf$POS1 + 1))
  
      for (j in 1:nrow(vcf)){
      print(j)
      if (vcf$REF[j] %in% c("G", "A")){
        vcf$tri_context_96[j] <- as.character(reverseComplement(DNAString(vcf$tri_context_ref[j])))
        vcf$change[j] <- paste0(as.character(complement(DNAString(vcf$REF[j]))), ">", as.character(complement(DNAString(vcf$ALT[j]))))
      } else {
        vcf$tri_context_96[j] <- vcf$tri_context_ref[j]
        vcf$change[j] <- paste0(vcf$REF[j], ">", vcf$ALT[j])
      }
  
      vcf$tri_context_96[j] <- paste0(substr(vcf$tri_context_96[j],1,1), "[", vcf$change[j], "]", substr(vcf$tri_context_96[j],3,3))
  
      }
    
    # save the data
    write.table(vcf, file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/prostate_exome/exome_vcfs/", cancer_type_code,"/annotated/", metadata_included_sub$sample_id[k], ".tsv"), row.names = F, quote = F, sep = "\t")
  
  }
}



# make a summary of sbs1 counts 



for (i in 1:length(cancer_type_codes_of_interest)){
  cancer_type_code <- cancer_type_codes_of_interest[i]
  
  sub_sample_metadata <- sample_metadata[sample_metadata$cancer_type_code == cancer_type_codes_of_interest,]
  
  for (i in 1:nrow(metadata_included_sub)){
  
    muts_df <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/prostate_exome/exome_vcfs/", cancer_type_code, "/annotated/", metadata_included_sub$sample_id[i], ".tsv"), stringsAsFactors = F, header = T, sep = "\t")
  
    tot_mut_no <- nrow(muts_df)
    
    # SBS1 context
    muts_df_sbs1 <- muts_df[grepl('\\w\\[C>T\\]G', muts_df$tri_context_96),]
    muts_df_sbs1 <- muts_df_sbs1[substr(muts_df_sbs1$tri_context_96, 1,1) !="T",]
  
    sbs1_mut_no <- nrow(muts_df_sbs1)
  
    metadata_included_sub$tot_count[i] <- tot_mut_no
    metadata_included_sub$sbs1_count[i] <- sbs1_mut_no
  
  
  }
  
  metadata_included_sub$sbs1_per_year <- metadata_included_sub$sbs1_count/metadata_included_sub$age
  
  
  write.table(metadata_included_sub, file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/prostate_exome/results/r-objects/prim_exom/metadata_included_", cancer_type_code, ".tsv"), sep = "\t", row.names = F, quote = F)
  
}



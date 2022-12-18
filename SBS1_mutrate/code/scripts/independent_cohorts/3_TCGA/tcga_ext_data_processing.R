
### In this script we process TCGA WES external data set 



## Init --------------------------------

expanding_path <- '/home/ali313/Documents/studies/master/umc-project'

path_prefix <- c(
  hpc  ='/hpc/',
  local = paste0(expanding_path, '/hpc/')
)
path_prefix <- path_prefix[min(which(dir.exists(path_prefix)))]
path_prefix <- sub('/hpc/','',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')



# define cancer types that we are interested in their WES mutation 
cancer_type_codes_of_interest <- c("BRCA", "COREAD", "OV", "KIRC", "PRAD", "THCA")


# load required packages
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
ref_genome <- getBSgenome(ref_genome)


for (i in 1:length(cancer_type_codes_of_interest)){
  
  cancer_type_code <- cancer_type_codes_of_interest[i]
  
  print(paste0("#####################################", cancer_type_code))
  
  # read in the meta
  meta_tcga <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/prostate_exome/tcga/metadata/", cancer_type_code, "/clinical.tsv"), sep = "\t", stringsAsFactors = F, header = T)
  
  # read in the mutation file
  mutations_tcga <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/prostate_exome/tcga/mutations/", cancer_type_code, ".maf.xz"), sep = "\t", stringsAsFactors = F, header = T)
  
  # subset for sbs
  mutations_tcga <- mutations_tcga[nchar(mutations_tcga$Reference_Allele) == 1 & nchar(mutations_tcga$Tumor_Seq_Allele2) == 1,]
  
  # format data
  mutations_tcga$Chromosome <- paste0("chr", mutations_tcga$Chromosome)
  mutations_tcga <- mutations_tcga[!is.na(mutations_tcga$Tumor_Sample_Barcode),]
  
  fi <- as.vector(sapply(str_split(mutations_tcga$Tumor_Sample_Barcode, pattern = "-"), "[[", 1))
  se <- sapply(str_split(mutations_tcga$Tumor_Sample_Barcode, pattern = "-"), "[[", 2)
  th <- sapply(str_split(mutations_tcga$Tumor_Sample_Barcode, pattern = "-"), "[[", 3)
  
  mutations_tcga$case_submitter_id <- paste(fi, se, th, sep = "-")
  
  # add tri-nucleotide context annotation
  mutations_tcga$tri_context_ref <- as.character(Biostrings::getSeq(ref_genome, mutations_tcga$Chromosome, mutations_tcga$Start_Position - 1, mutations_tcga$Start_Position + 1))
  
  
  for (j in 1:nrow(mutations_tcga)){
    print(j)
    if (mutations_tcga$Reference_Allele[j] %in% c("G", "A")){
      mutations_tcga$tri_context_96[j] <- as.character(reverseComplement(DNAString(mutations_tcga$tri_context_ref[j])))
      mutations_tcga$change[j] <- paste0(as.character(complement(DNAString(mutations_tcga$Reference_Allele[j]))), ">", as.character(complement(DNAString(mutations_tcga$Tumor_Seq_Allele2[j]))))
    } else {
      mutations_tcga$tri_context_96[j] <- mutations_tcga$tri_context_ref[j]
      mutations_tcga$change[j] <- paste0(mutations_tcga$Reference_Allele[j], ">", mutations_tcga$Tumor_Seq_Allele2[j])
    }
    
    mutations_tcga$tri_context_96[j] <- paste0(substr(mutations_tcga$tri_context_96[j],1,1), "[", mutations_tcga$change[j], "]", substr(mutations_tcga$tri_context_96[j],3,3))
    
  }
  
  
  
  # count and add number of sbs1 mutations to metadata
  
  
  for (i in 1:length(unique(meta_tcga$case_submitter_id))){
    print(i)
    sample <- unique(meta_tcga$case_submitter_id)[i]
    tot_mut_no <- nrow(mutations_tcga[mutations_tcga$case_submitter_id == sample,])
    
    # sbs1 context definition
    muts_df_sbs1 <- mutations_tcga[grepl('\\w\\[C>T\\]G', mutations_tcga$tri_context_96) & mutations_tcga$case_submitter_id == sample,]
    muts_df_sbs1 <- muts_df_sbs1[substr(muts_df_sbs1$tri_context_96, 1,1) !="T",]
    
    sbs1_mut_no <- nrow(muts_df_sbs1)
    
    # add counts to metadata
    meta_tcga$tot_count[meta_tcga$case_submitter_id == sample] <- tot_mut_no
    meta_tcga$sbs1_count[meta_tcga$case_submitter_id == sample] <- sbs1_mut_no
    
    
  }
  
  # sbs1 counts per year
  meta_tcga$sbs1_per_year <- meta_tcga$sbs1_count/meta_tcga$age_at_index
  
  meta_tcga[meta_tcga$tot_count == 0,159:161] <- NA
  
  
  # save the data
  write.table(meta_tcga, file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/prostate_exome/results/r-objects/tcga_exom/meta_", cancer_type_code, "_tcga_prim_with_counts.tsv"), sep = "\t", row.names = F, quote = F)
  
  
}


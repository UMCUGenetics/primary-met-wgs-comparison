
### In this script we process Diana Miao et al. WES external data set 



## Init --------------------------------

expanding_path <- '/home/ali313/Documents/studies/master/umc-project'

path_prefix <- c(
  hpc  ='/hpc/',
  local = paste0(expanding_path, '/hpc/')
)
path_prefix <- path_prefix[min(which(dir.exists(path_prefix)))]
path_prefix <- sub('/hpc/','',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')



# read in sample metadata of Diana Miao et al. external data set (Kidney wes metastatic)


meta_miao <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/kidney_exome/kidney_met_whole_exome_science2019/ccrcc_dfci_2019/data_clinical_sample.txt"), sep = "\t", stringsAsFactors = F, header = T, skip = 4) 



# load required packages
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
ref_genome <- getBSgenome(ref_genome)

# read in the mutation file
mutations_miao <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/kidney_exome/kidney_met_whole_exome_science2019/ccrcc_dfci_2019/data_mutations.txt"), sep = "\t", stringsAsFactors = F, header = T)


# subset for sbs
mutations_miao <- mutations_miao[nchar(mutations_miao$Reference_Allele) == 1 & nchar(mutations_miao$Tumor_Seq_Allele2) == 1,]
mutations_miao$Chromosome <- paste0("chr", mutations_miao$Chromosome)

# add tri-nucleotide context annotation
mutations_miao$tri_context_ref <- as.character(Biostrings::getSeq(ref_genome, mutations_miao$Chromosome, mutations_miao$Start_Position - 1, mutations_miao$Start_Position + 1))


for (j in 1:nrow(mutations_miao)){
  print(j)
  if (mutations_miao$Reference_Allele[j] %in% c("G", "A")){
    mutations_miao$tri_context_96[j] <- as.character(reverseComplement(DNAString(mutations_miao$tri_context_ref[j])))
    mutations_miao$change[j] <- paste0(as.character(complement(DNAString(mutations_miao$Reference_Allele[j]))), ">", as.character(complement(DNAString(mutations_miao$Tumor_Seq_Allele2[j]))))
  } else {
    mutations_miao$tri_context_96[j] <- mutations_miao$tri_context_ref[j]
    mutations_miao$change[j] <- paste0(mutations_miao$Reference_Allele[j], ">", mutations_miao$Tumor_Seq_Allele2[j])
  }

  mutations_miao$tri_context_96[j] <- paste0(substr(mutations_miao$tri_context_96[j],1,1), "[", mutations_miao$change[j], "]", substr(mutations_miao$tri_context_96[j],3,3))

}



# count and add number of sbs1 mutations to metadata

for (i in 1:nrow(meta_miao)){
  sample <- meta_miao$SAMPLE_ID[i]
  tot_mut_no <- nrow(mutations_miao[mutations_miao$Tumor_Sample_Barcode == sample,])


  muts_df_sbs1 <- mutations_miao[grepl('\\w\\[C>T\\]G', mutations_miao$tri_context_96) & mutations_miao$Tumor_Sample_Barcode == sample,]
  muts_df_sbs1 <- muts_df_sbs1[substr(muts_df_sbs1$tri_context_96, 1,1) !="T",]

  sbs1_mut_no <- nrow(muts_df_sbs1)


  meta_miao$tot_count[i] <- tot_mut_no
  meta_miao$sbs1_count[i] <- sbs1_mut_no


}

# read in patient meta to obtain age info and add it to sample meta
meta_kid_dfci_patient <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/kidney_exome/kidney_met_whole_exome_science2019/ccrcc_dfci_2019/data_clinical_patient.txt"), sep = "\t", stringsAsFactors = F, header = T, skip = 4)

meta_miao <- merge(meta_miao, meta_kid_dfci_patient[,c("PATIENT_ID", "AGE")], by = "PATIENT_ID")

# sbs1 counts per year
meta_miao$sbs1_per_year <- meta_miao$sbs1_count/meta_miao$AGE

# save the data
write.table(meta_miao, file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/kidney_exome/results/r-objects/meta_kid_dfci_with_counts.tsv"), sep = "\t", row.names = F, quote = F)






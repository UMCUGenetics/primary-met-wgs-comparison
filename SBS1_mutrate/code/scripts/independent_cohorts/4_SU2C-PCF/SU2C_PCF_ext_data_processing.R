
### In this script we process SU2C-PCF WES external data set 



## Init --------------------------------

expanding_path <- '/home/ali313/Documents/studies/master/umc-project'

path_prefix <- c(
  hpc  ='/hpc/',
  local = paste0(expanding_path, '/hpc/')
)
path_prefix <- path_prefix[min(which(dir.exists(path_prefix)))]
path_prefix <- sub('/hpc/','',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')



# read in sample metadata of SU2C-PCF external data set (prostate wes metastatic)


meta_SU2C_PCF <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/prostate_exome/prostate_met_whole_exome_pnas2019/metadata/data_clinical_sample.txt"), sep = "\t", stringsAsFactors = F, header = T, skip = 4) 



# load required packages
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
ref_genome <- getBSgenome(ref_genome)


# read in the mutation file
mutations_pros_dfci <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/prostate_exome/prostate_met_whole_exome_pnas2019/mutations/data_mutations.txt"), sep = "\t", stringsAsFactors = F, header = T)

# subset for sbs
mutations_pros_dfci <- mutations_pros_dfci[nchar(mutations_pros_dfci$Reference_Allele) == 1 & nchar(mutations_pros_dfci$Tumor_Seq_Allele2) == 1,]

# format the data
mutations_pros_dfci <- mutations_pros_dfci[substr(mutations_pros_dfci$Chromosome, 1, 1) != "G",]
mutations_pros_dfci$Chromosome <- paste0("chr", mutations_pros_dfci$Chromosome)

# add tri-nucleotide context annotation
mutations_pros_dfci$tri_context_ref <- as.character(Biostrings::getSeq(ref_genome, mutations_pros_dfci$Chromosome, mutations_pros_dfci$Start_Position - 1, mutations_pros_dfci$Start_Position + 1))


for (j in 1:nrow(mutations_pros_dfci)){
  print(j)
  if (mutations_pros_dfci$Reference_Allele[j] %in% c("G", "A")){
    mutations_pros_dfci$tri_context_96[j] <- as.character(reverseComplement(DNAString(mutations_pros_dfci$tri_context_ref[j])))
    mutations_pros_dfci$change[j] <- paste0(as.character(complement(DNAString(mutations_pros_dfci$Reference_Allele[j]))), ">", as.character(complement(DNAString(mutations_pros_dfci$Tumor_Seq_Allele2[j]))))
  } else {
    mutations_pros_dfci$tri_context_96[j] <- mutations_pros_dfci$tri_context_ref[j]
    mutations_pros_dfci$change[j] <- paste0(mutations_pros_dfci$Reference_Allele[j], ">", mutations_pros_dfci$Tumor_Seq_Allele2[j])
  }

  mutations_pros_dfci$tri_context_96[j] <- paste0(substr(mutations_pros_dfci$tri_context_96[j],1,1), "[", mutations_pros_dfci$change[j], "]", substr(mutations_pros_dfci$tri_context_96[j],3,3))

}


# count and add number of sbs1 mutations to metadata

for (i in 1:nrow(meta_SU2C_PCF)){
  print(i)
  sample <- meta_SU2C_PCF$SAMPLE_ID[i]

  tot_mut_no <- nrow(mutations_pros_dfci[mutations_pros_dfci$Tumor_Sample_Barcode == sample,])


  muts_df_sbs1 <- mutations_pros_dfci[grepl('\\w\\[C>T\\]G', mutations_pros_dfci$tri_context_96) & mutations_pros_dfci$Tumor_Sample_Barcode == sample,]
  muts_df_sbs1 <- muts_df_sbs1[substr(muts_df_sbs1$tri_context_96, 1,1) !="T",]

  sbs1_mut_no <- nrow(muts_df_sbs1)


  meta_SU2C_PCF$tot_count[i] <- tot_mut_no
  meta_SU2C_PCF$sbs1_count[i] <- sbs1_mut_no


}

# read in patient meta to obtain age info and add it to sample meta
meta_pros_dfci_patient <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/prostate_exome/prostate_met_whole_exome_pnas2019/metadata/data_clinical_patient.txt"), sep = "\t", stringsAsFactors = F, header = T, skip = 4)

meta_SU2C_PCF <- merge(meta_SU2C_PCF, meta_pros_dfci_patient[,c("PATIENT_ID", "AGE_AT_DIAGNOSIS")], by = "PATIENT_ID")

# sbs1 counts per year
meta_SU2C_PCF$sbs1_per_year <- meta_SU2C_PCF$sbs1_count/meta_SU2C_PCF$AGE_AT_DIAGNOSIS


# save the data
write.table(meta_SU2C_PCF, file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/prostate_exome/results/r-objects/meta_pros_dfci_with_counts.tsv"), sep = "\t", row.names = F, quote = F)




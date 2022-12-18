
### In this script we process David Quigley et al. WES external data set 



## Init --------------------------------

expanding_path <- '/home/ali313/Documents/studies/master/umc-project'

path_prefix <- c(
  hpc  ='/hpc/',
  local = paste0(expanding_path, '/hpc/')
)
path_prefix <- path_prefix[min(which(dir.exists(path_prefix)))]
path_prefix <- sub('/hpc/','',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')


# load required libraries


library(readxl)
library(stringr)
library(dplyr)
library(tidyverse)
library(stringr)
library(vcfR)
library(BSgenome)
library(Biostrings)





# read in sample metadata of David Quigley et al. external data set (prostate wgs metastatic)


meta_quigley <- as.data.frame(read_xlsx(path = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/prostate_meta_wgs_cell2018/mmc1.xlsx"), sheet = 1))





# the raw data of this cohort was ran by hartwig pipeline to normalize for biases. Below is a conversion table to associate original IDs with hartwig produced IDs 

transfer_table <- read.csv(file = paste0(base_dir, "external_datasets/prostate_meta_wgs_cell2018/transfer-table.tsv"), sep = "\t", stringsAsFactors = F, header = F)
colnames(transfer_table) <- c("IDENTIFIER", "type", "complete_id", "number_id", "type_number_id")
transfer_table <- transfer_table[transfer_table$type == "tumor",]
transfer_table <- transfer_table[-c(3,4),]
row.names(transfer_table) <- 1:nrow(transfer_table)

# annotate the meta with number of total mutations + SBS1 mutations

ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
ref_genome <- getBSgenome(ref_genome)

for (i in 1:nrow(meta_quigley)){
  print(paste0("===============================",i))
  
  # assign path to VCF
  path_to_vcf <- paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/prostate_meta_wgs_cell2018/VCF_Hartwig/p", transfer_table$number[i], "ref-p", transfer_table$number[i], "tumor-bamfilebucket", "/sage_somatic/p", transfer_table$number[i],"tumor.sage.somatic.filtered.vcf.gz")

  # read in vcf file
  vcf <- tryCatch(read.vcfR(path_to_vcf), error=function(e) NULL)

  
  if (!is.null(vcf)){

    df <- vcfR2tidy(vcf)
    muts_df <- data.frame(df$fix)
    muts_df <- muts_df[muts_df$FILTER == "PASS",]
    row.names(muts_df) <- 1:nrow(muts_df)

    # add tri-nucleotide context annotation
    muts_df$tri_context_ref <- as.character(Biostrings::getSeq(ref_genome, muts_df$CHROM, muts_df$POS - 1, muts_df$POS + 1))

      for (j in 1:nrow(muts_df)){
      print(j)
      if (muts_df$REF[j] %in% c("G", "A")){
        muts_df$tri_context_96[j] <- as.character(reverseComplement(DNAString(muts_df$tri_context_ref[j])))
        muts_df$change[j] <- paste0(as.character(complement(DNAString(muts_df$REF[j]))), ">", as.character(complement(DNAString(muts_df$ALT[j]))))
      } else {
        muts_df$tri_context_96[j] <- muts_df$tri_context_ref[j]
        muts_df$change[j] <- paste0(muts_df$REF[j], ">", muts_df$ALT[j])
      }

      muts_df$tri_context_96[j] <- paste0(substr(muts_df$tri_context_96[j],1,1), "[", muts_df$change[j], "]", substr(muts_df$tri_context_96[j],3,3))

    }

    # save the data
    write.table(muts_df, file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/prostate_meta_wgs_cell2018/VCF_Hartwig/annotated/", transfer_table$IDENTIFIER[i], ".sage.somatic.filtered.vcf.gz"), row.names = F, quote = F, sep = "\t")
  }

}



# annotate the meta with number of total mutations + SBS1 mutations



meta_quigley <- meta_quigley[meta_quigley$IDENTIFIER %in% transfer_table$IDENTIFIER,]
colnames(meta_quigley)[(ncol(meta_quigley)-1):ncol(meta_quigley)] <- c("tot_count_original", "sbs1_count_original")


for (i in 1:nrow(meta_quigley)){
  print(i)
  muts_df <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/prostate_meta_wgs_cell2018/VCF_Hartwig/annotated/", meta_quigley$IDENTIFIER[i], ".sage.somatic.filtered.vcf.gz"), header = T, sep = "\t", stringsAsFactors = F)

  tot_mut_no <- nrow(muts_df)

  # sbs1 context
  muts_df_sbs1 <- muts_df[grepl('\\w\\[C>T\\]G', muts_df$tri_context_96),]
  muts_df_sbs1 <- muts_df_sbs1[substr(muts_df_sbs1$tri_context_96, 1,1) !="T",]

  sbs1_mut_no <- nrow(muts_df_sbs1)

  meta_quigley$tot_count_hartwig[i] <- tot_mut_no
  meta_quigley$sbs1_count_hartwig[i] <- sbs1_mut_no


}


# save the data
write.table(meta_quigley, file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/prostate_meta_wgs_cell2018/results/r-objects/meta_quigley_hartwig_counts.tsv"), sep = "\t", row.names = F, quote = F)










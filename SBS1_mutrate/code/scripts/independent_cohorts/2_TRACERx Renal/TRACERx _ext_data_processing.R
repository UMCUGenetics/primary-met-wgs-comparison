
### In this script we process David TRACERx Renal WES external data set 



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
library(Biostrings)
library(BSgenome.Hsapiens.UCSC.hg19)

quals <- c("sage", "QUAL_50", "QUAL_75", "QUAL_100", "QUAL_120")
quals_names <- c(".sage", "_QUAL_50", "_QUAL_75", "_QUAL_100", "_QUAL_120")

sample_ids <- c("EV001_M1", "EV001_M2a", "EV001_M2b", "EV001_R2", "EV001_R3", "EV001_R4", "EV001_R5",
                "EV001_R8", 'EV001_R9', "EV002_M", "EV002_R1", "EV002_R3", "EV002_R4",
                "EV002_R6", "EV002_R7", "EV002_R9")




# annotate the meta with number of total mutations + SBS1 mutations



ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"
ref_genome <- getBSgenome(ref_genome)


# for different quality levels

for (j in 1:5){
  qual <- quals[j]
  print(qual)
  quals_name <- quals_names[j]
  
  #for each sample
  for (i in 1:16){ 
    print(paste0("===============================",i))

    # assign vcf path
    path_to_vcf <- paste0(wd, "external_datasets/renal_TRACERX/aspera/VCFs/", qual, "/", sample_ids[i], quals_name, ".vcf.gz")

    # read in vcf file
    vcf <- read.vcfR(path_to_vcf)


    df <- vcfR2tidy(vcf)
    muts_df <- data.frame(df$fix)
    muts_df <- muts_df[muts_df$FILTER == "PASS",]

    muts_df$CHROM <- paste0("chr", muts_df$CHROM)

    # add tri-nucleotide context annotation
    muts_df$tri_context_ref <- as.character(Biostrings::getSeq(ref_genome, muts_df$CHROM, muts_df$POS - 1, muts_df$POS + 1))

    for (k in 1:nrow(muts_df)){
      print(k)
      if (muts_df$REF[k] %in% c("G", "A")){
        muts_df$tri_context_96[k] <- as.character(reverseComplement(DNAString(muts_df$tri_context_ref[k])))
        muts_df$change[k] <- paste0(as.character(complement(DNAString(muts_df$REF[k]))), ">", as.character(complement(DNAString(muts_df$ALT[k]))))
      } else {
        muts_df$tri_context_96[k] <- muts_df$tri_context_ref[k]
        muts_df$change[k] <- paste0(muts_df$REF[k], ">", muts_df$ALT[k])
      }

      muts_df$tri_context_96[k] <- paste0(substr(muts_df$tri_context_96[k],1,1), "[", muts_df$change[k], "]", substr(muts_df$tri_context_96[k],3,3))

    }

    muts_df[nchar(muts_df$REF) > 1 | nchar(muts_df$ALT) > 1 ,c("tri_context_ref", "tri_context_96", "change")] <- c("NA", "NA", "NA")
    row.names(muts_df) <- 1:nrow(muts_df)

    # save the data
    write.table(muts_df, file = paste0(wd, "external_datasets/renal_TRACERX/processed/", qual, "/", sample_ids[i], quals_name, ".txt"), row.names = F, quote = F, sep = "\t")

  }

}


# summarize and save number of total/sbs1 mutations


df <- data.frame(sample_name = character(), quality = character(), tot_mut = numeric(), sbs1_mut = numeric())

for (j in 1:5){
  
  qual <- quals[j]
  quals_name <- quals_names[j]
  
  for (i in 1:16){
    sample_id <- sample_ids[i]

    print(paste0("########### ", sample_ids[i]))
    muts_df <- read.csv(file = paste0(wd, "external_datasets/renal_TRACERX/processed/", qual, "/", sample_ids[i], quals_name, ".txt"), header = T, sep = "\t", stringsAsFactors = F)

    tot_mut_no <- nrow(muts_df)

    muts_df_sbs1 <- muts_df[grepl('\\w\\[C>T\\]G', muts_df$tri_context_96),]
    muts_df_sbs1 <- muts_df_sbs1[substr(muts_df_sbs1$tri_context_96, 1,1) !="T",]

    sbs1_mut_no <- nrow(muts_df_sbs1)

    df[i + (j-1)*16,] <- c(sample_id, quals_name, tot_mut_no, sbs1_mut_no)

  }

}

# save the data
write.table(df, file = paste0(wd, "external_datasets/renal_TRACERX/processed/summary_matrix.txt"), row.names = F, quote = F, sep = "\t")

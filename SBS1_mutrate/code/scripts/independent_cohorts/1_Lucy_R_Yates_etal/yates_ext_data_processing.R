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

source(paste0(base_dir, 'drivers/analysis/dna-rep-ann/final-update/code_for_github/SBS1/independent_cohorts/1_Lucy_R_Yates_etal/yates_source.R'))



## Summarize total counts  --------------------------------


row_nr <- length(patsam_ids)

summary_df_tot <- data.frame(patient_id = character(row_nr), total = numeric(row_nr), a = numeric(row_nr), b = numeric(row_nr), c = numeric(row_nr), d = numeric(row_nr), e = numeric(row_nr),
                             a2 = numeric(row_nr), b2 = numeric(row_nr), c2 = numeric(row_nr), d2 = numeric(row_nr), e2 = numeric(row_nr))

summary_df_tot[summary_df_tot == 0] <- NA

for (i in 1:row_nr){
  print(i)
  
  muts_df <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/lyates_cancer_cell_2014/WGS_mutation_calls/SUPP_DATA_5/substitution_and_subclonality_data/", patient_ids[i], "_dirichlet_subs_plus_timing_v0.1.txt"), header = T, stringsAsFactors = F, sep = "\t")
  
  
  # in the tenth sample metadata there are "INCONCLUSIVE" annotations that needs to be removed
  
  if (i == 10){
    muts_df_simplified <- muts_df[,paste0("VAR_CATEGORY_", names(patsam_ids)[i], patsam_ids[[i]])]
    new.muts_df_simplified <- cbind(muts_df_simplified, muts_df_simplified%>%
                                      unite(together, colnames(muts_df_simplified), sep=""))
    muts_df_simplified <- new.muts_df_simplified %>% filter(!grepl('INCONCLUSIVE', together))
    muts_df_simplified <- muts_df_simplified[,-ncol(muts_df_simplified)]
  } else {
    muts_df_simplified <- muts_df[,paste0("dis_CLASS_", names(patsam_ids)[i], patsam_ids[[i]])]
  }
  
  
  muts_df_simplified[muts_df_simplified == "PRESENT"] <- 1
  muts_df_simplified[muts_df_simplified == "ABSENT"] <- 0
  
  
  muts_df_simplified <- muts_df_simplified[,sapply(str_split(colnames(muts_df_simplified), pattern = "_"), "[[", 3) %in% meta_lyates$sample]
  
  
  muts_df_simplified <- sapply(muts_df_simplified, as.numeric)
  
  
  # count number of mutations that are present in each sample
  
  counts <- colSums(muts_df_simplified)
  
  names(counts) <- str_sub(colnames(muts_df_simplified), start= -1)
  
  # count only unique mutations (exclusive mutations)
  muts_df_simplified_exclusives <- muts_df_simplified[rowSums(muts_df_simplified) == 1,]
  
  counts_exclusives <- colSums(muts_df_simplified_exclusives)
  
  names(counts_exclusives) <- paste0(str_sub(colnames(muts_df_simplified_exclusives), start= -1), "2")
  
  summary_df_tot[i, names(counts)] <- as.numeric(counts)
  summary_df_tot[i, names(counts_exclusives)] <- as.numeric(counts_exclusives)
  summary_df_tot[i, "patient_id"] <- patient_ids[i]
  summary_df_tot[i, "total"] <- nrow(muts_df_simplified)
  
}



# add the total summary to metadata

meta_lyates$tot_sample_count <- NA
meta_lyates$tot_sample_count_excl <- NA
meta_lyates$tot_patient_count <- NA

for (i in 1:nrow(summary_df_tot)){
  
  # add number of total sample/patient mutation counts
  vect <- summary_df_tot[i,3:7]
  vec <- as.numeric(vect)
  
  names(vec) <- colnames(vect)
  names(vec) <- paste0(summary_df_tot$patient_id[i], colnames(summary_df_tot)[3:7])
  
  vec <- vec[names(vec) %in% meta_lyates$sample[meta_lyates$sample %in% names(vec)]]
  vec <- vec[order(match(names(vec), meta_lyates$sample[meta_lyates$sample %in% names(vec)]))]
  
  meta_lyates$tot_sample_count[meta_lyates$sample %in% names(vec)] <- as.numeric(vec)
  meta_lyates$tot_patient_count[meta_lyates$sample %in% names(vec)] <- summary_df_tot[i,2]
  
  
  # add number of total exclusive mutation count
  vect <- summary_df_tot[i,8:12]
  vec <- as.numeric(vect)
  
  names(vec) <- colnames(vect)
  names(vec) <- paste0(summary_df_tot$patient_id[i], colnames(summary_df_tot)[3:7])
  
  vec <- vec[names(vec) %in% meta_lyates$sample[meta_lyates$sample %in% names(vec)]]
  vec <- vec[order(match(names(vec), meta_lyates$sample[meta_lyates$sample %in% names(vec)]))]
  
  meta_lyates$tot_sample_count_excl[meta_lyates$sample %in% names(vec)] <- as.numeric(vec)
  
  
}


## Add tri-nucleotide annotation --------------------------------

# read required packages
ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"

library(Biostrings)
library(ref_genome, character.only = TRUE)
library(BSgenome.Hsapiens.UCSC.hg19)

ref_genome <- getBSgenome(ref_genome)


for (i in 1:length(patsam_ids)){
  print(paste0("===============================",i))
  
  # read in data
  muts_df <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/lyates_cancer_cell_2014/WGS_mutation_calls/SUPP_DATA_5/substitution_and_subclonality_data/", patient_ids[i], "_dirichlet_subs_plus_timing_v0.1.txt"), header = T, stringsAsFactors = F, sep = "\t")
  
  muts_df$Chrom <- paste0("chr", muts_df$Chrom)
  muts_df$tri_context_ref <- as.character(Biostrings::getSeq(ref_genome, muts_df$Chrom, muts_df$Pos - 1, muts_df$Pos + 1))
  
  # add tri-nucleotide context
  for (j in 1:nrow(muts_df)){
    if (muts_df$Ref[j] %in% c("G", "A")){
      muts_df$tri_context_96[j] <- as.character(reverseComplement(DNAString(muts_df$tri_context_ref[j])))
    } else {
      muts_df$tri_context_96[j] <- muts_df$tri_context_ref[j]
    }
    muts_df$tri_context_96[j] <- paste0(substr(muts_df$tri_context_96[j],1,1), "[", muts_df$change[j], "]", substr(muts_df$tri_context_96[j],3,3))
  }
  
  # save the annotated data set
  write.table(muts_df, file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/lyates_cancer_cell_2014/WGS_mutation_calls/SUPP_DATA_5/substitution_and_subclonality_data/", patient_ids[i], "_dirichlet_subs_plus_timing_v0.2.txt"), row.names = F, quote = F, sep = "\t")

}



# summarize sbs1 counts

summary_df_sbs1 <- data.frame(patient_id = character(row_nr), total = numeric(row_nr), a = numeric(row_nr), b = numeric(row_nr), c = numeric(row_nr), d = numeric(row_nr), e = numeric(row_nr),
                             a2 = numeric(row_nr), b2 = numeric(row_nr), c2 = numeric(row_nr), d2 = numeric(row_nr), e2 = numeric(row_nr))
summary_df_sbs1[summary_df_sbs1 == 0] <- NA


for (i in 1:length(patsam_ids)){
  
  # read in data
  muts_df <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/lyates_cancer_cell_2014/WGS_mutation_calls/SUPP_DATA_5/substitution_and_subclonality_data/", patient_ids[i], "_dirichlet_subs_plus_timing_v0.2.txt"), header = T, stringsAsFactors = F, sep = "\t")
  
  # SBS1 context
  muts_df <- muts_df[grepl('\\w\\[C>T\\]G', muts_df$tri_context_96),]
  muts_df <- muts_df[substr(muts_df$tri_context_96, 1,1) !="T",]


  if (i == 10){
    muts_df_simplified <- muts_df[,paste0("VAR_CATEGORY_", names(patsam_ids)[i], patsam_ids[[i]])]
    new.muts_df_simplified <- cbind(muts_df_simplified, muts_df_simplified%>%
                                      unite(together, colnames(muts_df_simplified), sep=""))
    muts_df_simplified <- new.muts_df_simplified %>% filter(!grepl('INCONCLUSIVE', together))
    muts_df_simplified <- muts_df_simplified[,-ncol(muts_df_simplified)]
  } else {
    muts_df_simplified <- muts_df[,paste0("dis_CLASS_", names(patsam_ids)[i], patsam_ids[[i]])]
  }

  muts_df_simplified[muts_df_simplified == "PRESENT"] <- 1
  muts_df_simplified[muts_df_simplified == "ABSENT"] <- 0

  muts_df_simplified <- sapply(muts_df_simplified, as.numeric)


  muts_df_simplified <- muts_df_simplified[,sapply(str_split(colnames(muts_df_simplified), pattern = "_"), "[[", 3) %in% meta_lyates$sample]


  counts <- colSums(muts_df_simplified)

  names(counts) <- str_sub(colnames(muts_df_simplified), start= -1)


  muts_df_simplified2 <- muts_df_simplified[rowSums(muts_df_simplified) == 1,]

  counts2 <- colSums(muts_df_simplified2)

  names(counts2) <- paste0(str_sub(colnames(muts_df_simplified2), start= -1), "2")


  summary_df_sbs1[i, names(counts)] <- as.numeric(counts)
  summary_df_sbs1[i, names(counts2)] <- as.numeric(counts2)
  summary_df_sbs1[i, "patient_id"] <- patient_ids[i]
  summary_df_sbs1[i, "total"] <- nrow(muts_df_simplified)


}


# add the summary of SBS1 to metadata

meta_lyates$sbs1_sample_count <- NA
meta_lyates$sbs1_sample_count_excl <- NA
meta_lyates$sbs1_patient_count <- NA

for (i in 1:nrow(summary_df_sbs1)){

  # add number of sample/patient SBS1 mutation counts
  vect <- summary_df_sbs1[i,3:7]
  vec <- as.numeric(vect)
  
  names(vec) <- colnames(vect)
  names(vec) <- paste0(summary_df_sbs1$patient_id[i], colnames(summary_df_sbs1)[3:7])

  vec <- vec[names(vec) %in%meta_lyates$sample[meta_lyates$sample %in% names(vec)]]
  vec <- vec[order(match(names(vec), meta_lyates$sample[meta_lyates$sample %in% names(vec)]))]

  meta_lyates$sbs1_sample_count[meta_lyates$sample %in% names(vec)] <- as.numeric(vec)
  meta_lyates$sbs1_patient_count[meta_lyates$sample %in% names(vec)] <- summary_df_sbs1[i,2]


  # add number of exclusive SBS1 mutation count
  vect <- summary_df_sbs1[i,8:12]
  vec <- as.numeric(vect)

  names(vec) <- colnames(vect)
  names(vec) <- paste0(summary_df_sbs1$patient_id[i], colnames(summary_df_sbs1)[3:7])

  vec <- vec[names(vec) %in% meta_lyates$sample[meta_lyates$sample %in% names(vec)]]
  vec <- vec[order(match(names(vec), meta_lyates$sample[meta_lyates$sample %in% names(vec)]))]

  meta_lyates$sbs1_sample_count_excl[meta_lyates$sample %in% names(vec)] <- as.numeric(vec)

}




# format the annotated meta data
row.names(meta_lyates) <- 1:nrow(meta_lyates)
meta_lyates$days_to_biopsy <- as.numeric(meta_lyates$days_to_biopsy)
meta_lyates$Age_at_diagnosis <- as.numeric(meta_lyates$Age_at_diagnosis)
meta_lyates$age_at_biopsy <- meta_lyates$Age_at_diagnosis + round(meta_lyates$days_to_biopsy/365, digit = 1)



# save the data
write.table(meta_lyates, file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external_datasets/lyates_cancer_cell_2014/results/r-objects/meta_lyates_cpg_counts.tsv"), sep = "\t", row.names = F, quote = F)



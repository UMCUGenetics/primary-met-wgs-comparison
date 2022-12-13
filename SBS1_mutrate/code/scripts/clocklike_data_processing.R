
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




# Read in and process raw PCWG age information


pcawg_meta <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/PCAWG/Metadata/donor.info.followup.all_projects.csv"), header = T, stringsAsFactors = F, sep = "\t")


colnames(pcawg_meta)[c(1,9)] <- c("sample_id", "age")
pcawg_meta <- pcawg_meta[,c(1,9)]




# Merging the age info of both cohorts to the sample metadata object
age_meta_all <- rbind(hartwig_meta, pcawg_meta)
sample_metadata <- merge(sample_metadata, age_meta_all, by = "sample_id")



## Read in mutation counts--------------------------------

# read in SBS1 timing counts and merge with relevant columns from sample metadata to make a new data frame dubbed all_signatures
sbs1_timing <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/rebuttal/sbs1_timing.tsv.gz"), header = T, stringsAsFactors = F, sep = "\t")
all_signatures <- merge(sbs1_timing, sample_metadata[,c("sample_id", "age", "tmb_total", "clonal_tmb", "subclonal_tmb")], all.x = T)



# read in SBS5 and 40 timing counts and add it to all_signatures data frame
all_timing_sbs5 <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/rebuttal/sbs5_timing.tsv.gz"), sep = "\t", stringsAsFactors = F, header = T)
all_timing_sbs40 <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/rebuttal/sbs40_timing.tsv.gz"), sep = "\t", stringsAsFactors = F, header = T)

# NA to 0 for numeric calculations
all_timing_sbs5[all_timing_sbs5$sbs5_count == 0,3:7] <- 0
all_timing_sbs40[all_timing_sbs40$sbs40_count == 0,3:7] <- 0

# SBS5_40 = SBS5 + SBS40
all_timing_sbs5_40 <-all_timing_sbs5[,2:7] + all_timing_sbs40[,2:7]
all_timing_sbs5_40 <- cbind(all_timing_sbs5[,1], all_timing_sbs5_40)
colnames(all_timing_sbs5_40) <- c("sample_id", "sbs5_40_count", "sbs5_40_subcl_subclonal", "sbs5_40_clonal_late", "sbs5_40_clonal_early", "sbs5_40_clonal_na", "sbs5_40_mutationTimeR_subclonal")


# merge the data with all_signatures
all_signatures <- merge(merge(merge(all_signatures, all_timing_sbs5, by = "sample_id"), all_timing_sbs40, by = "sample_id"), all_timing_sbs5_40, by = "sample_id")


# read in total mutation timing counts and add it to all_signatures data frame
all_timing <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/rebuttal/all_timing.tsv.gz"), sep = "\t", stringsAsFactors = F, header = T)
all_signatures <- merge(all_signatures, all_timing[,c("sample_id", "tmb_clonal_late", "tmb_clonal_early", "tmb_clonal_na", "tmb_mutationTimeR_subclonal")], by = "sample_id", all.x = T)



# read in signature contribution matrix and add merge it with all_signatures data frame
cont_mat <- read.csv(file = paste0(base_dir, "passengers/processed/sigs_denovo/extractions/12_fixed_seeds/sig_contrib/fit_lsq.post_processed/denovo_contribs.lsq.post_processed.txt.gz"), header = T, stringsAsFactors = F, sep = "\t")
cont_mat$sample_id <- row.names(cont_mat)

# SBS5_40 = SBS5 + SBS40
cont_mat$SBS5_40 <- cont_mat$SBS5 + cont_mat$SBS40

cont_mat <- cont_mat[,c("sample_id", "SBS5", "SBS40", "SBS5_40")]
colnames(cont_mat)[2:4] <- c("sbs5_exposure", "sbs40_exposure", "sbs5_40_exposure") 

# merge the data with all_signatures
all_signatures <- merge(all_signatures, cont_mat, by = "sample_id")


# read in ploidy info and add it to all_signatures data frame
ploidy_prim <- as.data.frame(read_xlsx(path = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external-files/SuppTable2_karyotype _karyotype_271022.xlsx"), sheet = 1))
ploidy_met <- as.data.frame(read_xlsx(path = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external-files/SuppTable2_karyotype _karyotype_271022.xlsx"), sheet = 2))

ploidy <- rbind(ploidy_prim, ploidy_met)
colnames(ploidy)[1] <- "sample_id_2"


ploidy <- merge(ploidy, sample_metadata[,c("sample_id_2", "sample_id")])
colnames(ploidy)[2] <- "ploidy"

# merge the data with all_signatures
all_signatures <- merge(all_signatures, ploidy[,c("sample_id", "ploidy")], by = "sample_id")




## Compute clonality ratios --------------------------------



# compute the total clonality ratios
all_signatures$total_clonality_ratio <- all_signatures$clonal_tmb/all_signatures$tmb_total
all_signatures$total_subclonality_ratio <- all_signatures$subclonal_tmb/all_signatures$tmb_total


# compute the sbs1 clonality ratios
all_signatures$sbs1_clonality_ratio <- (all_signatures$sbs1_count - all_signatures$sbs1_subcl_subclonal)/all_signatures$sbs1_count
all_signatures$sbs1_normalized_clonality <- all_signatures$sbs1_clonality_ratio/all_signatures$total_clonality_ratio


# compute the sbs5 clonality ratios
all_signatures$sbs5_clonality_ratio <- (all_signatures$sbs5_count - all_signatures$sbs5_subcl_subclonal)/all_signatures$sbs5_count
all_signatures$sbs5_normalized_clonality <- all_signatures$sbs5_clonality_ratio/all_signatures$total_clonality_ratio

# compute the sbs40 clonality ratios
all_signatures$sbs40_clonality_ratio <- (all_signatures$sbs40_count - all_signatures$sbs40_subcl_subclonal)/all_signatures$sbs40_count
all_signatures$sbs40_normalized_clonality <- all_signatures$sbs40_clonality_ratio/all_signatures$total_clonality_ratio


# compute the sbs5_40 clonality ratios
all_signatures$sbs5_40_clonality_ratio <- (all_signatures$sbs5_40_count - all_signatures$sbs5_40_subcl_subclonal)/all_signatures$sbs5_40_count
all_signatures$sbs5_40_normalized_clonality <- all_signatures$sbs5_40_clonality_ratio/all_signatures$total_clonality_ratio



## Compute timing ratios --------------------------------



# compute the sbs1 timing ratios
all_signatures$sbs1_timing_normalized1 <- (all_signatures$sbs1_clonal_late/(all_signatures$sbs1_clonal_late + all_signatures$sbs1_clonal_early))/(all_signatures$tmb_clonal_late/(all_signatures$tmb_clonal_late + all_signatures$tmb_clonal_early))
all_signatures$sbs1_timing_normalized2 <- (all_signatures$sbs1_clonal_late/(all_signatures$sbs1_clonal_late + all_signatures$sbs1_clonal_early + all_signatures$sbs1_clonal_na))/(all_signatures$tmb_clonal_late/(all_signatures$tmb_clonal_late + all_signatures$tmb_clonal_early + all_signatures$tmb_clonal_na))

# compute the sbs5 timing ratios
all_signatures$sbs5_timing_normalized1 <- (all_signatures$sbs5_clonal_late/(all_signatures$sbs5_clonal_late + all_signatures$sbs5_clonal_early))/(all_signatures$tmb_clonal_late/(all_signatures$tmb_clonal_late + all_signatures$tmb_clonal_early))
all_signatures$sbs5_timing_normalized2 <- (all_signatures$sbs5_clonal_late/(all_signatures$sbs5_clonal_late + all_signatures$sbs5_clonal_early + all_signatures$sbs5_clonal_na))/(all_signatures$tmb_clonal_late/(all_signatures$tmb_clonal_late + all_signatures$tmb_clonal_early + all_signatures$tmb_clonal_na))

# compute the sbs40 timing ratios
all_signatures$sbs40_timing_normalized1 <- (all_signatures$sbs40_clonal_late/(all_signatures$sbs40_clonal_late + all_signatures$sbs40_clonal_early))/(all_signatures$tmb_clonal_late/(all_signatures$tmb_clonal_late + all_signatures$tmb_clonal_early))
all_signatures$sbs40_timing_normalized2 <- (all_signatures$sbs40_clonal_late/(all_signatures$sbs40_clonal_late + all_signatures$sbs40_clonal_early + all_signatures$sbs40_clonal_na))/(all_signatures$tmb_clonal_late/(all_signatures$tmb_clonal_late + all_signatures$tmb_clonal_early + all_signatures$tmb_clonal_na))

# compute the sbs5_40 timing ratios
all_signatures$sbs5_40_timing_normalized1 <- (all_signatures$sbs5_40_clonal_late/(all_signatures$sbs5_40_clonal_late + all_signatures$sbs5_40_clonal_early))/(all_signatures$tmb_clonal_late/(all_signatures$tmb_clonal_late + all_signatures$tmb_clonal_early))
all_signatures$sbs5_40_timing_normalized2 <- (all_signatures$sbs5_40_clonal_late/(all_signatures$sbs5_40_clonal_late + all_signatures$sbs5_40_clonal_early + all_signatures$sbs5_40_clonal_na))/(all_signatures$tmb_clonal_late/(all_signatures$tmb_clonal_late + all_signatures$tmb_clonal_early + all_signatures$tmb_clonal_na))







# convert NAN to NA
all_signatures <- all_signatures %>% mutate_all(~ifelse(is.nan(.), NA, .))




# set the right column order

all_signatures <- all_signatures[,c("sample_id", "sbs1_count", "sbs1_subcl_subclonal", "sbs1_clonal_late", "sbs1_clonal_early", "sbs1_clonal_na", "sbs1_mutationTimeR_subclonal",
                          "sbs5_exposure", "sbs5_count", "sbs5_subcl_subclonal", "sbs5_clonal_late", "sbs5_clonal_early", "sbs5_clonal_na", "sbs5_mutationTimeR_subclonal",
                          "sbs40_exposure", "sbs40_count", "sbs40_subcl_subclonal", "sbs40_clonal_late", "sbs40_clonal_early", "sbs40_clonal_na", "sbs40_mutationTimeR_subclonal",
                          "sbs5_40_exposure", "sbs5_40_count", "sbs5_40_subcl_subclonal", "sbs5_40_clonal_late", "sbs5_40_clonal_early", "sbs5_40_clonal_na",
                          "sbs5_40_mutationTimeR_subclonal", "tmb_total", "clonal_tmb", "subclonal_tmb", "tmb_clonal_late", "tmb_clonal_early", "tmb_clonal_na", "tmb_mutationTimeR_subclonal",
                          "age", "ploidy", "total_clonality_ratio", "total_subclonality_ratio", "sbs1_clonality_ratio", 
                          "sbs1_normalized_clonality", "sbs1_timing_normalized1",
                          "sbs1_timing_normalized2", "sbs5_clonality_ratio", "sbs5_normalized_clonality", "sbs40_clonality_ratio", "sbs40_normalized_clonality",
                          "sbs5_40_clonality_ratio", "sbs5_40_normalized_clonality", "sbs5_timing_normalized1", "sbs5_timing_normalized2", "sbs40_timing_normalized1",
                          "sbs40_timing_normalized2", "sbs5_40_timing_normalized1", "sbs5_40_timing_normalized2")]



# save the data
write.table(all_signatures, file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/rebuttal/complete1-5-40.tsv"), sep = "\t", row.names = F, quote = F)




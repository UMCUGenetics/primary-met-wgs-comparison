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
sample_metadata[,"cancer_type"] <- factor(sample_metadata[,"cancer_type"], levels = cancer_type_order_fig1)
sample_metadata[,"cancer_type_code"] <- factor(sample_metadata[,"cancer_type_code"], levels = cancer_type_code_order_fig1)
sample_metadata[,"gender"] <- factor(sample_metadata[,"gender"], levels = c("MALE", "FEMALE"))


## Define cancer types and their respective treatment groups in our dataset  --------------------------------

cancer_type_code_for_active_driver <- c("BLCA", "BRCA", "BRCA_ERpos", "BRCA_HER2pos", "BRCA_TNB", "CESC", "CHOL", "COREAD", 
                                   "COREAD_MSS", "DLBCL", "ESCA", "GBM", "HNSC", "KIRC", "LIHC", "LMS", "LPS", "LUAD",
                                   "LUSC", "OV", "PAAD", "PANET", "PRAD", "SKCM", "STAD", "THCA", "UCEC")

treatments_groups <- list(c("Platinum", "Pyrimidine_antagonist", "untreated"),
                       c("Alkaloid", "Alkylating", "Anthracycline", "Anti_AR__GnRH", "Anti_HER2", "Anti_VEGF",
                         "Aromatase_inhibitor", "CDK4__6_inhibitor", "Folate_antagonist", "Microtubule_inhibitor",
                         "mTOR_inhibitor", "Platinum", "Pyrimidine_antagonist", "Selective_ER_modulator", "Taxane",
                         "untreated"),
                       c("Alkaloid", "Alkylating", "Anthracycline", "Anti_AR__GnRH", "Anti_HER2", "Anti_VEGF",
                         "Aromatase_inhibitor", "CDK4__6_inhibitor", "Folate_antagonist", "Microtubule_inhibitor",
                         "mTOR_inhibitor", "Platinum", "Pyrimidine_antagonist", "Selective_ER_modulator", "Taxane",
                         "untreated"),
                       c("Alkylating", "Anthracycline", "Anti_AR__GnRH", "Anti_HER2", "Aromatase_inhibitor",
                         "Pyrimidine_antagonist", "Selective_ER_modulator", "Taxane", "untreated"),
                       c("Alkylating", "Anthracycline", "Platinum", "Pyrimidine_antagonist", "Taxane", "untreated"),
                       c("Platinum", "Taxane", "untreated"),
                       c("Platinum", "Pyrimidine_antagonist", "untreated"),
                       c("Anti_EGFR", "Anti_VEGF", "Platinum", "Pyrimidine_antagonist", "Topoisomerase_inhibitor", "untreated"),
                       c("Anti_EGFR", "Anti_VEGF", "Platinum", "Pyrimidine_antagonist", "Topoisomerase_inhibitor", "untreated"),
                       c("untreated"),
                       c("Platinum", "Pyrimidine_antagonist", "Taxane", "untreated"),
                       c("Alkylating", "untreated"),
                       c("Platinum", "untreated"),
                       c("Multikinase_inhibitor", "untreated"),
                       c("untreated"),
                       c("Anthracycline", "untreated"),
                       c("untreated"),
                       c("Anti_EGFR", "Folate_antagonist", "Immunotherapy", "Platinum", "untreated"),
                       c("Platinum", "Pyrimidine_antagonist", "untreated"),
                       c("Anthracycline", "Anti_VEGF", "Platinum", "Pyrimidine_antagonist", "Selective_ER_modulator",
                         "Taxane", "untreated"),
                       c("Platinum", "Pyrimidine_antagonist", "Topoisomerase_inhibitor", "untreated"),
                       c("untreated"),
                       c("Anti_AR__GnRH", "Immunotherapy", "Taxane", "untreated"),
                       c("BRAF_inhibitor", "Immunotherapy", "MEK_inhibitor", "untreated"),
                       c("Platinum", "Pyrimidine_antagonist", "untreated"),
                       c("untreated"),
                       c("Platinum", "Taxane", "untreated"))



genomic_elements <- c("3utr", "5utr", "lncrna_exons", "nc_elements", "proximal_promoters", "splice_sites")

names(treatments_groups) <- cancer_type_code_for_active_driver


## prepare the mutation input files (per cancer type per treatment) by modifying the format of mutation data frame and save them

for (i in 1:length(cancer_type_code_for_active_driver)) {
  
  cancer_type <- cancer_type_code_for_active_driver[i]
  associated_treatments <- treatments_groups[[cancer_type]]
  
  for (j in 1:length(associated_treatments)){
    
    treatment <- associated_treatments[j]
    
    # define the mutation file location
    path <- paste0(base_dir, "drivers/processed/resistance/dndscv/",
                   cancer_type, "/", treatment, ".dndscv.input.tsv.gz")
    
    # red in the data
    mutation_df <- tryCatch(read.csv(file = path, stringsAsFactors = F, header = T, sep = "\t"), error=function(e) NULL)

      if (!is.null(mutation_df)) {
      
      # removing redundant rows
      mutation_df <- mutation_df[mutation_df$sampleID != "sampleID",]

      # format the data frame
      mutation_df$chr <- paste0("chr", mutation_df$chr)
      mutation_df[,3] <- as.numeric(mutation_df[,3])
      mutation_df$pos1 <- mutation_df$pos
      mutation_df$pos2 <- mutation_df$pos + nchar(mutation_df$ref) - 1
      mutation_df <- mutation_df[,-3]
      colnames(mutation_df) <- c("patient", "chr", "ref", "alt", "pos1", "pos2")
      mutation_df <- mutation_df[,c(2,5,6,3,4,1)]
      mutation_df <- mutation_df[mutation_df$chr != "chrMT",]
      
      # define the saving location
      save.dir <- paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external-files/input-som-muts-ActiveDriverWGS/version-202211/", cancer_type)
      if (!dir.exists(save.dir)) {
        print("dir not exists!")
        dir.create(path = save.dir, recursive = T, mode = "0777")
        }

      save.path <- paste0(save.dir, "/", cancer_type, ".", treatment, ".dndscv.activedriverwgs.input.tsv.gz")
      
      # saving the modified data frame
      write.table(mutation_df, file = gzfile(save.path), quote = F, row.names = F, sep = "\t")
    } else {print(paste0("ALERT: ", cancer_type, " and ", treatment, " is absent!"))}
  }
}


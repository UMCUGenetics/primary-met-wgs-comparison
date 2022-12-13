
## Init --------------------------------

expanding_path <- '/home/ali313/Documents/studies/master/umc-project'

path_prefix <- c(
  hpc  ='/hpc/',
  local = paste0(expanding_path, '/hpc/')
)
path_prefix <- path_prefix[min(which(dir.exists(path_prefix)))]
path_prefix <- sub('/hpc/','',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')


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



# summarize the significantly mutated non-coding elements

datalist = list()

for (i in 1:length(cancer_type_code_for_active_driver)) {
  cancer_type_of_interest <- cancer_type_code_for_active_driver[i]
  print(i)
  for (j in 1:length(treatments_groups[[cancer_type_of_interest]])) {
    treatment_of_interest <- treatments_groups[[cancer_type_of_interest]][j]
    print(j)
    for (k in 1:length(genomic_elements)) {
      genomic_element_of_interest <- genomic_elements[k]
      print(k)

      
      # read in the output of activeDriverWGS
      file_path <- paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/non-coding-drivers/version-202211/", cancer_type_of_interest, "/", treatment_of_interest, "/",
                          genomic_element_of_interest, "/", cancer_type_of_interest, ".", treatment_of_interest, ".",
                          genomic_element_of_interest, ".results.activedriverwgs.tsv.gz")
      ff <- read.csv(file = file_path,
                     header = T, stringsAsFactors = F, sep = "\t")

      # only keeping elements with an FDR score of 0.05 or below
      ff <- ff[ff$fdr_element <= 0.05,]

      if (nrow(ff) > 0){
        attachment <- data.frame(cacnet_type_code = rep(cancer_type_of_interest, times = nrow(ff)), treatment = rep(treatment_of_interest, times = nrow(ff)),
                                 genomic_element = rep(genomic_element_of_interest, times = nrow(ff)))
        ff <- cbind(ff,attachment)
        datalist[[paste0(as.character(i),as.character(j),as.character(k))]] <- ff
      }
    }
  }
}

# convert list of all significantly mutated non-coding elements to a data frame
big_data = do.call(rbind, datalist)

# save the final data frame
write.table(big_data, file = gzfile(paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/non-coding-drivers/version-202211/active-driver-all-sig-hits-version-202211.txt.gz")), quote = F, row.names = F, sep = "\t")





## Init --------------------------------

expanding_path <- '/home/ali313/Documents/studies/master/umc-project'

path_prefix <- c(
  hpc  ='/hpc/',
  local = paste0(expanding_path, '/hpc/')
)
path_prefix <- path_prefix[min(which(dir.exists(path_prefix)))]
path_prefix <- sub('/hpc/','',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')

# load ActiveDriverWGSR package

devtools::load_all(paste0(base_dir,'drivers/analysis/dna-rep-ann/final-update/external-files/ActiveDriver/ActiveDriverWGSR/'))


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




# prepring activeDriverWGS input (this is script is run by ActiveDriverWGSDataProcessor.sh)
args <- commandArgs(trailingOnly=T)

i <- as.numeric(args[1])
cancer_type_of_interest <- cancer_type_code_for_active_driver[i]
print(i)
print(cancer_type_of_interest)
print("========================================================================")

j <- as.numeric(args[2])
treatment_of_interest <- treatments_groups[[cancer_type_of_interest]][j]
print(j)
print(treatment_of_interest)
print("========================================================================")


k <- as.numeric(args[3])
genomic_element_of_interest <- genomic_elements[k]
print(k)
print(genomic_element_of_interest)
print("========================================================================")


# defining the saving location
save.dir <- paste0(wd, "r-objects/non-coding-drivers/version-202211/", cancer_type_of_interest, "/", treatment_of_interest, "/",
                   genomic_element_of_interest)

if (!dir.exists(save.dir)) {
  print("dir not exists!")
  dir.create(path = save.dir, recursive = T, mode = "0777")
}

save.path <- paste0(save.dir, "/", cancer_type_of_interest, ".", treatment_of_interest, ".", genomic_element_of_interest,
                    ".results.activedriverwgs.tsv.gz")


# read in the mutation input data
mutation_input_path <- paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external-files/input-som-muts-ActiveDriverWGS/version-202211/", cancer_type_of_interest, "/", cancer_type_of_interest, ".", treatment_of_interest, ".dndscv.activedriverwgs.input.tsv.gz")
mutation_input <- read.csv(file = mutation_input_path, stringsAsFactors = F, header = T, sep = "\t")

# make sure there are not any duplicated rows
mutation_input <- mutation_input[!duplicated(mutation_input),]


# read in the mutation inpput data
element_input_path <- paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external-files/genomic-elements/", genomic_element_of_interest , ".tsv.gz")
element_input <- read.csv(file = element_input_path, stringsAsFactors = F, header = T, sep = "\t")

# make sure the rows corresponding to chromosome M are excluded
element_input <- element_input[element_input$chr != "chrM",]



## Run AcriveDriverWGS   --------------------------------

results = ActiveDriverWGS(mutations = mutation_input,
                          elements = element_input) 


# save AcriveDriverWGS result
write.table(results,
            file = gzfile(save.path),
            row.names = F, quote = F, sep = "\t")





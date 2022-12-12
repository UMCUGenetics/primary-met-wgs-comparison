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


## Sample metadata--------------------------------
sample_metadata <- read.delim(
  paste0(base_dir,'/passengers/processed/metadata/main/metadata_20221111_1251.txt.gz'),
  na.strings=c('NA','#N/A')
)


# Prepare the metadata
sample_metadata <- sample_metadata[!sample_metadata$is_blacklisted,]
row.names(sample_metadata) <- 1:nrow(sample_metadata)


sample_metadata[,"cohort"] <- factor(sample_metadata[,"cohort"], levels = cohort_order)
sample_metadata[,"tissue_group"] <- factor(sample_metadata[,"tissue_group"], levels = tissue_group_order)
sample_metadata[,"cancer_type"] <- factor(sample_metadata[,"cancer_type"], levels = cancer_type_order_fig1)
sample_metadata[,"cancer_type_code"] <- factor(sample_metadata[,"cancer_type_code"], levels = cancer_type_code_order_fig1)
sample_metadata[,"gender"] <- factor(sample_metadata[,"gender"], levels = c("MALE", "FEMALE"))


sample_metadata_hmf <- sample_metadata[sample_metadata$cohort == "Hartwig",]


## Define treatment category conversion table --------------------------------

treatment_conversion <- data.frame(treatment_mech = c("Aromatase inhibitor", "Selective ER modulator", "Anti-ER", "Pyrimidine (ant)agonist",
                                                   "Platinum", "Folate antagonist", "Experimental", "Alkaloid",
                                                   "Taxane", "Multikinase inhibitor", "GnRH (ant)agonist", "Anti-AR",
                                                   "Glucocorticoid", "Alkylating", "Antifolate", "Anthracycline",
                                                   "Anthracycline, alkylating", "Anti-HER2", "mTOR inhibitor", "Somatostatin analogue",
                                                   "Unknown", "Anti-EGFR", "Folinic acid, pyrimidine (ant)agonist, topoisomerase inhibitor", "Microtubule inhibitor",
                                                   "Anti-CTLA-4", "Anti-VEGF", "BRAF inhibitor", "Anti-PD-1",
                                                   "MEK inhibitor", "Vinca Alkaloid", "Folinic acid", "Radionuclide",
                                                   "HIPEC", "Topoisomerase inhibitor", "CDK4/6 inhibitor or placebo", "Antitumor antibiotic",
                                                   "Anti-PDGFR-?", "Progestogen", "Immune therapy or placebo", "PI3K inhibitor",
                                                   "CDK4/6 inhibitor", "Experimental or placebo", "Immunomodulator", "PARP inhibitor",
                                                   "FGFR inhibitor", "Monoclonal antibody", "Oncolytic virus", "ALK/ROS1 inhibitor",
                                                   "ALK inhibitor", "Anaplastic lymphoma kinase inhibitor, tyrosine kinase inhibitor", "Pyrimidine (ant)agonist, anthracycline, alkylating", "Anti-PD-L1",
                                                   "Pyrimidine (ant)agonist, other", "Bisphosphonate", "FAK inhibitor", "Anti-CD20,MABthera-therapy",
                                                   "Androgen inhibitor", "Other", "Radiotherapy", "Immunetherapy",
                                                   "Pyrimidine (ant)agonist, platinum", "PARP inhibitor or placebo", "", "Anti-RANK-L",
                                                   "Allosteric inhibitor", "Folinic acid, pyrimidine (ant)agonist, platinum", "Pyrimidine (ant)agonist, folinic acid", "Alkylating, vinca alkaloid, alkylating, glucocorticoid",
                                                   "Anthracycline, antitumor antibiotic, vinca alkaloid, alkylating", "Vinca Alkaloid, alkaloid, glucocorticoid, anthracycline", "Alkylating, alkaloid, pyrimidine (ant)agonist, alkylating", "Glucocorticoid, pyrimidine (ant)agonist, platinum, alkylating, anthracycline, alkaloid",
                                                   "PIPAC", "Nucleoside-analogon", "Alkylating, antifolate, pyrimidine (ant)agonist", "Chemoradiotherapy",
                                                   "Taxane, anthracycline, alkylating", "Anti-PD-L1 or placebo", "Pyrimidine (ant)agonist, anthracycline, alkylating, taxane", "Anti-SMO",
                                                   "Platinum, taxane", "Anti-GITR", "Anti-mitotic vinca alkaloid", "Hyperthermia",
                                                   "Anti-PDGFR-? or placebo", "Nuclear therapy or placebo", "Anti-PSMA", "Immunotherapy",
                                                   "Folate antagonist, Platinum", "Vinca Alkaloid, antitumor antibiotic, alkylating", "Anti-PD-1 or placebo", "Alkylating, alkylating, vinca alkaloid",
                                                   "Platinum, Pyrimidine (ant)agonist", "Radiofrequency ablation", "Orteronel or placebo", "Taxane, Bisphosphonate",
                                                   "Anti-CD20, alkylating, anthracycline, vinca alkaloid, glucocorticoid", "Proteasome inhibitor", "Immune modulator", "olinic acid, pyrimidine (ant)agonist, topoisomerase inhibitor, platinum",
                                                   "Aromatase inhibitor or placebo", "PI3K inhibitor or placebo"),
                                treatment_group = c("Hormone_therapy", "Hormone_therapy", "Hormone_therapy", "Chemotherapy",
                                                    "Chemotherapy", "Chemotherapy", "Other", "Chemotherapy",
                                                    "Chemotherapy", "Targeted_therapy", "Hormone_therapy", "Hormone_therapy",
                                                    "Other", "Chemotherapy", "Chemotherapy", "Chemotherapy",
                                                    "Chemotherapy", "Targeted_therapy", "Targeted_therapy", "Hormone_therapy",
                                                    "Unknown", "Targeted_therapy", "Chemotherapy", "Chemotherapy",
                                                    "Immunotherapy", "Targeted_therapy", "Targeted_therapy", "Immunotherapy",
                                                    "Targeted_therapy", "Chemotherapy", "Chemotherapy", "Radiotherapy",
                                                    "Chemotherapy", "Chemotherapy", "Targeted_therapy", "Chemotherapy",
                                                    "Targeted_therapy", "Hormone_therapy", "Immunotherapy", "Targeted_therapy",
                                                    "Targeted_therapy", "Other", "Immunotherapy", "Targeted_therapy",
                                                    "Targeted_therapy", "Targeted_therapy", "Immunotherapy", "Targeted_therapy",
                                                    "Targeted_therapy", "Targeted_therapy", "Chemotherapy", "Immunotherapy",
                                                    "Chemotherapy", "Other", "Targeted_therapy", "Targeted_therapy",
                                                    "Hormone_therapy", "Other", "Radiotherapy", "Immunotherapy",
                                                    "Chemotherapy", "Targeted_therapy", "", "Targeted_therapy",
                                                    "Targeted_therapy", "Chemotherapy", "Chemotherapy", "Chemotherapy",
                                                    "Chemotherapy", "Chemotherapy", "Chemotherapy", "Chemotherapy",
                                                    "Chemotherapy", "Chemotherapy", "Chemotherapy", "Chemotherapy|Radiotherapy",
                                                    "Chemotherapy", "Immunotherapy", "Chemotherapy", "Targeted_therapy",
                                                    "Chemotherapy", "Targeted_therapy", "Chemotherapy", "Other",
                                                    "Targeted_therapy", "Radiotherapy", "Targeted_therapy", "Immunotherapy",
                                                    "Chemotherapy", "Chemotherapy", "Immunotherapy", "Chemotherapy",
                                                    "Chemotherapy", "Radiotherapy", "Hormone_therapy", "Chemotherapy",
                                                    "Chemotherapy|Targeted_therapy", "Targeted_therapy", "Immunotherapy", "Chemotherapy",
                                                    "Hormone_therapy", "Targeted_therapy"))




# Read in the treatment data 
treatment_data <- read.csv(file = paste0(wd, "external-files/treatment_info.tsv"), sep = "\t", header = T, stringsAsFactors = F)

# add the simplified mecahnism to the sample treatment data
treatment_data <- merge(treatment_data, treatment_conversion, by.x = "mechanism", by.y = "treatment_mech", sort = F)
treatment_data <- treatment_data[,c("sample_id", "treatment_group")]
agg_treatment_data <- aggregate(data=treatment_data,treatment_group~.,FUN=paste,collapse="|")


# only keep samples that are included in the main analysis (whitelisted)
agg_treatment_data <- agg_treatment_data[agg_treatment_data$sample_id %in% sample_metadata_hmf$sample_id,]



# radiotherapy info has been incorporated into the metadata beforehand
treatment_metadata <- merge(agg_treatment_data, sample_metadata_hmf[,c("sample_id", "had_radiotherapy")], by = "sample_id", all.y = T)


# add the treatment info to the sample metadata of Hartwig dataset
treatment_metadata$had_other_treatment <- F
treatment_metadata$had_chemotherapy <- F
treatment_metadata$had_hormone_therapy <- F
treatment_metadata$had_targeted_therapy <- F
treatment_metadata$had_immunotherapy <- F


for (i in 1:nrow(treatment_metadata)){
  
  sample <- treatment_metadata$sample_id[i]
  
  if (!is.na(treatment_metadata$treatment_group[treatment_metadata$sample_id == sample])){
    
    if (str_detect(treatment_metadata$treatment_group[treatment_metadata$sample_id == sample], pattern = "Other")) {
      treatment_metadata$had_other_treatment[treatment_metadata$sample_id == sample] <- T
    }
    if (str_detect(treatment_metadata$treatment_group[treatment_metadata$sample_id == sample], pattern = "Chemotherapy")) {
      treatment_metadata$had_chemotherapy[treatment_metadata$sample_id == sample] <- T
    }
    if (str_detect(treatment_metadata$treatment_group[treatment_metadata$sample_id == sample], pattern = "Hormone_therapy")) {
      treatment_metadata$had_hormone_therapy[treatment_metadata$sample_id == sample] <- T
    }
    if (str_detect(treatment_metadata$treatment_group[treatment_metadata$sample_id == sample], pattern = "Targeted_therapy")) {
      treatment_metadata$had_targeted_therapy[treatment_metadata$sample_id == sample] <- T
    }
    if (str_detect(treatment_metadata$treatment_group[treatment_metadata$sample_id == sample], pattern = "Immunotherapy")) {
      treatment_metadata$had_immunotherapy[treatment_metadata$sample_id == sample] <- T
    }
  } else {
    treatment_metadata[treatment_metadata$sample_id == sample,c("had_other_treatment", "had_chemotherapy", "had_hormone_therapy", "had_targeted_therapy", "had_immunotherapy")] <- c(NA, NA, NA, NA, NA)
  }
  
}


# determine for which samples the treatment data is available
treatment_metadata$treatment_info_available <- ifelse(rowSums(is.na(treatment_metadata[,3:8])) == 6, F, T)


# save the data
path_to_directory <- '...'
write.table(treatment_metadata, file = path_to_directory, sep = "\t", row.names = F, quote = F)



## Make the summary plot using ggplot2 --------------------------------

# summarize the treament data
treatment_summary <- data.frame()

for (cancer_type in unique(sample_metadata_hmf$cancer_type)){
  tmp <- sample_metadata_hmf[sample_metadata_hmf$cancer_type == cancer_type,]
  treatment_summary <- rbind(treatment_summary, c(cancer_type, nrow(tmp), sum(tmp$treatment_info_available), sum(tmp$had_chemotherapy, na.rm = T), sum(tmp$had_hormone_therapy, na.rm = T), sum(tmp$had_targeted_therapy, na.rm = T), sum(tmp$had_immunotherapy, na.rm = T), sum(tmp$had_radiotherapy, na.rm = T)))
}


colnames(treatment_summary) <- c("cancer_type", "tot_sample", "tot_treatment", "Chemotherapy", "Hormone_therapy", "Targeted_therapy", "Immunotherapy", "Radiotherapy")
treatment_summary[,2:8] <- apply(treatment_summary[,2:8], 2, as.numeric)
treatment_summary[,3:8] <- 100*treatment_summary[,3:8]/treatment_summary[,2]



# turn the data frame into a tibble to be inputted to ggplot
treatment_summary_tibb <- gather(treatment_summary, key = "treatment_group", value = "percent", 3:8)
treatment_summary_tibb$cancer_type <- factor(treatment_summary_tibb$cancer_type, levels = included_cancer_type_fig1)
treatment_summary_tibb$treatment_group <- factor(treatment_summary_tibb$treatment_group, levels = rev(c("tot_treatment", "Chemotherapy", "Radiotherapy", "Targeted_therapy", "Immunotherapy","Hormone_therapy")))




# make the treatment summary plot
treatment_plot <- ggplot(treatment_summary_tibb, aes(x = treatment_group,y = percent, fill = treatment_group)) + facet_wrap(~cancer_type, nrow = 24) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_manual(values = c("#cc0000", "#999999", "#67a9cf", "#7FC97F", "#ef8a62", "#000000")) +
  ylim(c(0,100)) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  guides(fill = F)
  



# save the plot


path_to_directory <- '...'
pdf(file = path_to_directory, height = 16.75, width = 14)
print(treatment_plot)
dev.off()




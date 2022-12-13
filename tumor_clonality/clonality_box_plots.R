
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

sample_metadata$total_tmb <- rowSums(sample_metadata[,c("sbs_load", "dbs_load", "indel_load")])

sample_metadata_hmf <- sample_metadata[sample_metadata$cohort == "Hartwig",]

# we will plot the clonality based on the biopsy location, therefore samples for which the metstatic location is unknown can be excluded
sample_metadata_hmf <- sample_metadata_hmf[!is.na(sample_metadata_hmf$metastatic_location),]


## Read in the clonality data  --------------------------------
clonality_data <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/all-purple-timing.txt.gz"), stringsAsFactors = F, header = T, sep = "\t")

clonality_data <- merge(clonality_data, sample_metadata[,c("sample_id", "tumor_purity", "cancer_type","cancer_type_code", "cohort", "total_tmb", "metastatic_location")], by = "sample_id")
clonality_data$metastatic_location <- factor(clonality_data$metastatic_location, levels = c("Local", "Lymph", "Distant"))

# low purity samples can complicate our analysis and clonality values
clonality_data <- clonality_data[clonality_data$tumor_purity >= 0.5,]


# remove cancer types that don't reach the minimum sample of 5 for one of the metastatic location groups
for (cancer_type in cancer_type_order[-1]){
  
  if (any(as.numeric(table(clonality_data$metastatic_location[clonality_data$cancer_type == cancer_type])) < 5)){
    clonality_data <- clonality_data[clonality_data$cancer_type != cancer_type,]
  }
  
}


clonality_data$cancer_type <- as.character(clonality_data$cancer_type)


# count the number of samples for each metastatic location of each cancer type
sample_nr <- data.frame()

for (cct in unique(clonality_data$cancer_type)){
  for (bps in c("Local", "Lymph", "Distant")){
    sam_nr <- nrow(clonality_data[clonality_data$cancer_type == cct & clonality_data$metastatic_location == bps,])
    sample_nr <- rbind(sample_nr, c(cct, bps, sam_nr))
  }
}

colnames(sample_nr) <- c("cancer_type", "metastatic_location", "sample_size")



# summarize statistical test results into a data frame for annotating the plot 

signif_df <- data.frame(metastatic_location = c(rep(c("Local", "Lymph", "Distant"), times = 7 )),
                        cancer_type = c(rep(c("Breast carcinoma", "Colorectal carcinoma", "Esophageal carcinoma", "Kidney renal clear cell carcinoma", "Lung adenocarcinoma", "Prostate carcinoma", "Skin melanoma"), each = 3)),
                        significance = c(NA, "p.val=0.01", "p.val<0.001", NA, "p.val=0.03", "p.val=0.06", NA, "ns", "p.val=0.01",
                                         NA , "ns", "p.val=0.03", NA, "ns", "ns", NA, "ns", "ns", NA, "na", "na"))




# make the plot using ggplot2
clonality_boxplot <- clonality_data %>% ggplot(aes(x = metastatic_location, y = 100*(clonal_tmb_80_cutoff/total_tmb), color = metastatic_location)) +
  facet_wrap(~cancer_type, scales = "free") +
  geom_boxplot(outlier.shape = NA) +
  stat_boxplot(geom = "errorbar",
               width = 0.25) +
  geom_jitter(alpha = 0.4) +
  scale_y_continuous(breaks = seq(0, 100, by = 25), labels = c("0", "25", "50", "75", "100"), limits = c(0,110))+
  stat_compare_means(aes(group=simplified_biopsiy_site), ref.group = "Local", label = "p.format", method = 'wilcox', label.y = 100) +
  scale_color_manual(values = c("#8f6014ff", "#14368fff", "#8f1460ff"), name = "Biopsy Site")+
  geom_text(data = sample_nr, y = 10, aes(label = paste0("N=", sample_size)), inherit.aes = T, size = 4, color = "black") +
  theme_bw() +
  xlab("\n Metastatic Tumor Biopsy Site") +
  ylab("Clonality (%)") +
  ggtitle("Clonality Based on Metastatic Site \n") +
  theme(plot.title = element_text(size = 35, hjust = 0.5),
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20),
        axis.text.x = element_text(angle = 90, size = 15, vjust = 0.6, color = "black"),
        axis.text.y = element_text(angle = 90, size = 15, vjust = 0.6, color = "black"),
        legend.title=element_text(size=20),
        legend.text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 15),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        strip.background = element_blank(),
        plot.margin = margin(2,2,2.2,2, "cm")) 




# save the plot
path_to_directory <- '...'
pdf(file = path_to_directory, height = 16.75, width = 14)
print(age_plot)
dev.off()


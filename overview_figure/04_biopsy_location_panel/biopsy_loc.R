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

sample_metadata_hmf <- sample_metadata[sample_metadata$cohort == "Hartwig",]


# read in the biopsy data and merge with sample metadata of Hartwig cohort
vv <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external-files/biopsy_loc_info_raw.tsv"),
               stringsAsFactors = F, header = T, sep = "\t")

sample_metadata_hmf <- merge(sample_metadata_hmf, vv[,c("sampleId", "biopsySite")], by.x = "sample_id", by.y = "sampleId", all.x = T)


# summarize the biopsy data by binning it into 4 groups: local, lymph, distant, and unknown based on the specific sites annotated in the clinical data
sample_metadata_hmf$simplified_biopsiy_site <- NA


for (i in 2:24){
  cc_type <- cancer_type_order_fig1[i]
  
  if ( i == 2){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite != "ovarium"] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "ovarium"] <- "Distant"
  }
  
  if ( i == 3){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Primary"] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "null"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  
  if ( i == 4){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "null"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 5){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite %in% c("Primary","Lung")] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "null"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 6){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite %in% c("Primary","Lung")] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "null"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 7){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "null"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 8){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite %in% c("Primary", "Breast")] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "null"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 9){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Primary"] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite %in% c("null", "other")] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 10){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite %in% c("Primary", "Liver")] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 11){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite %in% c("Primary", "pancreas")] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "null"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 12){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "null"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 13){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Colorectum"] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite %in% c("Unknown", "Other site")] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 14){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite %in% c("Primary", "oesophagus")] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "null"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 15){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite %in% c("Primary", "stomach", "Stomach")] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "null"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 16){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite %in% c("Primary", "Kidney")] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "null"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 17){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite %in% c("Primary", "cervix")] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "null"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 18){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite %in% c("Primary", "Ovarium CA", "ovary")] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "null"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 19){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite %in% c("Primary")] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "null"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 20){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite %in% c("Primary")] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite %in% c("null", "unknown")] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 21){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite %in% c("Primary", "prostate")] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "null"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 22){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "null"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 23){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "null"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
  
  if ( i == 24){
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Skin/Subcutaneous"] <- "Local"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Unknown"] <- "Unknown"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & sample_metadata_hmf$biopsySite == "Lymph node"] <- "Lymph"
    sample_metadata_hmf$simplified_biopsiy_site[sample_metadata_hmf$cancer_type == cc_type & is.na(sample_metadata_hmf$simplified_biopsiy_site)] <- "Distant"
  }
}


# only keep the relevant information
sample_metadata_hmf <- sample_metadata_hmf[,c(1, 3, 11, 51:52)]

# save the data
path_to_directory <- '...'
write.table(sample_metadata_hmf, file = path_to_directory, sep = "\t", row.names = F, quote = F)



# summarize the biopsy data for plotting
biopsy_summary <- table(sample_metadata_hmf$cancer_type, sample_metadata_hmf$simplified_biopsiy_site)
biopsy_summary <- as.data.frame.matrix(biopsy_summary) 

# set the order
biopsy_summary$cancer_type_code <- cancer_type_code_order_fig1[match(row.names(biopsy_summary), cancer_type_order_fig1)]
biopsy_summary <- biopsy_summary[,c(5,1:4)]
row.names(biopsy_summary) <- 1:24





biopsy_summary <- cbind(biopsy_summary, Total = rowSums(biopsy_summary[,2:5]))

# add the pan-cacner information
biopsy_summary_com <- rbind(c("PAN", sum(biopsy_summary$Local), sum(biopsy_summary$Lymph), sum(biopsy_summary$Distant), sum(biopsy_summary$Unknown), sum(biopsy_summary$Total)),biopsy_summary)

# ensure the numbers are numeric
biopsy_summary_com[,2:6] <- apply(biopsy_summary_com[,2:6], 2, as.numeric)

# calculate the percentages
biopsy_summary_com[,2:5] <- 100*biopsy_summary_com[,2:5]/biopsy_summary_com[,6]

# turn the data frame into a tibble 
biopsy_summary_tibb <- gather(biopsy_summary_com, key = "biopsy_site", value = "percent", 2:5)

biopsy_summary_tibb$cancer_type_code <- factor(biopsy_summary_tibb$cancer_type_code, levels = cancer_type_code_order_fig1)
biopsy_summary_tibb$biopsy_site <- factor(biopsy_summary_tibb$biopsy_site, levels = rev(c("Local", "Lymph", "Distant", "Unknown")))

# this dummy variable is used only for representation purposes
biopsy_summary_tibb$dummy <- 1



# make the plot using ggplot2
biopsy_plot <- ggplot(biopsy_summary_tibb, aes(x = dummy,y = percent, fill = biopsy_site)) + facet_wrap(~cancer_type_code, nrow = 24) +
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  scale_fill_manual(values = rev(c("#8f6014", "#14368f", "#8f1460", "#99a3a4"))) +
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
print(biopsy_plot)
dev.off()



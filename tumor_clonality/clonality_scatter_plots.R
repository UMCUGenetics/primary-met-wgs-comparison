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

## Read in the clonality data  --------------------------------

clonality_data <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/all-purple-timing.txt.gz"), stringsAsFactors = F, header = T, sep = "\t")



clonality_data <- merge(clonality_data, sample_metadata[,c("sample_id", "cancer_type","cancer_type_code", "cancer_subtype", "cohort", "total_tmb", "smnv_load", "metastatic_location")], by = "sample_id")


# select the portion of the Hartwig data that we would like to compare against PCAWG

met_sites <- c("All_Hartwig", "Local", "Lymph", "Distant")

met_site_of_interest <- met_sites[1]

# only keep the selected samples
if (met_site_of_interest != "All_Hartwig"){
  
  clonality_data <- clonality_data[clonality_data$metastatic_location %in% c("Primary", met_site_of_interest),]
  
}

# format the data frame
clonality_data <- clonality_data[!is.na(clonality_data$subclonal),]
clonality_data <- clonality_data[clonality_data$subclonal_tmb_80_cutoff != 0,]
cancer_types <- as.character(unique(clonality_data$cancer_type))
cancer_type_abb <- as.character(unique(clonality_data$cancer_type_code))
row_nr <- length(cancer_types)


# summarize the test results of comparisons
clonality_comparison_summary <- data.frame(significance = character(row_nr), cancer_type_abb = character(row_nr), cancer_type_complete = character(row_nr), mean_clona_met = numeric(row_nr), min_clona_met = numeric(row_nr), max_clona_met = numeric(row_nr), mean_clona_prim = numeric(row_nr), min_clona_prim = numeric(row_nr), max_clona_prim = numeric(row_nr), p_val = numeric(row_nr), cohort_size = numeric(row_nr))
for (i in 1:row_nr){
  print(i)
  cancer_type <- cancer_types[i]
  cancer_type_abb <- as.character(unique(sample_metadata$cancer_type_code[sample_metadata$cancer_type == cancer_type]))
  cancer_type_complete <- paste0(cancer_type, " (", cancer_type_abb , ")")
  tmp_df <- clonality_data[clonality_data$cancer_type == cancer_type,]
  
  cohort_size <- nrow(tmp_df)
  clonality_comparison_summary[i,c(1:3,11)] <- c(NA, cancer_type_abb, cancer_type_complete,cohort_size)

  
    if (nrow(tmp_df[tmp_df$cohort == "PCAWG",]) >= 5 & nrow(tmp_df[tmp_df$cohort == "Hartwig",]) >= 5){
    tmp_df_pcawg <- tmp_df$clonal_tmb_80_cutoff[tmp_df$cohort == "PCAWG"]/tmp_df$smnv_load[tmp_df$cohort == "PCAWG"]
    pcawg.model <- lm(tmp_df_pcawg ~ 1)
    conf.int.pcwg <- confint(pcawg.model, level=0.95)
    
    
    tmp_df_hmf <- tmp_df$clonal_tmb_80_cutoff[tmp_df$cohort == "Hartwig"]/tmp_df$smnv_load[tmp_df$cohort == "Hartwig"]
    hmf.model <- lm(tmp_df_hmf ~ 1)
    conf.int.hmf <- confint(hmf.model, level=0.95)
    
    res <- wilcox.test(tmp_df_pcawg, tmp_df_hmf)
    
    clonality_comparison_summary[i,4:10] <- c(mean(tmp_df_hmf, na.rm = T), as.numeric(conf.int.hmf), mean(tmp_df_pcawg, na.rm = T), as.numeric(conf.int.pcwg), res$p.value)
    
  } else {
    clonality_comparison_summary[i,4:9] <- rep(NA, times = 6)
    clonality_comparison_summary[i,10] <- NA
  }
}


# apply the multiple test correction
clonality_comparison_summary[,10] <- p.adjust(clonality_comparison_summary[,10], method = "hochberg")

# post-processing the data frame
clonality_comparison_summary[clonality_comparison_summary$p_val > 0.05 & !is.na(clonality_comparison_summary$p_val),1] <- "black"
clonality_comparison_summary[clonality_comparison_summary$p_val <= 0.05 & !is.na(clonality_comparison_summary$p_val),1] <- "red"

clonality_comparison_summary$cancer_type_abb <- factor(clonality_comparison_summary$cancer_type_abb, levels = cancer_type_code_order)
clonality_comparison_summary$significance <- factor(clonality_comparison_summary$significance, levels = c("red", "black"))
clonality_comparison_summary$cancer_type_complete <- factor(clonality_comparison_summary$cancer_type_complete)
clonality_comparison_summary$cohort_size <- as.numeric(clonality_comparison_summary$cohort_size)
clonality_comparison_summary$enrichment <- clonality_comparison_summary$mean_clona_met > clonality_comparison_summary$mean_clona_prim
clonality_comparison_summary$enrichment <- ifelse(clonality_comparison_summary$enrichment, "Metastatic", "Primary")
clonality_comparison_summary$enrichment <- factor(clonality_comparison_summary$enrichment, levels = c("Primary", "Metastatic"))


# make the scatter plot
plloott_mean2 <- ggplot(clonality_comparison_summary, aes(x = 100*mean_clona_prim, y = 100*mean_clona_met, fill = log2(mean_clona_met/mean_clona_prim))) +
  geom_point(aes(size = cohort_size), color = clonality_comparison_summary$significance, pch=21, stroke = 2) +
  scale_fill_gradient2(low = cohort_color_palette[1],
                       mid = "white",
                       high = cohort_color_palette[2],
                       midpoint = 0,
                       space = "Lab",
                       na.value = "transparent",
                       guide = "colourbar",
                       aesthetics = "fill",
                       breaks=c(-1,-0.5,0,0.5,1),labels=c("-1","-0.5","0","0.5","1"),
                       limits=c(-1,1),
                       name = "Log2 Clonality Ratio (Metastatic/Primary)")+
  scale_size_continuous(range = c(6, 20), name = "Cohort Size", breaks = c(64,128,256,512), labels = c("64", "128", "256", "512")) +
  ggrepel::geom_text_repel(aes(label = cancer_type_abb, color = "black"), max.overlaps = 40, size = 7, box.padding = unit(1.5, "lines")) + 
  scale_color_manual(values = "black") +
  xlim(c(50,100)) +
  xlab("Primary (%)") +
  ylim(c(50,100)) +
  ylab(paste0(met_site_of_interest," Metastatic (%)")) +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw()+
  theme(plot.title = element_text(size = 45, hjust = 0.5),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(angle = 0, size = 25, vjust = 0.6, color = "black"),
        axis.text.y = element_text(angle = 0, size = 25, vjust = 0.6, color = "black"),
        legend.title=element_text(size=20),
        legend.text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        plot.margin = margin(0.5,0.5,3,3, "cm"))  +
  guides(color = FALSE) + theme(legend.position="none")






##### ***********************************************************************************************************

## Comparison of clonality for BRCA - COREAD - UCEC subtype groups   --------------------------------



sample_metadata_subtypes <- sample_metadata
sample_metadata_subtypes$cancer_subtype[is.na(sample_metadata_subtypes$cancer_subtype)] <- "NA"


sample_metadata_subtypes$cancer_subtype[sample_metadata_subtypes$cancer_type_code == "BRCA"] <- paste0("BRCA_", sample_metadata_subtypes$cancer_subtype[sample_metadata_subtypes$cancer_type_code == "BRCA"])
sample_metadata_subtypes$cancer_subtype[sample_metadata_subtypes$cancer_type_code == "COREAD"] <- paste0("COREAD_", sample_metadata_subtypes$cancer_subtype[sample_metadata_subtypes$cancer_type_code == "COREAD"])
sample_metadata_subtypes$cancer_subtype[sample_metadata_subtypes$cancer_type_code == "UCEC"] <- paste0("UCEC_", sample_metadata_subtypes$cancer_subtype[sample_metadata_subtypes$cancer_type_code == "UCEC"])


## Read in the clonality data  --------------------------------

clonality_data <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/all-purple-timing.txt.gz"), stringsAsFactors = F, header = T, sep = "\t")
clonality_data <- merge(clonality_data, sample_metadata_subtypes[,c("sample_id", "cancer_type","cancer_type_code", "cohort", "smnv_load", "cancer_subtype")], by = "sample_id")

# make  copy of clonality data to contain the cancer type groups of interest and combine it with cancer subtype groups
clonality_data_cp <- clonality_data

clonality_data_cp <- clonality_data_cp[clonality_data_cp$cancer_type_code %in% c("BRCA", "COREAD", "UCEC"),] 
clonality_data_cp <- clonality_data_cp[clonality_data_cp$cancer_subtype != "BRCA_NA",]

clonality_data_cp$cancer_type_code <- clonality_data_cp$cancer_subtype

clonality_data <- clonality_data[clonality_data$cancer_type_code %in% c("BRCA", "COREAD", "UCEC"),] 

clonality_data <- rbind(clonality_data, clonality_data_cp)

clonality_data$cancer_type_code <- as.character(clonality_data$cancer_type_code)
clonality_data$cancer_type_code <- factor(clonality_data$cancer_type_code)



# summarize the test results of comparisons
row_nr <- 11
cancer_types <- levels(clonality_data$cancer_type_code)

clonality_comparison_summary <- data.frame(significance = character(row_nr), cancer_type_abb = character(row_nr), cancer_type_complete = character(row_nr), mean_clona_met = numeric(row_nr), min_clona_met = numeric(row_nr), max_clona_met = numeric(row_nr), mean_clona_prim = numeric(row_nr), min_clona_prim = numeric(row_nr), max_clona_prim = numeric(row_nr), p_val = numeric(row_nr), cohort_size = numeric(row_nr))
for (i in 1:row_nr){
  print(i)
  cancer_type_abb <- cancer_types[i]
  cancer_type <- cancer_type_abb
  cancer_type_complete <- paste0(cancer_type, " (", cancer_type_abb , ")")
  tmp_df <- clonality_data[clonality_data$cancer_type_code == cancer_type_abb & clonality_data$subclonal_tmb_80_cutoff != 0,]
  
  cohort_size <- nrow(tmp_df)
  clonality_comparison_summary[i,c(1:3,11)] <- c(NA, cancer_type_abb, cancer_type_complete,cohort_size)
  if (nrow(tmp_df[tmp_df$cohort == "PCAWG",]) >= 5 & nrow(tmp_df[tmp_df$cohort == "Hartwig",]) >= 5){
    tmp_df_pcawg <- tmp_df$clonal_tmb_80_cutoff[tmp_df$cohort == "PCAWG"]/tmp_df$smnv_load[tmp_df$cohort == "PCAWG"]
    pcawg.model <- lm(tmp_df_pcawg ~ 1)
    conf.int.pcwg <- confint(pcawg.model, level=0.95)
    
    
    tmp_df_hmf <- tmp_df$clonal_tmb_80_cutoff[tmp_df$cohort == "Hartwig"]/tmp_df$smnv_load[tmp_df$cohort == "Hartwig"]
    hmf.model <- lm(tmp_df_hmf ~ 1)
    conf.int.hmf <- confint(hmf.model, level=0.95)
    
    res <- wilcox.test(tmp_df_pcawg, tmp_df_hmf)
    clonality_comparison_summary[i,4:10] <- c(mean(tmp_df_hmf, na.rm = T), as.numeric(conf.int.hmf), mean(tmp_df_pcawg, na.rm = T), as.numeric(conf.int.pcwg), res$p.value)
  } else {
    clonality_comparison_summary[i,4:9] <- rep(NA, times = 2)
    clonality_comparison_summary[i,10] <- res$p.value
  }
}

# apply the multiple test correction
clonality_comparison_summary[,10] <- p.adjust(clonality_comparison_summary[,10], method = "hochberg")

# post-processing the data frame
clonality_comparison_summary[clonality_comparison_summary$p_val > 0.05 & !is.na(clonality_comparison_summary$p_val),1] <- "black"
clonality_comparison_summary[clonality_comparison_summary$p_val <= 0.05 & !is.na(clonality_comparison_summary$p_val),1] <- "red"


clonality_comparison_summary$cancer_type_abb <- factor(clonality_comparison_summary$cancer_type_abb, levels = levels(clonality_data$cancer_type_code))
clonality_comparison_summary$significance <- factor(clonality_comparison_summary$significance, levels = c("red", "black"))
clonality_comparison_summary$cancer_type_complete <- factor(clonality_comparison_summary$cancer_type_complete)
clonality_comparison_summary$cohort_size <- as.numeric(clonality_comparison_summary$cohort_size)



clonality_comparison_summary$enrichment <- clonality_comparison_summary$mean_clona_met > clonality_comparison_summary$mean_clona_prim
clonality_comparison_summary$enrichment <- ifelse(clonality_comparison_summary$enrichment, "Metastatic", "Primary")

clonality_comparison_summary$enrichment <- factor(clonality_comparison_summary$enrichment, levels = c("Primary", "Metastatic"))

# make the plot
plloott_mean2 <- clonality_comparison_summary %>% ggplot(aes(x = 100*mean_clona_prim, y = 100*mean_clona_met, fill = log2(mean_clona_met/mean_clona_prim))) +
  geom_point(aes(size = cohort_size), color = clonality_comparison_summary$significance, pch=21, stroke = 2) +
  scale_fill_gradient2(low = cohort_color_palette[1],
                       mid = "white",
                       high = cohort_color_palette[2],
                       midpoint = 0,
                       space = "Lab",
                       na.value = "transparent",
                       guide = "colourbar",
                       aesthetics = "fill",
                       breaks=c(-1,-0.5,0,0.5,1),labels=c("-1","-0.5","0","0.5","1"),
                       limits=c(-1,1),
                       name = "Log2 Clonality Ratio (Metastatic/Primary)")+
  scale_size_continuous(range = c(6, 20), name = "Cohort Size", breaks = c(64,128,256,512), labels = c("64", "128", "256", "512")) +
  ggrepel::geom_text_repel(aes(label = cancer_type_abb, color = "black"), max.overlaps = 40, size = 7, box.padding = unit(1.5, "lines")) +
  scale_color_manual(values = "black") +
  xlim(c(50,100)) +
  xlab("Primary (%)") +
  ylim(c(50,100)) +
  ylab("Metastatic (%)") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw()+
  theme(plot.title = element_text(size = 45, hjust = 0.5),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(angle = 0, size = 25, vjust = 0.6, color = "black"),
        axis.text.y = element_text(angle = 0, size = 25, vjust = 0.6, color = "black"),
        legend.title=element_text(size=20),
        legend.text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        plot.margin = margin(0.5,0.5,3,3, "cm"))  +
  guides(color = FALSE) + theme(legend.position="none")





##### ***********************************************************************************************************

## Comparison of clonality for PANET - PRAD progression states  --------------------------------

sample_metadata_progression <- sample_metadata[sample_metadata$cancer_type_code %in% c("PANET", "PRAD"),]


# keep the two cancer types with progression states and make a data frame containing both the sample cancer type and progression state
sample_metadata_progression$progression_status[is.na(sample_metadata_progression$progression_status) & sample_metadata_progression$cancer_type_code == "PRAD"] <- "PRAD_NA"
sample_metadata_progression$progression_status[sample_metadata_progression$progression_status == "metastatic" & sample_metadata_progression$cancer_type_code == "PRAD"] <- "PRAD_metastatic"
sample_metadata_progression$progression_status[sample_metadata_progression$progression_status == "metastatic" & sample_metadata_progression$cancer_type_code == "PANET"] <- "PANET_metastatic"


sample_metadata_progression$progression_status[sample_metadata_progression$progression_status == "progression/relapse" & sample_metadata_progression$cancer_type_code == "PRAD"] <- "PRAD_progression/relapse"
sample_metadata_progression$progression_status[sample_metadata_progression$progression_status == "stable/remission" & sample_metadata_progression$cancer_type_code == "PRAD"] <- "PRAD_stable/remission"



sample_metadata_progression$progression_status[sample_metadata_progression$progression_status == "progression/relapse" & sample_metadata_progression$cancer_type_code == "PANET"] <- "PANET_progression/relapse"
sample_metadata_progression$progression_status[sample_metadata_progression$progression_status == "stable/remission" & sample_metadata_progression$cancer_type_code == "PANET"] <- "PANET_stable/remission"


sample_metadata_progression <- sample_metadata_progression[sample_metadata_progression$progression_status != "PRAD_NA",]


sample_metadata_progression$progression_status <- factor(sample_metadata_progression$progression_status, levels = c(
  "PRAD_metastatic", "PRAD_stable/remission", "PRAD_progression/relapse",
  "PANET_metastatic", "PANET_stable/remission", "PANET_progression/relapse"))




clonality_data <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/r-objects/all-purple-timing.txt.gz"), stringsAsFactors = F, header = T, sep = "\t")


clonality_data <- merge(clonality_data, sample_metadata_progression[,c("sample_id", "cancer_type","cancer_type_code", "cohort", "progression_status", "smnv_load", "cancer_subtype")], by = "sample_id")



clonality_data$cancer_type_code <- as.character(clonality_data$cancer_type_code)
clonality_data$cancer_type_code <- as.factor(clonality_data$cancer_type_code)





# summarize the test results of comparisons
row_nr <- 4
progress_stats <- c("PANET_progression/relapse", "PANET_stable/remission", "PRAD_progression/relapse", "PRAD_stable/remission")

clonality_comparison_summary <- data.frame(significance = character(row_nr), progress_stat = character(row_nr), cancer_type_complete = character(row_nr), mean_clona_met = numeric(row_nr), min_clona_met = numeric(row_nr), max_clona_met = numeric(row_nr), mean_clona_prim = numeric(row_nr), min_clona_prim = numeric(row_nr), max_clona_prim = numeric(row_nr), p_val = numeric(row_nr), cohort_size = numeric(row_nr))

for (i in 1:4){
  print(i)
  progress_stat <- progress_stats[i]
  cancer_type_code <- unique(clonality_data$cancer_type_code[clonality_data$progression_status == progress_stat])
  cancer_type_complete <- paste0(progress_stat, " (", progress_stat , ")")
  tmp_df <- clonality_data[clonality_data$cancer_type_code ==  cancer_type_code & clonality_data$subclonal_tmb_80_cutoff != 0,]
  
  cohort_size <- nrow(tmp_df)
  clonality_comparison_summary[i,c(1:3,11)] <- c(NA, progress_stat, cancer_type_complete,cohort_size)
  if (cancer_type_code == "PRAD"){
    tmp_df_pcawg <- tmp_df$clonal_tmb_80_cutoff[tmp_df$progression_status == progress_stat]/tmp_df$smnv_load[tmp_df$progression_status == progress_stat]
    pcawg.model <- lm(tmp_df_pcawg ~ 1)
    conf.int.pcwg <- confint(pcawg.model, level=0.95)
    
    
    tmp_df_hmf <- tmp_df$clonal_tmb_80_cutoff[tmp_df$progression_status == "PRAD_metastatic"]/tmp_df$smnv_load[tmp_df$progression_status == "PRAD_metastatic"]
    hmf.model <- lm(tmp_df_hmf ~ 1)
    conf.int.hmf <- confint(hmf.model, level=0.95)
    
    res <- wilcox.test(tmp_df_pcawg, tmp_df_hmf)
    clonality_comparison_summary[i,4:10] <- c(mean(tmp_df_hmf, na.rm = T), as.numeric(conf.int.hmf), mean(tmp_df_pcawg, na.rm = T), as.numeric(conf.int.pcwg), res$p.value)
  } else {
    
    tmp_df_pcawg <- tmp_df$clonal_tmb_80_cutoff[tmp_df$progression_status == progress_stat]/tmp_df$smnv_load[tmp_df$progression_status == progress_stat]
    pcawg.model <- lm(tmp_df_pcawg ~ 1)
    conf.int.pcwg <- confint(pcawg.model, level=0.95)
    
    
    tmp_df_hmf <- tmp_df$clonal_tmb_80_cutoff[tmp_df$progression_status == "PANET_metastatic"]/tmp_df$smnv_load[tmp_df$progression_status == "PANET_metastatic"]
    hmf.model <- lm(tmp_df_hmf ~ 1)
    conf.int.hmf <- confint(hmf.model, level=0.95)
    
    res <- wilcox.test(tmp_df_pcawg, tmp_df_hmf)
    clonality_comparison_summary[i,4:10] <- c(mean(tmp_df_hmf, na.rm = T), as.numeric(conf.int.hmf), mean(tmp_df_pcawg, na.rm = T), as.numeric(conf.int.pcwg), res$p.value)
  }
}

# apply the multiple test correction
clonality_comparison_summary[,10] <- p.adjust(clonality_comparison_summary[,10], method = "hochberg")

# post-processing the data frame
clonality_comparison_summary[clonality_comparison_summary$p_val > 0.05 & !is.na(clonality_comparison_summary$p_val),1] <- "black"
clonality_comparison_summary[clonality_comparison_summary$p_val <= 0.05 & !is.na(clonality_comparison_summary$p_val),1] <- "red"

clonality_comparison_summary$progress_stat <- factor(clonality_comparison_summary$progress_stat, levels = progress_stats)
clonality_comparison_summary$significance <- factor(clonality_comparison_summary$significance, levels = c("red", "black"))
clonality_comparison_summary$cancer_type_complete <- factor(clonality_comparison_summary$cancer_type_complete)
clonality_comparison_summary$cohort_size <- as.numeric(clonality_comparison_summary$cohort_size)

clonality_comparison_summary$enrichment <- clonality_comparison_summary$mean_clona_met > clonality_comparison_summary$mean_clona_prim
clonality_comparison_summary$enrichment <- ifelse(clonality_comparison_summary$enrichment, "Metastatic", "Primary")
clonality_comparison_summary$enrichment <- factor(clonality_comparison_summary$enrichment, levels = c("Primary", "Metastatic"))

# make the plot
plloott_mean2 <- clonality_comparison_summary %>% ggplot(aes(x = 100*mean_clona_prim, y = 100*mean_clona_met, fill = log2(mean_clona_met/mean_clona_prim))) +
  geom_point(aes(size = cohort_size), color = clonality_comparison_summary$significance, pch=21, stroke = 2) +
  scale_fill_gradient2(low = cohort_color_palette[1],
                       mid = "white",
                       high = cohort_color_palette[2],
                       midpoint = 0,
                       space = "Lab",
                       na.value = "transparent",
                       guide = "colourbar",
                       aesthetics = "fill",
                       breaks=c(-1,-0.5,0,0.5,1),labels=c("-1","-0.5","0","0.5","1"),
                       limits=c(-1,1),
                       name = "Log2 Clonality Ratio (Metastatic/Primary)")+
  scale_size_continuous(range = c(6, 20), name = "Cohort Size", breaks = c(64,128,256,512), labels = c("64", "128", "256", "512")) +
  ggrepel::geom_text_repel(aes(label = progress_stat, color = "black"), max.overlaps = 40, size = 7, box.padding = unit(1.5, "lines")) + 
  scale_color_manual(values = "black") +
  xlim(c(50,100)) +
  xlab("Primary (%)") +
  ylim(c(50,100)) +
  ylab("Metastatic (%)") +
  geom_abline(slope = 1,intercept = 0, lty = 2) +
  theme_bw()+
  theme(plot.title = element_text(size = 45, hjust = 0.5),
        axis.title.x = element_text(size = 30),
        axis.title.y = element_text(size = 30),
        axis.text.x = element_text(angle = 0, size = 25, vjust = 0.6, color = "black"),
        axis.text.y = element_text(angle = 0, size = 25, vjust = 0.6, color = "black"),
        legend.title=element_text(size=20),
        legend.text = element_text(size=15),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(color = 'black'),
        plot.margin = margin(0.5,0.5,3,3, "cm"))  +
  guides(color = FALSE) + theme(legend.position="none")







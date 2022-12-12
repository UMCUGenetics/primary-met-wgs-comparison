## Init --------------------------------

expanding_path <- '/home/ali313/Documents/studies/master/umc-project'

path_prefix <- c(
  hpc  ='/hpc/',
  local = paste0(expanding_path, '/hpc/')
)
path_prefix <- path_prefix[min(which(dir.exists(path_prefix)))]
path_prefix <- sub('/hpc/','',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')

## Dependencies + saved objects --------------------------------

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


# summarize sex proportions in a new data frame
row_nr <- length(unique(cohort_order))*length(cancer_type_order_fig1)
sex_df <- data.frame(cancer_type = character(row_nr), cohort = character(row_nr), male = numeric(row_nr), female = numeric(row_nr))

for (i in 1:length(cancer_type_order_fig1)){
  for (j in 1:length(cohort_order)){
    
    cancer_type <- cancer_type_order_fig1[i]
    cohort <- cohort_order[j]
    
    selected_row <- i+(i-1) + (j-1)
    sex_df[selected_row,1:2] <- c(cancer_type, cohort)
    
    if (i != 1){
      sex_df[selected_row,3:4] <- 100*c(nrow(sample_metadata[sample_metadata$cancer_type == cancer_type & sample_metadata$cohort == cohort &  sample_metadata$gender == "MALE",]),
                                           nrow(sample_metadata[sample_metadata$cancer_type == cancer_type & sample_metadata$cohort == cohort &  sample_metadata$gender == "FEMALE",]))/nrow(sample_metadata[sample_metadata$cancer_type == cancer_type & sample_metadata$cohort == cohort,])
    } else {
      sex_df[selected_row,3:4] <- 100*c(nrow(sample_metadata[sample_metadata$cohort == cohort &  sample_metadata$gender == "MALE",]),
                                           nrow(sample_metadata[sample_metadata$cohort == cohort &  sample_metadata$gender == "FEMALE",]))/nrow(sample_metadata[sample_metadata$cohort == cohort,])
    }
  }
}


# prepare the data frame as an input for ggplot2
sex_df <- tidyr::gather(sex_df, key = "gender", value = "proportion", 3:4)
sex_df$cancer_type <- factor(sex_df$cancer_type, levels = cancer_type_order_fig1)
sex_df$cohort <- factor(sex_df$cohort, levels = cohort_order)
sex_df$gender <- factor(sex_df$gender, levels = c("male", "female"))


# make the plot using ggplot2
sex_plot <- ggplot(sex_df, aes(x = cohort,y = proportion, fill = gender)) + facet_wrap(~cancer_type, nrow = 48) +
  geom_bar(stat="identity") +
  scale_color_manual(values = cohort_color_palette, labels = cohort_order) +
  scale_fill_manual(values = c("#8498f8", "#fc7571")) +
  guides(fill=FALSE) +
  coord_flip() +
  scale_y_continuous(breaks = seq(0, 100, by = 50), labels = c("0", "50", "100%")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'white'))



# save the plot


path_to_directory <- '...'
pdf(file = path_to_directory, height = 16.75, width = 14)
print(sex_plot)
dev.off()


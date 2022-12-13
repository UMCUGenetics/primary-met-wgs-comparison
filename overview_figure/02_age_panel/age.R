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


# prepare the metadata
sample_metadata <- sample_metadata[!sample_metadata$is_blacklisted,]
row.names(sample_metadata) <- 1:nrow(sample_metadata)

sample_metadata[,"cohort"] <- factor(sample_metadata[,"cohort"], levels = cohort_order)
sample_metadata[,"tissue_group"] <- factor(sample_metadata[,"tissue_group"], levels = tissue_group_order)
sample_metadata[,"cancer_type"] <- factor(sample_metadata[,"cancer_type"], levels = cancer_type_order_fig1)
sample_metadata[,"cancer_type_code"] <- factor(sample_metadata[,"cancer_type_code"], levels = cancer_type_code_order_fig1)
sample_metadata[,"gender"] <- factor(sample_metadata[,"gender"], levels = c("MALE", "FEMALE"))


# read in and process raw Hartwig age information
hartwig_meta <- read.csv(file = paste0(base_dir, "drivers/analysis/dna-rep-ann/final-update/external-files/hmf-clinical-metadata-20211112.tsv"), stringsAsFactors = F, header = T, sep = "\t")

hartwig_meta$biopsyDateProcessed <- NA
hartwig_meta$biopsyDateProcessed <- unlist(lapply(str_split(hartwig_meta$biopsyDate, pattern = "-"), "[[", 1))
hartwig_meta$age_at_biopsy <- NA

hartwig_meta[,c("birthYear","biopsyDateProcessed")] <- apply(hartwig_meta[,c("birthYear","biopsyDateProcessed")], as.numeric, MARGIN = 2)
hartwig_meta$age <- hartwig_meta$biopsyDateProcessed - hartwig_meta$birthYear

colnames(hartwig_meta)[2] <- "sample_id"
hartwig_meta <- hartwig_meta[,c(2,ncol(hartwig_meta))]


# read in and process raw PCWG age information
pcawg_meta <- read.csv(file = paste0("/home/ali313/Documents/studies/master/umc-project/hpc/cuppen/shared_resources/PCAWG/Metadata/donor.info.followup.all_projects.csv"), header = T, stringsAsFactors = F, sep = "\t")


colnames(pcawg_meta)[c(1,9)] <- c("sample_id", "age")
pcawg_meta <- pcawg_meta[,c(1,9)]


# merge the age info of both cohorts to the sample metadata object
age_meta_all <- rbind(hartwig_meta, pcawg_meta)
sample_metadata <- merge(sample_metadata, age_meta_all, by = "sample_id")


# incorporate the "pan_cancer" group into the sample metadata object
sample_metadata_pan <- sample_metadata
sample_metadata_pan$cancer_type <- "Pan-cancer"
sample_metadata <- rbind(sample_metadata, sample_metadata_pan)


# make the v line position dataset
p_meds <- ddply(sample_metadata, .(cancer_type, cohort), summarise, med = median(age, na.rm = T))


# make the plot using ggplot2
age_plot <- ggplot(sample_metadata, aes(x = age, color = cohort))+ facet_wrap(~ cancer_type, ncol = 1)+
  geom_density(aes(color = cohort), size = 3) +
  scale_color_manual(values = c("#E66711","#7217C5")) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_blank()) +
  theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
  geom_vline(data=p_meds, aes(xintercept=med, colour=cohort),
             linetype=1, size=5) +
  guides(color = F)


# save the plot


path_to_directory <- '...'
pdf(file = path_to_directory, height = 16.75, width = 14)
print(age_plot)
dev.off()


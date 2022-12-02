############### Extract the subclonal likelihood ratio from hartwig processed VCF files ############### 
# author: remy (sascha)
# date: 08/06/2021
# last updated: 16/02/2022

### Description
# This script extracts the subclonal likelihood ratio given by purple output
# then categorizes variants into clonal (likelihood < 0.8) and subclonal (likelihood >= 0.8)
# to calculate the amount of clonal and subclonal variants in each sample.
# Each VCF is processed separately.

### Input
# The manifest file of postprocessed VCF file paths you generated in step 9

### Output
# Table containing the sample_id, the mutation type, clonality and the number of variants in each category.

### Usage
# sbatch calc_tmb_from_purple_vcf.sh -m <manifest_file>

# libs
source(paste0(here::here(), '/code/r_objects/libs.R'))

#========= Path prefixes =========#
base_dir <- list(
  path=paste0(here::here())
)

for(i in base_dir){
  if(dir.exists(i)){
    base_dir <- i
    break
  }
}

# pass commandline directory arguments to 'args' vector
args = commandArgs(trailingOnly=TRUE)

# Directories
sample_id <- args[1]
purple_hartwig_input <- args[2]
output_dir <- args[3]

# on hpc
vcf <- mutSigExtractor::variantsFromVcf(vcf.file = purple_hartwig_input, 
                                        merge.consecutive = T,
                                        vcf.filter = 'PASS', 
                                        vcf.fields = c('CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'INFO'))



# add the subclonal likelihood ratio as a separate column
subcl <- mutSigExtractor::getInfoValues(vcf$info, keys = 'SUBCL')
subcl <- as.numeric(subcl)

# combine the subcl and vcf
vcf <- cbind(vcf, subcl)

# fill NA values with zeros (these are considered clonal; they are not reported)
vcf <- vcf %>% tidyr::replace_na(., replace = list(subcl = 0))

# rename that column, then add another column indicating whether this variant is clonal or subclonal or NA
vcf <- vcf %>%
  dplyr::rename(subclonal_likelihood = 'subcl') %>%
  mutate(clonality = if_else(subclonal_likelihood >= 0.8, 'subclonal', 'clonal')) %>%
  select(-info)

# add another column 'mut_type' that shows the mutation type label
vcf$mut_type <- detSmnvType(vcf$ref, vcf$alt)

# count the TMB of SNVs, MNVs and Indels, each split up clonality
vcf_tmb_per_mut_type_per_clonality <- vcf %>%
  group_by(mut_type, clonality) %>%
  summarise(n_per_mut_type_and_clonality = n()) %>%
  ungroup()

# add the sample name as a separate column
vcf_tmb_per_mut_type_per_clonality <- vcf_tmb_per_mut_type_per_clonality %>%
  mutate(sample_id = sample_id, .before = mut_type)

# write to disk
write.table(vcf_tmb_per_mut_type_per_clonality, 
            file = paste0(output_dir, sample_id, '_clonality_annotated.temp'),
            sep = '\t', quote = FALSE, append = FALSE, row.names = FALSE,
            col.names = FALSE)

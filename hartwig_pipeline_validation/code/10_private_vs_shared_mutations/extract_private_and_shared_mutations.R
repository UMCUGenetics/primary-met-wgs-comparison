############### Extract the subclonal likelihood ratio from Hartwig processed VCF files ############### 
# author: remy (sascha)
# date: 01/11/2021
# last updated: 02/11/2021

### Description
# This script extracts the private and shared variant calls between the PCAWG and Hartwig (HMF) pipeline on a per sample level.
# Now also extracting clonality of HMF VCFs based on the SUBCL field.

### Input
# The manifest file you generated in step 1.

### Output
# Table of variant counts that are exclusively called by PCAWG or Hartwig pipeline, or shared per sample.

### Usage
# sbatch extract_private_and_shared_mutations.sh

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
sample <- args[1]
hartwig_vcf_file_input <- args[2]
pcawg_vcf_file_input <- args[3]
pcawg_vcf_indel_file_input <- args[4]
output_dir <- args[5]

# read in the VCF
pcawg_vcf <- mutSigExtractor::variantsFromVcf(vcf.file = pcawg_vcf_file_input,
                                              merge.consecutive = T,
                                              vcf.fields = c('CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'INFO'))

pcawg_indel_vcf <- mutSigExtractor::variantsFromVcf(vcf.file = pcawg_vcf_indel_file_input,
                                                    merge.consecutive = T,
                                                    vcf.fields = c('CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'INFO'))

hartwig_vcf <- mutSigExtractor::variantsFromVcf(vcf.file = hartwig_vcf_file_input,
                                            merge.consecutive = T,
                                            vcf.filter = 'PASS',
                                            vcf.fields = c('CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'INFO'))

# combine first four columns into one to create a unique ID for each variant
hartwig_vcf <- hartwig_vcf %>%
  unite(col = 'var_id', c(chrom, pos, ref, alt))
pcawg_vcf <- pcawg_vcf %>%
  unite(col = 'var_id', c(chrom, pos, ref, alt))
pcawg_indel_vcf <- pcawg_indel_vcf %>%
  unite(col = 'var_id', c(chrom, pos, ref, alt))

# combine indel and SNV PCAWG VCFs
pcawg_vcf <- rbind(pcawg_vcf, pcawg_indel_vcf)

# drop filter and info column, add pipeline columns for both vcfs
pcawg_vcf <- pcawg_vcf %>% select(-filter, -info) %>% mutate(pipeline1 = 'PCAWG')
hartwig_subcl <- getInfoValues(hartwig_vcf$info, keys = 'SUBCL')
hartwig_vcf$hartwig_subclonal_likelihood <- as.numeric(hartwig_subcl)
hartwig_vcf <- hartwig_vcf %>% select(-info) %>% mutate(pipeline2 = 'Hartwig') %>%
  mutate(hartwig_subclonal_likelihood = if_else(is.na(hartwig_subclonal_likelihood), 0, hartwig_subclonal_likelihood)) %>%
  mutate(hartwig_clonality = if_else(hartwig_subclonal_likelihood >= 0.8, 'subclonal', 'clonal')) %>%
  select(-hartwig_subclonal_likelihood)

# full join to get all the shared and private variants
vcf <- pcawg_vcf %>%
  full_join(hartwig_vcf, by = 'var_id')

# mark private and shared variants based on missing values in pipeline columns
vcf <- vcf %>%
  mutate(hartwig_clonality = if_else(is.na(hartwig_clonality), 'none', hartwig_clonality)) %>%
  mutate(privacy = case_when(is.na(pipeline1) & pipeline2 == 'Hartwig' ~ 'Hartwig_only',
                             pipeline1 == 'PCAWG' & is.na(pipeline2) ~ 'PCAWG_only',
                             TRUE ~ 'shared')) %>%
  select(-pipeline1, -pipeline2) %>%
  separate(col = 'var_id', into = c('chrom', 'pos', 'ref', 'alt'), sep = '_')

# get the mutation type based on ref and alt column
vcf$mut_type <- detSmnvType(vcf$ref, vcf$alt)

# summarise the table, count the number of variants per mut_type that are private to either pipeline and shared
vcf_summarised <- vcf %>%
  group_by(privacy, mut_type, hartwig_clonality) %>%
  summarise(n = n(), .groups = 'drop') %>%
  complete(privacy, mut_type, hartwig_clonality, fill = list(n = 0)) %>%
  mutate(sample_id = sample, .before = 'privacy') %>%
  # post filtering
  filter(privacy != 'Hartwig_only' | hartwig_clonality != 'none',
         privacy != 'PCAWG_only' | hartwig_clonality != 'clonal',
         privacy != 'PCAWG_only' | hartwig_clonality != 'subclonal')

# write to disk
write.table(vcf_summarised, 
            file = paste0(output_dir, sample, '_variant_privacy.temp'),
            sep = '\t', quote = FALSE, append = FALSE, row.names = FALSE,
            col.names = FALSE)
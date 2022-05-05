############### Extract the subclonal likelihood ratio & gene name from HMF processed VCF files ############### 
# author: remy (sascha)
# date: 07/10/2021
# last updated: 11/10/2021

### Description
# input: needs a 
# This script extracts the subclonal likelihood ratio and gene name (SEC field) given by purple VCF output
# then categorizes variants into clonal (likelihood < 0.8) and subclonal (likelihood >= 0.8)
# Each VCF is processed and stored separately.
# Only genes in the HMF gene panel list are considered.

# set options
options(stringsAsFactors = F)

# libs
library(dplyr) # data manipulation
library(tidyr) # data cleaning
library(readr) # reading / writing data
library(mutSigExtractor) # reading and manipulating VCF files

#========= Path prefixes =========#
base_dir <- list(
  hpc='/base/path',
  mnt='/base/path',
  umc='/base/path'
)

for(i in base_dir){
  if(dir.exists(i)){
    base_dir <- i
    break
  }
}

# pass commandline directory arguments to 'args' vector
args = commandArgs(trailingOnly=TRUE)

# input variables
sample <- args[1]
purple_hmf_input <- args[2]
output_dir <- args[3]
linx_catalog_date <- args[4]

# read in the newest linx driver catalog
linx_drivers <- read_tsv(file = paste0(base_dir, '/path/to/', linx_catalog_date, '/linx_drivers.tsv'))

# on hpc
vcf <- mutSigExtractor::variantsFromVcf(vcf.file = purple_hmf_input, 
                                        merge.consecutive = T, 
                                        vcf.filter = 'PASS', 
                                        vcf.fields = c('CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'INFO'))



# add the subclonal likelihood ratio as a separate column
subcl <- mutSigExtractor::getInfoValues(vcf$info, keys = c('SEC', 'SUBCL', 'TIER'))
subcl$SUBCL <- as.double(subcl$SUBCL)

# combine the subcl and vcf
vcf <- cbind(vcf, subcl)

# fill NA values with zeros (these are considered clonal; they are not reported)
vcf <- vcf %>% tidyr::replace_na(., replace = list(SUBCL = 0))

# rename columns, then add another column indicating whether this variant is clonal or subclonal or NA
vcf <- vcf %>%
  dplyr::rename(subclonal_likelihood = 'SUBCL',
                gene = 'SEC',
                tier = 'TIER')

# drop the NA rows in the gene column, clean gene column, add clonality column based on subcl cutoff threshold
vcf <- drop_na(vcf, gene) %>%
  mutate(gene = str_replace(gene, pattern = '\\,.+', replacement = '')) %>%
  mutate(clonality = if_else(subclonal_likelihood >= 0.8, 'subclonal', 'clonal'))

# filter only for those genes in the driver gene panel, choose relevant columns
vcf <- vcf %>%
  # filter(gene %in% gene_panel) %>%
  select(gene, subclonal_likelihood, clonality, tier) %>%
  mutate(driver = 'MUTATION')

# add the sample name as a separate column
vcf <- vcf %>%
  mutate(sample_id = sample, .before = gene)

# only retain variants that were identified as drivers by linx
vcf_final <- vcf %>%
  semi_join(linx_drivers, by = c('sample_id', 'gene', 'driver'))

# on hpc
write.table(vcf_final, 
            file = paste0(output_dir, sample, '_clonality_annotated.temp'),
            sep = '\t', quote = FALSE, append = FALSE, row.names = FALSE,
            col.names = FALSE)
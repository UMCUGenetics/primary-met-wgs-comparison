############### Postprocess PCAWG VCF files ############### 
# author: remy (sascha)
# date: 11/02/2022
# last updated: 14/02/2022

### Description
# This script loads the individual VCF files of PCAWG samples, removes duplicated variants at same position, 
# and prioritizes the most base substitutions at each position 
# (e.g. if on chrom 1 pos 5000 there is a SBS and a DBS reported then only the DBS is retained).

### Input
# The manifest file you generated in step 9.
# The corresponding VCF file per sample.

### Output
# A post processed VCF file

### Usage
# sbatch calc_tmb_from_purple_vcf.sh

# libs
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(mutSigExtractor)

#========= Path prefixes =========#
base_dir <- list(
  hpc='/hpc/cuppen/projects/P0022_HMF_validation',
  mnt='/home/sascha/hpc_mnt/projects/P0022_HMF_validation',
  umc='/home/cog/sbrunner/hpc/cuppen/projects/P0022_HMF_validation'
)

for(i in base_dir){
  if(dir.exists(i)){
    base_dir <- i
    break
  }
}

# print(base_dir)

# pass commandline directory arguments to 'args' vector
args = commandArgs(trailingOnly=TRUE)

# Directories
sample <- args[1]
vcf_body <- args[2]
purple_hmf_input <- args[3]
output_dir <- args[4]

print(c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
        str_replace(sample, pattern = 'T', replacement = 'R'), sample))

# read in the VCF
# # local
# vcf_file <- '/home/sascha/hpc_mnt/shared_resources/PCAWG/pipeline5/per-donor/DO218351-from-jar/purplesoft3.3/DO218351T.purple.somatic.vcf.gz'
# vcf_file <- '/home/cog/sbrunner/hpc/cuppen/shared_resources/PCAWG/pipeline5/per-donor/DO8047-from-jar/purplesoft3.3/DO8047T.purple.somatic.vcf.gz'
# 
# vcf <- mutSigExtractor::readVcfFields(vcf.file = vcf_file,
#                                       fields = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
#                                                  str_replace('DO218351T', pattern = 'T', replacement = 'R'), 'DO218351T'))
# vcf <- mutSigExtractor::readVcfFields(vcf.file = vcf_file,
#                                       fields = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
#                                                  str_replace('DO8047T', pattern = 'T', replacement = 'R'), 'DO8047T'))

# on hpc
vcf <- mutSigExtractor::readVcfFields(vcf.file = purple_hmf_input,
                                      fields = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
                                                 str_c(sample, 'R'), str_c(sample, 'T')))
# vcf <- mutSigExtractor::variantsFromVcf(vcf.file = purple_hmf_input, 
#                                         merge.consecutive = T,
#                                         vcf.filter = 'PASS', 
#                                         vcf.fields = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
#                                                        str_c(sample, 'R'), str_c(sample, 'T')))

chromosome_order <- c(1:22, 'X', 'Y')

# remove duplicates, prioritize highest number of base substitutions per chromosome and position
vcf_filter <- vcf %>%
  # filter(FILTER == 'PASS') %>%
  distinct() %>%
  group_by(CHROM, POS) %>%
  # prioritize longest indel / base substitution event
  slice_max(., order_by = REF, n = 1) %>%
  # prioritize longest indel event
  slice_max(., order_by = ALT, n = 1) %>%
  # prioritize PON over PASS variants
  slice_max(., order_by = FILTER, n = 1) %>%
  arrange(factor(CHROM, levels = chromosome_order))

# vcf_count <- vcf %>%
#   # filter(FILTER == 'PASS') %>%
#   # distinct() %>%
#   group_by(CHROM, POS) %>%
#   count()

# output_dir <- base_dir

# write this table to disk
# local
# write.table(vcf_filter, 
#             file = paste0(output_dir, '/snv_mnv/processed/calc_tmb_from_purple_vcf/results/DO38635_clonality_annotated.tsv'),
#             sep = '\t', quote = FALSE, append = FALSE, row.names = FALSE,
#             col.names = FALSE)

print(paste0(output_dir, '/', vcf_body))

# on hpc
# write.table(vcf_filter, 
#             file = paste0(output_dir, '/', vcf_body),
#             sep = '\t', quote = FALSE, append = TRUE, row.names = FALSE,
#             col.names = FALSE)

write_tsv(vcf_filter,
          file = paste0(output_dir, '/', vcf_body), 
          col_names = FALSE, append = FALSE)

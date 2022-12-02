############### Postprocess downsampled Hartwig VCF files ############### 
# author: remy (sascha)
# date: 23/03/2022
# last updated: 24/03/2022

### Description
# This script loads the individual VCF files of downsampled Hartwig samples (at 30X, 38X and 60X), 
# removes duplicated variants at same position, 
# prioritizes the highest base substitution at each position, (e.g. if on chrom 1 pos 5000 there is a SBS and a DBS reported then only the DBS is retained),
# and also removes low tumor quality variants in close proximity to a variant of higher priority (e.g. if at position 5000 there is a DBS reported, at position 5001 a MNV and at position 5002 a DBS then only the MNV is retained)

### Input
# The manifest file you generated in step 15.1.

### Output
# A postprocessed VCF file per sample with the above stated filters applied to it.

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
sample <- args[1]
vcf_body <- args[2]
purple_hartwig_input <- args[3]
output_dir <- args[4]

# read in the VCF
vcf <- mutSigExtractor::readVcfFields(vcf.file = purple_hartwig_input,
                                      fields = c('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT',
                                                 str_replace(sample, pattern = regex('[A-Z]$'), replacement = 'R'), str_c(sample, 'T')))
# def chrom order
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
  ungroup() %>%
  arrange(factor(CHROM, levels = chromosome_order))

drop_variant_table <- vcf_filter %>%
  mutate(lag_pos_diff = POS - lag(POS), .after = POS) %>%
  mutate(lead_pos_diff = lead(POS) - POS, .after = POS) %>%
  mutate(ref_length = nchar(REF), .after = lag_pos_diff) %>%
  filter(lag_pos_diff %in% c(1,2,3) | lead_pos_diff %in% c(1,2,3))

# logical column that indicates whether a variant should be dropped or not
drop_variant_table <- drop_variant_table %>%
  mutate(drop_variant = case_when(
    lead_pos_diff == 1 & ref_length > lead(ref_length) & lag(ref_length) < lag_pos_diff ~ FALSE,
    lead_pos_diff == 1 & ref_length > lead(ref_length) ~ FALSE, # catch the first row too (missed by the line above)
    lead_pos_diff == 2 & ref_length > 2 & ref_length >= lead(ref_length)  ~ FALSE,
    lead_pos_diff == 3 & ref_length > 3 & ref_length > lead(ref_length) ~ FALSE,
    lead_pos_diff == 1 & lead_pos_diff >= ref_length ~ FALSE,
    lead_pos_diff == 2 & lead_pos_diff >= ref_length & lag(ref_length) < lag_pos_diff ~ FALSE,
    lead_pos_diff == 3 & lead_pos_diff >= ref_length & lag(ref_length) < lag_pos_diff ~ FALSE,
    lag_pos_diff %in% c(1,2,3) & lag_pos_diff >= lag(ref_length) ~ FALSE,
    ref_length > lag_pos_diff & ref_length > lag(ref_length) ~ FALSE,
    TRUE ~ TRUE
  ), .after = ref_length)

drop_variant_table <- drop_variant_table %>%
  filter(drop_variant)

# anti_join the drop_variant_table to the original vcf_filter table to remove unwanted variants
vcf_filter <- vcf_filter %>%
  anti_join(drop_variant_table, by = c('CHROM', 'POS', 'REF', 'ALT'))

# write to disk
write_tsv(vcf_filter,
          file = paste0(output_dir, '/', vcf_body), 
          col_names = FALSE, append = FALSE)

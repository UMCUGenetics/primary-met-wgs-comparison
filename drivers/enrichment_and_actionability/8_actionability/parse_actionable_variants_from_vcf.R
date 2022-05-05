############### Extract potentially actionable variants from VCF files ############### 
# author: remy (sascha)
# date: 16/12/2021
# last updated: 26/01/2021

### Description
# input: VCF file
# This script parses the SNPEff annotated VCF INFO field. It extracts the Gene that is affected by the 
# mutation and further extract the base change and amino acid residue change.

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

# --------------------------------------------- AMINO ACID CODES (needed for string cleaning of the protein residue changes column)

# read in amino acid code table
amino_acids <- read_tsv(file = paste0(base_dir, '/path/to/amino_acids.tsv'))

# create names vector of three letter + single letter aa code pairs
# (the substitute string must be in the value slot, not the name slot)
amino_acids <- setNames(amino_acids$short, amino_acids$abbr)

# on hpc
vcf <- mutSigExtractor::variantsFromVcf(vcf.file = purple_hmf_input, 
                                        merge.consecutive = T, 
                                        vcf.filter = 'PASS', 
                                        vcf.fields = c('CHROM', 'POS', 'REF', 'ALT', 'FILTER', 'INFO'))



# extract SNPEff gene annotation from the info column
info <- mutSigExtractor::getInfoValues(vcf$info, keys = c('ANN', 'TIER'))

# test del split ANN
info <- info %>%
  separate(., col = 'ANN', into = str_c('transcript_', 1:1000), sep = ',')

# combine the info and vcf tables
vcf <- cbind(vcf, info)
rm(info)

# rename columns
vcf <- vcf %>%
  dplyr::rename(
    tier = 'TIER')

# split the gene column by comma into six columns
vcf <- drop_na(vcf, transcript_1) %>%
  pivot_longer(., cols = str_c('transcript_', 1:1000), names_to = 'names', 
               values_to = 'snpeff_ann', values_drop_na = TRUE) %>%
  separate(snpeff_ann, into = c('drop_1', 'var_type', 'drop_2', 'gene', 
                                str_c('drop_', 3:6), 'intron_exon_rank', 'drop_7', 'protein_residue_change'), sep = '\\|') %>%
  filter(protein_residue_change != '', var_type != 'synonymous_variant')

# clean protein_residue_change column
vcf <- vcf %>%
  mutate(protein_residue_change = str_remove(protein_residue_change, pattern = regex('^p.'))) %>%
  mutate(protein_residue_change = str_replace_all(protein_residue_change,
                                                  amino_acids))

# add the sample name as a separate column
vcf <- vcf %>%
  mutate(sample_id = sample, .before = gene)

# only keep relevant columns
vcf_final <- vcf %>%
  select(sample_id, gene, var_type, intron_exon_rank, protein_residue_change) %>%
  distinct()

# on hpc
write.table(vcf_final, 
            file = paste0(output_dir, sample, '_actionable_variants.temp'),
            sep = '\t', quote = FALSE, append = FALSE, row.names = FALSE,
            col.names = FALSE)
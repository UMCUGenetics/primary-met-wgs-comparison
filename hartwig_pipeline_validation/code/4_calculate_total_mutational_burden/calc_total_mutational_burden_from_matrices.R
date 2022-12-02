############### Calculate total mutational burden for extracted mutational signatures ###############
# author: remy (sascha)
# date: 24/05/2021
# last updated: 01/10/2021

### Description
# This script calculates the total mutational burden for both the PCAWG and hartwig processed 
# PCAWG data and adds them to a table called 'total_mutational_burden.tsv'
# Total mutational burden can be calculated for any mutational signature (SNV, Indel, DBS...)

### Input
# The manifest file that contains the paths to the generate mutational context matrices.

### Output
# One table for each mutation type (SBS, DBS, INDEL, SV) that contains the sample_id and the mutational burden derived from the respective
# PCAWG and Hartwig context matrices.

### Usage
# sbatch calc_total_mutational_burden_from_matrices.sh

# libraries
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

# parse commandline input
sample <- args[1]
snv_hartwig_input <- args[2]
dbs_hartwig_input <- args[3]
indel_hartwig_input <- args[4]
sv_hartwig_input <- args[5]
sv_hartwig_len_cutoff_input <- args[6]
snv_pcawg_input <- args[7]
dbs_pcawg_input <- args[8]
indel_pcawg_input <- args[9]
sv_pcawg_input <- args[10]
sv_pcawg_len_cutoff_input <- args[11]
output_dir <- args[12]

# read.delim function with fixed arguments
# function to write table with same arguments
read_delim_fixed <- function(...) { read.delim(..., sep = '\t', header = T) }

# read in the hartwig matrices
snv_hartwig_matrix <- read_delim_fixed(file = snv_hartwig_input)
dbs_hartwig_matrix <- read_delim_fixed(file = dbs_hartwig_input)
indel_hartwig_matrix <- read_delim_fixed(file = indel_hartwig_input)
sv_hartwig_matrix <- read_delim_fixed(file = sv_hartwig_input)
sv_hartwig_len_cutoff_matrix <- read_delim_fixed(file = sv_hartwig_len_cutoff_input)
# read in the pcawg matrices
snv_pcawg_matrix <- read_delim_fixed(file = snv_pcawg_input)
dbs_pcawg_matrix <- read_delim_fixed(file = dbs_pcawg_input)
indel_pcawg_matrix <- read_delim_fixed(file = indel_pcawg_input)
sv_pcawg_matrix <- read_delim_fixed(file = sv_pcawg_input)
sv_pcawg_len_cutoff_matrix <- read_delim_fixed(file = sv_pcawg_len_cutoff_input)

# function to calculate and extract the TMB from context matrices for any mutation type,
# pipeline results are stored in separate columns side by side.
# Table is written to disk
calc_tmb <- function(hartwig_matrix, pcawg_matrix, sample = 'DO1000', output_path) {
  
  hartwig_tot_mut_burden <- hartwig_matrix %>%
    summarise(across(everything(), sum, .names = 'hartwig_total_mut_burden'))
  
  pcawg_tot_mut_burden <- pcawg_matrix %>%
    summarise(across(everything(), sum, .names = 'pcawg_total_mut_burden'))
  
  # bind
  hartwig_pcawg_total_mut_burden <- cbind(hartwig_tot_mut_burden, pcawg_tot_mut_burden)
  
  # add patient_id column
  hartwig_pcawg_total_mut_burden_df <- hartwig_pcawg_total_mut_burden %>%
    mutate(patient_id = sample, .before = hartwig_total_mut_burden)
  
  # write total mutational burden table to disk
  write.table(hartwig_pcawg_total_mut_burden_df,
              file = paste0(output_dir, output_path),
              sep = '\t', quote = FALSE, append = TRUE, row.names = FALSE,
              col.names = FALSE)
  
}

### SNV
calc_tmb(snv_hartwig_matrix, snv_pcawg_matrix, sample = sample, 
         output_path = paste0('/snv_total_mutational_burden.tsv'))

### DBS
calc_tmb(dbs_hartwig_matrix, dbs_pcawg_matrix, sample = sample, 
         output_path = paste0('/dbs_total_mutational_burden.tsv'))

### Indel
calc_tmb(indel_hartwig_matrix, indel_pcawg_matrix, sample = sample, 
         output_path = paste0('/indel_total_mutational_burden.tsv'))

### SV
calc_tmb(sv_hartwig_matrix, sv_pcawg_matrix, sample = sample, 
         output_path = paste0('/sv_total_mutational_burden.tsv'))

### SV length cutoff
# filter out small SVs
drop_short_svs <- str_detect(rownames(sv_hartwig_len_cutoff_matrix), pattern = '\\_1e02_bp$|\\_3e02_bp$|\\_5e02_bp$', negate = TRUE)

sv_hartwig_len_cutoff_matrix <- sv_hartwig_len_cutoff_matrix %>%
  filter(drop_short_svs)
sv_pcawg_len_cutoff_matrix <- sv_pcawg_len_cutoff_matrix %>%
  filter(drop_short_svs)

calc_tmb(sv_hartwig_len_cutoff_matrix, sv_pcawg_len_cutoff_matrix, sample = sample, 
         output_path = paste0('/sv_total_mutational_burden_len_cutoff.tsv'))
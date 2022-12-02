############### Calculate cosine similarity from extracted signature matrices ############
# author: remy (sascha)
# date: 24/05/2021
# last updated: 20/09/2021

### Description
# This script loads the extracted mutational signature matrices per donor for both
# PCAWG and Hartwig (hartwig) pipeline binds them together then calculates the cosine similarity
# between the two.

### Input
# The manifest file that contains the paths to the generate mutational context matrices.

### Output
# One table for each mutation type (SBS, DBS, INDEL, SV) that contains the sample_id and the cosine similarity of the respective
# PCAWG and Hartwig context matrices.

### Usage
# sbatch calc_cosine_similarity_from_matrices.sh -m <manifest_file>

# load libraries
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

# parsing commandline input
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

# function to calculate cosine similarity, extract the value and storeit in a dataframe, then write dataframe to disk
calc_cosine_similarity <- function(hartwig_matrix, pcawg_matrix, sample = 'DO1000', column_name = 'cosine_similarity', output_path) {
  
  cosine_similarity <- lsa::cosine(hartwig_matrix[,1], pcawg_matrix[,1])
  cosine_similarity <- cosine_similarity[1,1]
  cosine_similarity_df <- data.frame(patient_id = sample, column_name = cosine_similarity)
  
  # write cosine similarity table to disk
  write.table(cosine_similarity_df, 
              file = paste0(output_dir, output_path),
              sep = '\t', quote = FALSE, append = TRUE, row.names = FALSE,
              col.names = FALSE)
  
}

# calculate and store the cosine similarity values for each mutation type 
# calculated between pipelines per sample

### SNV
calc_cosine_similarity(snv_hartwig_matrix, snv_pcawg_matrix, sample = sample, output_path = paste0('/snv_cosine_similarity.tsv'))

### DBS
calc_cosine_similarity(dbs_hartwig_matrix, dbs_pcawg_matrix, sample = sample, output_path = paste0('/dbs_cosine_similarity.tsv'))

### Indel
calc_cosine_similarity(indel_hartwig_matrix, indel_pcawg_matrix, sample = sample, output_path = paste0('/indel_cosine_similarity.tsv'))

### SV
calc_cosine_similarity(sv_hartwig_matrix, sv_pcawg_matrix, sample = sample, output_path = paste0('/sv_cosine_similarity.tsv'))

### SV length cutoff
# filter out small SVs 
drop_short_svs <- str_detect(rownames(sv_hartwig_len_cutoff_matrix), pattern = '\\_1e02_bp$|\\_3e02_bp$|\\_5e02_bp$', negate = TRUE)

sv_hartwig_len_cutoff_matrix <- sv_hartwig_len_cutoff_matrix %>%
  filter(drop_short_svs)
sv_pcawg_len_cutoff_matrix <- sv_pcawg_len_cutoff_matrix %>%
  filter(drop_short_svs)

calc_cosine_similarity(sv_hartwig_len_cutoff_matrix, sv_pcawg_len_cutoff_matrix, sample = sample, output_path = paste0('/sv_cosine_similarity_len_cutoff.tsv'))
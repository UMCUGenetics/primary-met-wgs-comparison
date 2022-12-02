############### Extract the SV context profile with different length bins ############
# author: remy (sascha)
# date: 08/07/2021
# last updated: 12/07/2021

### Description
# This script extracts the SV context (DEL, DUP, INV, TRA) based on the size of the SV and places them
# in length bins for each mutation type.

### Input
# a manifest file created in step 1

### Output
# sv context matrices

### Usage
# either run all context matrix extractions at once with: sbatch extract_sigs.sh -m <manifest_file> -l <sv_length_cutoff>
# or just the SVs context matrix extraction with: sbatch extract_sv_sigs_length_cutoff.sh -m <manifest_file> -l <sv_length_cutoff>


# global options
options(stringsAsFactors=F)

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

# commandline argument parsing
args <- commandArgs(trailingOnly=TRUE)
sample_name <- args[1]
vcf_sv_hmf_len_cutoff <- args[2]
vcf_sv_pcawg_len_cutoff <- args[3]
length_cutoff <- args[4]
out_dir <- args[5]

# parse to correct type
length_cutoff <- as.numeric(length_cutoff)

#========= Main =========#
# function that locks the same arguments for writing tables
write.tsv <- function(...){ write.table(..., sep='\t', quote=F) }

# file appendices
mut_types <- c('sv_hmf_len_cutoff', 'sv_pcawg_len_cutoff')

# function to extract the SV context profile for each sample per pipeline
main <- function(sample_name, vcf_sv_hmf_len_cutoff, vcf_sv_pcawg_len_cutoff){
  
  # create output paths for both matrices and name the vector
  out_paths <- paste0(out_dir,'/',mut_types,'/',sample_name,'_',mut_types,'.txt')
  names(out_paths) <- mut_types
  
  # initiate list to store matrices
  contexts <- list()
  
  # extract the SV context profile matrices from VCF files and write to disk
  message('## Extracting HMF SV contexts...')
  if(!file.exists(out_paths['sv_hmf_len_cutoff'])){
    dir.create(out_paths['sv_hmf_len_cutoff'])
    contexts$sv_hmf_len_cutoff <- extractSigsSv(vcf_sv_hmf_len_cutoff, output='contexts', vcf.filter='PASS', 
                                                sample.name=sample_name, sv.len.cutoffs = c(0, length_cutoff, 10^c(3:7), Inf))
    write.tsv(contexts$sv_hmf, out_paths['sv_hmf_len_cutoff'])
  }
  
  message('## Extracting PCAWG SV contexts...')
  if(!file.exists(out_paths['sv_pcawg_len_cutoff'])){
    dir.create(out_paths['sv_pcawg_len_cutoff'])
    contexts$sv_pcawg_len_cutoff <- extractSigsSv(vcf_sv_pcawg_len_cutoff, output='contexts', vcf.filter='PASS', 
                                                  sample.name=sample_name, sv.len.cutoffs = c(0, length_cutoff, 10^c(3:7), Inf), sv.caller = "pcawg") ## change the sv caller for PCAWG SV VCFs to pcawg
    write.tsv(contexts$sv_pcawg, out_paths['sv_pcawg_len_cutoff'])
  }
  
}

# execute function
main(sample_name, vcf_sv_hmf_len_cutoff, vcf_sv_pcawg_len_cutoff)

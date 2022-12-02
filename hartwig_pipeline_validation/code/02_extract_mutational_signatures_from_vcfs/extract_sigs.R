############### Extract the mutation context profile for SNVs, MNVs, indels and SVs ############
# author: remy (sascha)
# date: 07/05/2021
# last updated: 20/09/2021

### Description
# This script extracts the mutation context matrices based raw PCAWG/HMF processed VCF files.
# The outputs are one matrix for every mutation type and pipeline per sample (8 matrices total per sample).

### Input
# vcf_paths.txt.gz (manifest file you generated in step 1)

### Output
# Extracted mutational context matrices. One matrix per Hartwig (HMF) and PCAWG cohort and mutation type (SNV; DBS, Indel, SV).
# Output matrices are placed in separate folders.

### Usage
# sbatch extract_sigs.sh -l <sv_langth_cutoff> (e.g. sbatch extract_sigs.sh -l 500)
# for <sv_langth_cutoff> 500 is recommended, because that is the minimum SV length that the PCAWG pipeline calls

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
vcf_som_hmf <- args[2]
vcf_som_pcawg_snv <- args[3]
vcf_som_pcawg_indel <- args[4]
vcf_sv_hmf <- args[5]
vcf_sv_pcawg <- args[6]
output_dir <- args[7]

#========= Main =========#
# function to write table with same arguments
write.tsv <- function(...){ write.table(..., sep='\t', quote=F) }

# define subfolder names where matrices will be stored
mut_types <- c('snv_hmf', 'snv_pcawg', 
               'dbs_hmf', 'dbs_pcawg',
               'indel_hmf','indel_pcawg',
               'sv_hmf', 'sv_pcawg')

# function to extract SNV, DBS, indel and SV context profile matrices in one go from one sample
# for both pipelines
main <- function(sample_name, vcf_som_hmf, vcf_som_pcawg_snv, vcf_som_pcawg_indel, vcf_sv_hmf, vcf_sv_pcawg){
  
  # construct output paths to matrices and name this vector accordingly
  out_paths <- paste0(output_dir,'/',mut_types,'/',sample_name,'_',mut_types,'.txt')
  names(out_paths) <- mut_types
  
  # initate list where matrix paths will be stored
  contexts <- list()
  
  # extract SNV context for both pipelines
  message('## Extracting HMF SNV contexts...')
  if(!file.exists(out_paths['snv_hmf'])){
    dir.create(dirname(out_paths['snv_hmf']))
    contexts$snv_hmf <- extractSigsSnv(vcf_som_hmf, output='contexts', vcf.filter='PASS', sample.name=sample_name)
    write.tsv(contexts$snv_hmf, out_paths['snv_hmf'])
  }
  
  message('## Extracting PCAWG SNV contexts...')
  if(!file.exists(out_paths['snv_pcawg'])){
    dir.create(dirname(out_paths['snv_pcawg']))
    contexts$snv_pcawg <- extractSigsSnv(vcf_som_pcawg_snv, output='contexts', sample.name=sample_name, merge.consecutive = T) # these are already filtered for "PASS"
    write.tsv(contexts$snv_pcawg, out_paths['snv_pcawg'])
  }
  
  # extract DBS context for both pipelines
  message('## Extracting HMF DBS contexts...')
  if(!file.exists(out_paths['dbs_hmf'])){
    dir.create(dirname(out_paths['dbs_hmf']))
    contexts$dbs_hmf <- extractSigsDbs(vcf_som_hmf, output='contexts', vcf.filter='PASS', sample.name=sample_name)
    write.tsv(contexts$dbs_hmf, out_paths['dbs_hmf'])
  }
  
  message('## Extracting PCAWG DBS contexts...')
  if(!file.exists(out_paths['dbs_pcawg'])){
    dir.create(dirname(out_paths['dbs_pcawg']))
    contexts$dbs_pcawg <- extractSigsDbs(vcf_som_pcawg_snv, output='contexts', sample.name=sample_name, merge.consecutive = T) ## vcfs already filtered for 'PSS' & PCAWG vcf files only show SBS for each position, so consecutive SBS need to be merged
    write.tsv(contexts$dbs_pcawg, out_paths['dbs_pcawg'])
  }
  
  # extract Indel context for both pipelines
  message('## Extracting HMF indel contexts...')
  if(!file.exists(out_paths['indel_hmf'])){
    dir.create(dirname(out_paths['indel_hmf']))
    contexts$indel_hmf <- extractSigsIndel(vcf_som_hmf, vcf.filter='PASS', sample.name=sample_name, method = 'PCAWG') ## 'PCAWG' is the newer method
    write.tsv(contexts$indel_hmf, out_paths['indel_hmf'])
  }
  
  message('## Extracting PCAWG indel contexts...')
  if(!file.exists(out_paths['indel_pcawg'])){
    dir.create(dirname(out_paths['indel_pcawg']))
    contexts$indel_pcawg <- extractSigsIndel(vcf_som_pcawg_indel, sample.name=sample_name, method = 'PCAWG') ## 'PCAWG' is the newer method
    write.tsv(contexts$indel_pcawg, out_paths['indel_pcawg'])
  }
  
  # extract SV context for both pipelines
  message('## Extracting HMF SV contexts...')
  if(!file.exists(out_paths['sv_hmf'])){
    dir.create(dirname(out_paths['sv_hmf']))
    contexts$sv_hmf <- extractSigsSv(vcf_sv_hmf, output='contexts', vcf.filter='PASS', sample.name=sample_name)
    write.tsv(contexts$sv_hmf, out_paths['sv_hmf'])
  }
  
  message('## Extracting PCAWG SV contexts...')
  if(!file.exists(out_paths['sv_pcawg'])){
    dir.create(dirname(out_paths['sv_pcawg']))
    contexts$sv_pcawg <- extractSigsSv(vcf_sv_pcawg, output='contexts', vcf.filter='PASS', sample.name=sample_name, sv.caller = "pcawg") ## change the sv caller for PCAWG SV VCFs to manta
    write.tsv(contexts$sv_pcawg, out_paths['sv_pcawg'])
  }
  
}

# execute function with commandline input
main(sample_name, vcf_som_hmf, vcf_som_pcawg_snv, vcf_som_pcawg_indel, vcf_sv_hmf, vcf_sv_pcawg)
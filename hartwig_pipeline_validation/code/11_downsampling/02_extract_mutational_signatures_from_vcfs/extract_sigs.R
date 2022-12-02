############### Extract the mutation context profile for SNVs, MNVs and indels ############
# author: luan

### Description
# This script extracts the mutation context matrices based downsampled Hartwig VCF files.
# The outputs are one context matrix for every mutation type (SNV, DBS, INDEL).

### Input
# downsampled_vcf_paths.txt.gz (manifest file)

### Output
# Extracted mutational context matrices. One matrix per downsampled VCF and mutation type (SNV, DBS, Indel).
# Output matrices are placed in separate folders.

### Usage
# sbatch extract_sigs.sh

# set global options
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

# command line input
args <- commandArgs(trailingOnly=TRUE)
sample_id <- args[1]
output_dir <- args[2]
som_vcf <- args[3]

## Main ================================
write.tsv <- function(...){ write.table(..., sep='\t', quote=F) }

# contruct out paths
mut_types <- c('snv','dbs','indel_chord','indel')
out_paths <- paste0(output_dir,'/',mut_types,'/',sample_id,'_',mut_types,'.txt')
names(out_paths) <- mut_types

# create subdirs and extract signature context per mutation type
for(i in paste0(output_dir,'/',mut_types,'/')){ dir.create(i, showWarnings=F) }

main <- function(sample_id, som_vcf){

   contexts <- list()

   message('## Extracting SNV contexts...')
   if(!file.exists(out_paths['snv'])){
      contexts$snv <- extractSigsSnv(som_vcf, output='contexts', vcf.filter='PASS', sample.name=sample_id)
      write.tsv(contexts$snv, out_paths['snv'])
   }

   message('## Extracting indel contexts (CHORD)...')
   if(!file.exists(out_paths['indel_chord'])){
      contexts$indel_chord <- extractSigsIndel(som_vcf, vcf.filter='PASS', sample.name=sample_id, method='CHORD')
      write.tsv(contexts$indel_chord, out_paths['indel_chord'])
   }
   
   message('## Extracting indel contexts (PCAWG)...')
   if(!file.exists(out_paths['indel'])){
      contexts$indel <- extractSigsIndel(som_vcf, vcf.filter='PASS', sample.name=sample_id, method='PCAWG')
      write.tsv(contexts$indel, out_paths['indel'])
   }

   message('## Extracting DBS contexts...')
   if(!file.exists(out_paths['dbs'])){
      contexts$dbs <- extractSigsDbs(som_vcf, output='contexts', vcf.filter='PASS', sample.name=sample_id)
      write.tsv(contexts$dbs, out_paths['dbs'])
   }
}

main(sample_id, som_vcf)
############### Merge SMNV signature matrices ############
# author: luan

### Description
# This script extracts the SV context matrices based downsampled Hartwig VCF files.

### Input
# individual SMNV context matrices

### Output
# combined SMNV matrix for all downsampled Hartwig samples

### Usage
# just run it

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

input_date <- '2022_11_11'

## Main ================================
output_dir <- paste0(base_dir, '/results/11_downsampling/2_extract_mutational_signatures_from_vcfs/', input_date, '/results/matrices')
sub_dirs <- list.dirs(output_dir, recursive=F)

l_paths <- lapply(sub_dirs, list.files, full.names=T)
names(l_paths) <- basename(sub_dirs)

readContextsFromPaths <- function(paths){
   counter <- 0
   pb <- txtProgressBar(max=length(paths), style=3)
   l <- lapply(paths, function(i){
      counter <<- counter+1
      setTxtProgressBar(pb, counter)
      read.delim(i)
   })
   message('\n')
   return(l)
}

ll_contexts <- list()
ll_contexts$snv <- readContextsFromPaths(l_paths$snv)
ll_contexts$dbs <- readContextsFromPaths(l_paths$dbs)
ll_contexts$indel_chord <- readContextsFromPaths(l_paths$indel_chord)
ll_contexts$indel <- readContextsFromPaths(l_paths$indel)

## Merge into one matrix per mut type
l_contexts <- lapply(ll_contexts, function(i){
   t(do.call(cbind,i))
})

## Export
saveRDS(l_contexts, paste0(output_dir,'/smnv_contexts.rds'))
write.table(
   do.call(cbind,l_contexts), gzfile(paste0(output_dir,'/smnv_contexts.txt.gz')),
   sep='\t',quote=F
)


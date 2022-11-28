options(stringsAsFactors=F)

## Path prefixes ================================
path_prefix <- c(
   hpc  ='/hpc/',
   local=path.expand('~/hpc/')
)
path_prefix <- path_prefix[ min(which(dir.exists(path_prefix))) ]
path_prefix <- sub('/hpc/','/',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')

devtools::load_all(paste0(path_prefix,'/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor'))

## Main ================================
out_dir <- paste0(base_dir,'/passengers/processed/mut_contexts/04_fixed_smnvs/matrices/')
sub_dirs <- list.dirs(out_dir, recursive=F)

l_paths <- lapply(sub_dirs, list.files, full.names=T)
names(l_paths) <- basename(sub_dirs)

# l_sample_names <- lapply(l_paths,basename)
# l_sample_names <- lapply(l_sample_names, function(i){
#    sapply(strsplit(i,'_'),`[[`,1)
# })

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
ll_contexts$sv <- readContextsFromPaths(l_paths$sv)

## Merge into one matrix per mut type
l_contexts <- lapply(ll_contexts, function(i){
   t(do.call(cbind,i))
})

## Export
saveRDS(l_contexts, paste0(out_dir,'/contexts_merged.rds'))
write.table(
   do.call(cbind,l_contexts), gzfile(paste0(out_dir,'/contexts_merged.txt.gz')),
   sep='\t',quote=F
)


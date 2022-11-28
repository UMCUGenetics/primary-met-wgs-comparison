options(stringsAsFactors=F)

## Paths ================================
path_prefix <- c(
   hpc  ='/hpc/',
   local=path.expand('~/hpc/')
)
path_prefix <- path_prefix[ min(which(dir.exists(path_prefix))) ]
path_prefix <- sub('/hpc/','/',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')
wd <- paste0(base_dir,'/passengers/processed/sigs_denovo/extractions/12_fixed_seeds/')

## Load data ================================
out_dir <- paste0(wd,'/sig_contrib/fit_lsq.by_tissue_type/')
dir.create(out_dir, showWarnings=F)

contribs_raw_path <- paste0(out_dir,'/denovo_contribs.lsq.by_sample.rds')
if(!file.exists(contribs_raw_path)){
   fit_out_dir <- paste0(wd,'/sig_contrib/fit_lsq/')
   txt_files <- list.files(fit_out_dir, full.names=T)
   names(txt_files) <- sub('[.]txt$','',basename(txt_files))
   
   counter <- 0
   contribs_raw <- lapply(txt_files, function(i){
      counter <<- counter + 1
      if(counter %% 50 == 0){
         message('Reading [',counter,']: ',basename(i))
      }
      read.delim(i)
   })
   
   saveRDS(contribs_raw, contribs_raw_path)
} else {
   contribs_raw <- readRDS(contribs_raw_path)
}

## Merge signature contrib vectors ================================
fit_metadata <- read.delim(paste0(wd,'/metadata/fit_metadata.txt.gz'))

sample_groups <- split(fit_metadata$sample, fit_metadata$tissue_type_group)
contribs_by_tt <- lapply(sample_groups, function(i){
   df <- do.call(rbind, contribs_raw[i])
   rownames(df) <- i
   return(df)
})
#lapply(contribs_by_tt, colnames)

## RDS
saveRDS(contribs_by_tt, paste0(out_dir,'/denovo_contribs.lsq.by_tissue_type.rds'))

## Txt per tissue type
for(i in names(contribs_by_tt)){
   #i=names(contribs_by_tt)[1]
   write.table(
      contribs_by_tt[[i]],
      gzfile(paste0(out_dir,'/',i,'.txt.gz')),
      sep='\t',quote=F
   )
}

# ## Txt for all samples
# contribs_all <- do.call(function(...) {
#    
#    tmp <- plyr::rbind.fill(...)
#    rownames(tmp) <- lapply(contribs_by_tt, rownames)
#    return(tmp)
#    
# }, contribs_by_tt)
# 
# 
# contribs_all <- do.call(
#    plyr::rbind.fill,
#    contribs_by_tt
# )
# rownames(contribs_all) <- unlist(lapply(contribs_by_tt, rownames), use.names=F)
# 
# write.table(
#    contribs_all,
#    gzfile(paste0(out_dir,'/all_samples.txt.gz')),
#    sep='\t',quote=F
# )

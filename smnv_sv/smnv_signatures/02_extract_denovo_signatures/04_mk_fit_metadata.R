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
manifest <- read.delim(paste0(base_dir,'/passengers/processed/manifest/05_20220216/manifest_HMF_PCAWG.txt.gz'))
metadata <- read.delim(paste0(wd,'/metadata/nmf_metadata.txt'))

## Sig profiles
sig_profile_metadata <- (function(){
   dir <- paste0(wd,'/sig_profiles/')
   paths <- list.files(dir, pattern='[.]txt$', full.names=F)
   
   paths_split <- strsplit(
      sub('[.]txt$','',basename(paths)),
      '[.]'
   )
   df <- as.data.frame(do.call(rbind, paths_split))
   colnames(df) <- c('tissue_type_group','mut_type','rank')
   df$rank <- as.integer(df$rank)
   df$rank <- NULL
   
   df$path <- sub(path_prefix,'',paths)
   
   df <- reshape2::dcast(df, tissue_type_group~mut_type, value.var='path')
   colnames(df) <- c('tissue_type_group','sig_path_DBS','sig_path_ID','sig_path_SBS')
   
   return(df)
})()

## Make manifest for signature extraction ================================
## Init
fit_metadata <- manifest[,c('cohort','sample')]
fit_metadata$tissue_type_group <- metadata$tissue_type_group[match(fit_metadata$sample, metadata$sample_id)]

## Vcf
fit_metadata$sample_dir <- manifest$dir
fit_metadata$som_vcf_path <- manifest$som_vcf
fit_metadata <- fit_metadata[!is.na(fit_metadata$som_vcf_path),]

## Sig profile paths
fit_metadata$sig_profile_dir <- sub(path_prefix,'',paste0(wd,'/sig_profiles/'))
fit_metadata <- cbind(
   fit_metadata,
   sig_profile_metadata[
      match(fit_metadata$tissue_type_group, sig_profile_metadata$tissue_type_group),
      -1
   ]
)

## Export
write.table(
   fit_metadata, gzfile(paste0(wd,'/metadata/fit_metadata.txt.gz')),
   sep='\t',quote=F,row.names=F
)

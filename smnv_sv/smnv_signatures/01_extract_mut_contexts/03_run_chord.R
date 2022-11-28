options(stringsAsFactors=F)

## Path prefixes ================================
path_prefix <- c(
   hpc  ='/hpc/',
   local=path.expand('~/hpc/')
)
path_prefix <- path_prefix[ min(which(dir.exists(path_prefix))) ]
path_prefix <- sub('/hpc/','/',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/')

devtools::load_all(paste0(path_prefix,'/hpc/cuppen/projects/P0013_WGS_patterns_Diagn//CHORD/processed/scripts_main/mutSigExtractor/'))
devtools::load_all(paste0(path_prefix,'/hpc/cuppen/projects/P0013_WGS_patterns_Diagn//CHORD/processed/scripts_main/CHORD/'))

## Main ================================
wd <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/mut_contexts/04_fixed_smnvs/')

contexts_raw <- readRDS(paste0(wd,'/matrices/contexts_merged.rds'))

contexts <- contexts_raw[c('snv','indel_chord','sv')]
names(contexts) <- c('snv','indel','sv')

chord_preds <- chordPredict(contexts)
write.table(chord_preds, paste0(wd,'chord_preds.txt'), sep='\t', row.names=F, quote=F)

## Init ================================
options(stringsAsFactors=F)

WRITE_OUTPUT <- FALSE

## Paths --------------------------------
path_prefix <- c(
   hpc  ='/hpc/',
   local=path.expand('~/hpc/')
)
path_prefix <- path_prefix[ min(which(dir.exists(path_prefix))) ]
path_prefix <- sub('/hpc/','/',path_prefix)

base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')

## Dependencies --------------------------------
source(paste0(base_dir,'/passengers/analysis/common_functions/common_functions.R'))

## Main ================================
purple_purity <- read.delim(paste0(base_dir,'/passengers/processed/metadata/purity/05_20220216/purple.purity.txt.gz'))
loh_prop_in_diploid_samples <- read.delim(paste0(base_dir,'/passengers/processed/chrom_arm_comparison/loh_proportion_diploid.tsv'))
aneuploidy_score <- read.delim(paste0(base_dir,'/passengers/processed/chrom_arm_comparison/aneuploidy_score.tsv'))
tp53_mut <- read.delim(paste0(base_dir,'/passengers/processed/chrom_arm_comparison/tp53_mut_tf.tsv'))

selectSampleRows <- function(df, sample.id.col='sample_id', sel.samples=sample_whitelist, show.sample.ids=F, return.cols=NULL, drop=T){
   #df=loh_prop_in_diploid_samples
   
   out <- df[
      match(sel.samples, df[[sample.id.col]]),
   ,]
   
   missing_samples <- sel.samples[ !(sel.samples %in% df[[sample.id.col]]) ]
   if(length(missing_samples)>0){
      warning(length(missing_samples),' from `sel.samples` were missing from `df`')
      out$sample_id <- sample_whitelist
   }
   
   if(!is.null(return.cols)){
      out <- out[,return.cols, drop=drop]
   }
   
   if(show.sample.ids){
      out <- data.frame(sample_id=sample_whitelist, out)
   }
   
   return(out)
}

#selectSampleRows(aneuploidy_score, sample.id.col='sample_id', return.cols='aneuploidy_score', show.sample.ids=T) |> head()
#selectSampleRows(loh_prop_in_diploid_samples, sample.id.col='sample_id', return.cols='lohProportion', show.sample.ids=T) |> head()
# selectSampleRows(purple_purity, sample.id.col='sample', return.cols='diploidProportion', show.sample.ids=T) |> nrow()
# selectSampleRows(purple_purity, sample.id.col='sample', return.cols='wholeGenomeDuplication', show.sample.ids=T) |> nrow()
# selectSampleRows(tp53_mut, sample.id.col='sample_id', return.cols='tp53_mut', show.sample.ids=T) |> nrow()

df <- data.frame(
   sample_id=sample_whitelist,
   aneuploidy_score=selectSampleRows(aneuploidy_score, sample.id.col='sample_id', return.cols='aneuploidy_score'),
   loh_prop_in_diploid_samples=selectSampleRows(loh_prop_in_diploid_samples, sample.id.col='sample_id', return.cols='lohProportion'),
   has_wgd=as.logical( selectSampleRows(purple_purity, sample.id.col='sample', return.cols='wholeGenomeDuplication') ),
   has_tp53_mut=selectSampleRows(tp53_mut, sample.id.col='sample_id', return.cols='tp53_mut')
)

write.table(
   df, paste0(base_dir,'/passengers/analysis/sv_comparison/08_20221111/ext_data/genome_status.txt'),
   sep='\t', quote=F, row.names=F
)



# df_sacha <- read.delim('/Users/lnguyen/hpc/cuppen/projects/P0025_PCAWG_HMF/drivers/analysis/dna-rep-ann/final-update/r-objects/purple-purity-fraction-changed.txt.gz')
# df_luan <- purple_purity
# 
# tmp1 <- selectSampleRows(df_sacha, sample.id.col='sample_id', show.sample.ids=T, return.cols='diploidProportion')
# tmp2 <- selectSampleRows(df_luan, sample.id.col='sample', show.sample.ids=T, return.cols='diploidProportion')
# tmp <- data.frame(sample_id=sample_whitelist, diploidProportion_sacha=tmp1$out, diploidProportion_luan=tmp2$out)
# subset(tmp, diploidProportion_sacha!=diploidProportion_luan)

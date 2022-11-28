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

nmf_dir <- paste0(base_dir,'/passengers/processed/sigs_denovo/extractions/12_fixed_seeds/')
wd <- paste0(nmf_dir,'/scripts/08_etiology/')
#plots_dir <- paste0(wd,'/plots/')
#dir.create(plots_dir, showWarnings=F)

## Dependencies --------------------------------
#source(paste0(base_dir,'/passengers/analysis/common_functions/common_functions.R'))
#devtools::load_all(paste0(path_prefix,'/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor'))
library(mutSigExtractor)
library(ggplot2)

## Load data ================================
## Metadata 
sample_metadata <- read.delim(
  paste0(base_dir,'/passengers/processed/metadata/main/metadata_20220929_1720.txt.gz'),
  na.strings=c('NA','#N/A')
)

sample_whitelist <- subset(sample_metadata, !is_blacklisted, sample_id, drop=T)

## Signatures
sig_contribs_raw <- read.delim(paste0(nmf_dir,'/sig_contrib/fit_lsq.post_processed/denovo_contribs.lsq.post_processed.txt.gz'), check.names=F)
sig_contribs <- sig_contribs_raw[sample_whitelist,]

sig_metadata <- read.delim(paste0(nmf_dir,'/sig_contrib/fit_lsq.post_processed/sig_metadata.post_processed.txt.gz'))


#data.frame(colnames(sig_contribs))

## Signature cross correlation ================================
plotSigCrossCorr <- function(m, legend.title='correlation'){

  ## Prep data
  m_cor <- cor(m)
  pd <- reshape2::melt(m_cor)
  colnames(pd) <- c('x','y','cor')

  pd$cor[pd$cor<=0] <- 0
  pd$cor[pd$x==pd$y] <- NA

  ## Div lines
  sig_type <- stringr::str_extract(colnames(m),'[A-Z]+(_denovo_clust)*')
  rle_out <- rle(sig_type)
  current_pos <- 0
  div_pos <- sapply(rle_out$lengths, function(i){
    out <- current_pos+i
    current_pos <<- out
    return(out)
  })
  #div_pos <- ncol(m) - div_pos + 0.5

  ## Plot
  ggplot(pd, aes(x=x, y=y, fill=cor)) +
    geom_tile(color='grey') +
    scale_fill_distiller(palette='Spectral', name=legend.title) +

    geom_hline(yintercept=ncol(m)-div_pos+0.5, size=0.3) +
    geom_vline(xintercept=div_pos+0.5, size=0.3) +

    scale_x_discrete(expand=c(0,0), position='top') +
    scale_y_discrete(expand=c(0,0), limits=rev) +
    theme_bw() +
    theme(
      panel.background=element_rect(fill=NA),
      panel.grid=element_blank(),
      axis.text.x.top=element_text(angle=90, hjust=0, vjust=0.5),
      axis.title=element_blank()
    )
}

if(WRITE_OUTPUT){
  m <- sig_contribs
  m <- m/rowSums(m)
  m[is.na(m)] <- 0

  p <- plotSigCrossCorr(m, legend.title='Relative\ncontribution\ncorrelation')

  pdf(paste0(wd,'/sig_cross_corr.pdf'), 12, 10)
  plot(p)
  dev.off()
}

## Manual checks ================================
# df <- sig_contribs['DBS_denovo_1']
# df <- data.frame(sample=rownames(df), contrib=df[,1], row.names=NULL)
# #df <- subset(df, sample %in% sample_whitelist)
# 
# df$msi_status <- sample_metadata$msi_status[ match(df$sample, sample_metadata$sample_id) ]
# 
# gene_status <- read.delim(
#    paste0(base_dir,'/passengers/processed/metadata/gene_status/has_gene_mut.txt.gz'),
#    check.names=F
# )
# df$POLE_mut <- gene_status[match(df$sample, rownames(gene_status)),'POLE']
# df[order(df$contrib, decreasing=T),]

## Add signature etiologies ================================
sig_etiologies <- read.delim(paste0(wd,'/sig_etiologies.txt'))
sig_metadata$etiology <- NULL
sig_metadata <- sig_metadata[!is.na(sig_metadata$sig_name),] ## Rm rank 0 signatures

# sig_metadata[
#   match(sig_metadata$sig_name, sig_etiologies$sig_name) |> is.na() |> which()
# ,]

sig_etiologies_matched <- sig_etiologies[
  match(sig_metadata$sig_name, sig_etiologies$sig_name),
  c('etiology','etiology_group_1','etiology_group_2')
]

sig_metadata <- cbind(sig_metadata, sig_etiologies_matched)
write.table(
  sig_metadata,
  paste0(nmf_dir,'/sig_contrib/fit_lsq.post_processed/sig_metadata.post_processed.with_etiologies.txt.gz'),
  sep='\t',row.names=F, quote=F
)






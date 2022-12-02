############### Extract the mutation context profile for SVs ############
# author: luan

### Description
# This script extracts the SV context matrices based downsampled Hartwig VCF files.

### Input
# downsampled_sv_paths.txt.gz (manifest file)

### Output
# Extracted mutational context matrices. One matrix per downsampled Hartwig sample for the SVs.
# Output matrices are placed in separate folders.

### Usage
# sbatch extract_sv_type_load.sh

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

## Args ================================
write.tsv <- function(...){ write.table(..., sep='\t', quote=F) }

args <- commandArgs(trailingOnly=TRUE)
sample_id <- args[1]
linx_vis_sv_data <- args[2]
output_dir <- args[3]

## Main ================================
linx_svs <- read.delim(linx_vis_sv_data)
linx_svs$sample <- as.factor(sample_id)

## DEL, DUP --------------------------------
countDelDupByLen <- function(linx_svs, bin.breaks=c(0,10000,Inf), resolved.types=c('DEL','DUP'), output='contexts'){
   
   if(!(output %in% c('contexts','raw'))){
      stop("`output` must be 'contexts' or 'raw'")
   }
   
   ## --------------------------------
   bins <- levels(
      cut(0, bin.breaks, right=FALSE, include.lowest=FALSE)
   )
   bins <- unlist(lapply(resolved.types, function(i){
      paste0(i,'_',bins)
   }))
   
   template <- matrix(
      rep(0, length(bins)),
      nrow=1
   )
   colnames(template) <- bins
   
   ## --------------------------------
   df <- subset(
      linx_svs, 
      ResolvedType %in% resolved.types, 
      c(sample, ClusterId, ResolvedType, ChrStart, ChrEnd, PosStart, PosEnd)
   )
   
   if(nrow(df)==0){
      return(template)
   }
   
   df$ClusterId <- paste0(as.integer(df$sample),'_',df$ClusterId)
   
   ## --------------------------------
   ## Select intrachromosomal SVs
   df$is_same_chrom <- df$ChrStart==df$ChrEnd
   
   intrachrom_clusters <- unlist(lapply(split(df$is_same_chrom, df$ClusterId), all))
   intrachrom_clusters <- names(intrachrom_clusters)[intrachrom_clusters]
   df <- df[df$ClusterId %in% intrachrom_clusters,]
   
   NULL -> df$ChrStart -> df$ChrEnd -> df$is_same_chrom -> intrachrom_clusters
   
   ## Flatten clusters of SVs into one row
   df$is_start_fragment <- !duplicated(df$ClusterId, fromLast=F)
   df$is_end_fragment <- !duplicated(df$ClusterId, fromLast=T)
   
   df_flat <- data.frame(
      sample=df$sample[df$is_start_fragment],
      ClusterId=df$ClusterId[df$is_start_fragment],
      ResolvedType=df$ResolvedType[df$is_start_fragment],
      start=df$PosStart[df$is_start_fragment],
      end=df$PosEnd[df$is_end_fragment]
   )
   
   ## Calculate SV lengths
   # ## start coords are 1-based
   # all(df_flat$start <= df_flat$end)
   # df_flat[df_flat$start == df_flat$end,]
   df_flat$length <- 1 + df_flat$end - df_flat$start
   
   ## Bin SVs by length
   df_flat$length_bin <- cut(df_flat$length, bin.breaks, right=FALSE, include.lowest=FALSE)
   
   ## Remove SVs not falling into the bins specified in bin.breaks
   df_flat <- df_flat[!is.na(df_flat$length_bin),]
   
   ## Make SV type/length contexts
   df_flat$context <- paste0(df_flat$ResolvedType, '_', df_flat$length_bin)
   df_flat <- df_flat[order(df_flat$ResolvedType, df_flat$length_bin),]
   df_flat$context <- factor(df_flat$context, unique(df_flat$context))
   
   if(output=='raw'){ return(df_flat) }
   
   ## Make context matrix -------------------------------
   tab <- table(df_flat$context)
   
   m <- template[1,]
   m[names(tab)] <- tab
   m <- unclass(m)
   m <- matrix(m, nrow=1)
   colnames(m) <- colnames(template)
   
   return(m)
}

dels <- countDelDupByLen(linx_svs, resolved.types='DEL', bin.breaks=c(0,10000,Inf))
dups <- countDelDupByLen(linx_svs, resolved.types='DUP', bin.breaks=c(0,10000,Inf))

## Complex len --------------------------------
countComplexByLen <- function(linx_svs, bin.breaks=c(0,20,Inf), include.inv.trans=T){
   
   ## --------------------------------
   bins <- levels(
      cut(0, bin.breaks, right=FALSE, include.lowest=FALSE)
   )
   bins <- paste0('COMPLEX_',bins)
   
   template <- matrix(
      rep(0, length(bins)),
      nrow=1
   )
   colnames(template) <- bins
   
   ## --------------------------------
   df <- linx_svs[,c('sample', 'ClusterId', 'ResolvedType')]
   
   if(!include.inv.trans){
      df <- df[df$ResolvedType=='COMPLEX',]
   } else {
      df <- df[df$ResolvedType=='COMPLEX' | grepl('INV|TRANS',df$ResolvedType),]
   }
   
   if(nrow(df)==0){
      return(template)
   }
   
   ## Bin complex clusters by length --------------------------------
   complex_len <- with(df,{
      out <- aggregate(
         ClusterId,
         list(sample=sample, ClusterId=ClusterId),
         function(x){ length(x) }
      )
      colnames(out)[length(out)] <- 'len'
      return(out)
   })
   
   complex_len$len_bin <- cut(
      complex_len$len, 
      bin.breaks, 
      right=FALSE, 
      include.lowest=FALSE
   )
   
   ## Make context matrix -------------------------------
   tab <- table(complex_len$len_bin)
   names(tab) <- paste0('COMPLEX_',names(tab))
   
   m <- template[1,]
   m[names(tab)] <- tab
   m <- unclass(m)
   m <- matrix(m, nrow=1)
   colnames(m) <- colnames(template)
   
   return(m)
}

complex_svs <- countComplexByLen(linx_svs, bin.breaks=c(0,20,Inf), include.inv.trans=T)

## SV type load ================================
sv_type_load <- (function(){
   df <- data.frame(
      #SV_LOAD=sv_load,
      dels,
      dups,
      complex_svs,
      LINE=length(unique(linx_svs$ClusterId[linx_svs$ResolvedType=='LINE'])),
      check.names=F
   )
   
   df <- cbind(sample=sample_id, SV_LOAD=rowSums(df), df)
   rownames(df) <- NULL
   
   return(df)
})()

# write to disk
write.table(sv_type_load, output_dir, sep='\t', quote=F, row.names=F)

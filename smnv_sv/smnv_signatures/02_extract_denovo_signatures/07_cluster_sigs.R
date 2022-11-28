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
wd <- paste0(base_dir,'/passengers/processed/sigs_denovo/extractions/12_fixed_seeds/')

## Dependencies --------------------------------
#source(paste0(base_dir,'/passengers/analysis/common_functions/common_functions.R'))
#devtools::load_all(paste0(path_prefix,'/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor'))
library(mutSigExtractor)
library(ggplot2)

## Load sample metadata ================================
sample_metadata <- read.delim(
  paste0(base_dir,'/passengers/processed/metadata/main/metadata_20220929_1720.txt.gz'),
  na.strings=c('NA','#N/A')
)

## Primary/metastatic
sample_metadata$is_hmf_sample <- sample_metadata$cohort=='Hartwig'

## HRD
sample_metadata$is_hrd <- sample_metadata$hr_status=='HR_deficient'
sample_metadata$is_hrd[sample_metadata$hr_status=='cannot_be_determined'] <- NA

## MSI
sample_metadata$is_msi <- sample_metadata$msi_status=='MSI'
sample_metadata$is_msi[is.na(sample_metadata$msi_status)] <- NA

##
getSampleMetadata <- function(sample.names, keys=NULL, drop=T){
  #sample.names=c('CPCT02030502','CPCT02180026')
  if(is.null(keys)){ keys <- colnames(sample_metadata) }
  sample_metadata[match(sample.names, sample_metadata$sample_id), keys, drop=drop]
}

## Load signature data ================================
## Profiles --------------------------------
loadSigProfiles <- function(dir){
   #dir=paste0(wd,'/sig_profiles/')
   paths <- list.files(
      dir,
      pattern='[.]txt$', full.names=T
   )
   l <- lapply(paths, read.delim, check.names=F)
   names(l) <- sub('[.]txt$','',basename(paths))
   return(l)
}

sig_profiles <- loadSigProfiles(paste0(wd,'/sig_profiles/'))

## Contribs --------------------------------
sig_contribs <- (function(){
   paths <- list.files(
      paste0(wd,'/sig_contrib/fit_lsq.by_tissue_type/'),
      pattern='[.]txt[.]gz$', full.names=T
   )
   l <- lapply(paths, read.delim, check.names=F)
   names(l) <- sub('[.]txt[.]gz$','',basename(paths))
   return(l)
})()

## Make unique names for signatures
sig_contribs <- (function(){
   l_names <- names(sig_contribs)
   l <- lapply(names(sig_contribs), function(i){
      #i='Biliary'
      m <- sig_contribs[[i]]
      colnames(m) <- paste0( i,'.',sub('[.].*$','',colnames(m)))
      return(m)
   })
   names(l) <- l_names
   return(l)
})()


## Parse sig names ================================
sig_metadata <- reshape2::melt( lapply(sig_profiles, colnames) )
colnames(sig_metadata) <- c('denovo_name_raw','nmf_group')

## Parse extraction name
sig_metadata <- (function(){
   l <- strsplit(sig_metadata$nmf_group,'[.]')
   
   #l <- lapply(l, function(i){ as.data.frame(t(i)) })
   df <- do.call(rbind, l)
   colnames(df) <- c('tissue_type','mut_type','denovo_sigs_in_group')
   cbind(sig_metadata, df)
})()
sig_metadata$nmf_group <- NULL

## Parse denovo sig names
sig_metadata <- (function(){
   l <- strsplit(sig_metadata$denovo_name_raw,'[.]')
   l <- lapply(l, function(i){ as.data.frame(t(i)) })
   df <- do.call(plyr::rbind.fill, l)
   colnames(df) <- c('denovo_name_short','matched_ref_name','etiology','matched_ref_cos_sim')
   cbind(sig_metadata, df)
})()
sig_metadata$matched_ref_cos_sim <- as.numeric(sig_metadata$matched_ref_cos_sim)/100

##
sig_metadata$denovo_name_uniq  <- with(sig_metadata, paste0(tissue_type,'.',denovo_name_short))

## Remove unneeded columns
NULL -> sig_metadata$denovo_sigs_in_group -> sig_metadata$denovo_name_short

## Use most up to date etiologies
mutSigExtractor.sig_metadata <- read.delim(paste0(path_prefix,'/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor/inst/sigs_v3.2_metadata.txt'))
sig_metadata$etiology <- mutSigExtractor.sig_metadata$sig_etiology[ match(sig_metadata$matched_ref_name, mutSigExtractor.sig_metadata$sig_name) ]
sig_metadata$etiology[is.na(sig_metadata$etiology)] <- 'Unknown'

##
getSigMetadata <- function(denovo.name.uniq, keys=NULL, drop=T){
  if(is.null(keys)){ keys <- colnames(sig_metadata) }
  sig_metadata[match(denovo.name.uniq, sig_metadata$denovo_name_uniq), keys, drop=drop]
}

## Novel signatures ================================
mergeSigProfileMatrices <- function(sig.profiles, use.denovo.uniq.name=T, by='mut.type'){
   #l=sig_profiles
   
   if(!(by %in% c('mut.type','tissue.type'))){
      stop("`by` must be 'mut.type' or 'tissue.type'")
   }
   
   ##
   l <- sig.profiles
   l <- l[!grepl('^Other',names(l))]  ## Remove Other group
   extraction_names <- names(l)
   
   ##
   tissue_types <- gsub('[.].+$','',extraction_names)
   l <- lapply(1:length(l), function(i){
      #i=1
      m <- l[[i]]
      if(use.denovo.uniq.name){ colnames(m) <- sub('[.].+$','',colnames(m)) }
      colnames(m) <- paste0(tissue_types[[i]],'.',colnames(m))
      return(m)
   })
   names(l) <- extraction_names
   
   mut_types_uniq <- c('SBS','DBS','ID')
   
   ##
   if(by=='tissue.type'){
      names(l) <- stringr::str_extract(names(l), '[:alnum:]+[.][:alnum:]+')
      return(l)
   }
   
   ##
   if(by=='mut.type'){
      
      mut_types <- sapply(strsplit(extraction_names,'[.]'),`[[`,2)
      out <- lapply(mut_types_uniq,function(i){
         #i='SBS'
         l[mut_types==i] |> unname() |> (\(x) do.call(cbind, x))()
      })
      names(out) <- mut_types_uniq
   }
   
   return(out)
}

sig_profiles.merged <- mergeSigProfileMatrices(sig_profiles)


## Functions --------------------------------
clusterDenovoSigs <- function(profiles, k.override=NULL){
  
  if(F){
    profiles=l_profiles$ID
    k.override=7
    
    profiles=l_profiles$SBS
    k.override=27
  }
  
  ## Hierachical clustering
  cos_sim <- mutSigExtractor::cosSim(profiles, profiles, all.vs.all=T, by.row=F)
  distances <- 1-cos_sim
  #distances <- 1-cor(profiles)
  hc <- hclust(as.dist(distances))
  
  ## Determine optimal number of clusters
  k_range <- 2:(nrow(distances)-1)
  avg_sil_width <- sapply(k_range, function(i){
    clusters <- cutree(hc, k=i)
    sil <- cluster::silhouette(clusters, distances)
    mean(sil[,3])
  })
  
  best_k <- k_range[which.max(avg_sil_width)]
  summ <- data.frame(k=k_range, avg_sil_width, is_best=FALSE)
  
  ## Select k
  summ$is_chosen_k <- FALSE
  chosen_k <- best_k
  if(is.null(k.override)){
    summ$is_chosen_k[summ$k==best_k] <- TRUE
  } else {
    if(k.override>max(k_range)){
      chosen_k <- max(k_range)
      warning('`k` is larger than the max k: ',max(k_range),'. Max will be chosen instead')
    } else {
      chosen_k <- k.override
    }
    summ$is_chosen_k[summ$k==chosen_k] <- TRUE
  }
  
  clusters_final <- cutree(hc, k=chosen_k)
  
  ## Calculate within cluster cos sims   
  cluster_members <- split(names(clusters_final), clusters_final)
  mean_within_clust_dist <- sapply(cluster_members, function(i){
    #i=cluster_members[[2]]
    ## With this method, one dissimilar profile can skew the measure
    #mean(as.dist(distances[i,i])) 

    ## This method is less sensitive to outliers
    distances_ss <- distances[i,i]
    if(length(distances_ss)==1){ return(NA) }
    diag(distances_ss) <- NA
    #median(rowMeans(distances_ss, na.rm=T))
    mean(rowMeans(distances_ss, na.rm=T))
  })
  
  ## Clustering summary table
  summ <- getSigMetadata(names(clusters_final))
  summ$cluster <- unname(clusters_final)
  summ <- summ[order(summ$cluster),]
  summ$silhouette_score <- cluster::silhouette(clusters_final, distances)[,3]
  summ$mean_within_clust_cossim <- 1-mean_within_clust_dist[summ$cluster]
  
  class(summ) <- c(class(summ), 'clust_summ')
  
  ##
  sil_out <- data.frame(
    k=k_range, 
    avg_sil_width, 
    is_best_k=k_range==best_k,
    is_chosen_k=k_range==chosen_k
  )
  class(sil_out) <- c(class(sil_out), 'sil')
  
  ##
  list(summ=summ, silhouette=sil_out)
}
#clust_out <- clusterDenovoSigs(profiles)

plotSilhouette <- function(x){
  #profiles=l_profiles[['ID']]
  require(ggplot2)
  
  if(inherits(x,'sil')){
    pd <- x
  } else {
    pd <- clusterDenovoSigs(x)$silhouette
  }
  
  chosen_k <- pd$k[pd$is_chosen_k]
  best_k <- pd$k[pd$is_best_k]
  
  ggplot(pd, aes(x=k, y=avg_sil_width)) +
    geom_point() +
    geom_line() +
    
    geom_vline(xintercept=best_k, color='blue', linetype='dotted') +
    annotate(geom='text', x=best_k, y=max(pd$avg_sil_width), vjust=1, hjust=0, color='blue', label=paste0('Best k = ',best_k)) +
    
    geom_vline(xintercept=chosen_k, color='red', linetype='dashed') +
    annotate(geom='text', x=chosen_k, y=min(pd$avg_sil_width), vjust=1, hjust=0, color='red', label=paste0('Chosen k = ',chosen_k)) +
    
    labs(x='No. of clusters', y='Average silhouette width') +
    theme_bw() +
    theme(
      panel.grid.minor=element_blank()
    )
}

clusterReport <- function(summ, profiles, plot.title=NULL){
  #summ=subset(clust_out$summ, cluster==1)
  
  require(ggplot2)
  
  ## Contribs --------------------------------
  ## Get contribs
  tmp <- summ[,c('tissue_type','denovo_name_uniq')]
  pd_contribs <- lapply(1:nrow(tmp), function(i){
    #print(i)
    #i=1
    tissue_type <- tmp$tissue_type[i]
    denovo_name_uniq <- tmp$denovo_name_uniq[i]
    df <- sig_contribs[[tissue_type]][denovo_name_uniq]
    data.frame(sample=rownames(df), signature=colnames(df), contrib=df[,1])
  })
  pd_contribs <- do.call(rbind, pd_contribs)
  pd_contribs$signature <- factor(pd_contribs$signature, unique(pd_contribs$signature))
  
  ## Calc ordered x position
  pd_contribs_split <- split(pd_contribs, pd_contribs$signature)
  pd_contribs <- do.call(rbind, lapply(pd_contribs_split, function(j){
    #j=pd_sig_contribs_split[[1]]
    j <- j[order(j$contrib),]
    j$index <- 1:nrow(j)
    #j$xpos <- (1:nrow(j)) / nrow(j)
    return(j)
  })); rownames(pd_contribs) <- NULL
  rm(pd_contribs_split)
  
  ## Add cohort
  pd_contribs$cohort <- getSampleMetadata(pd_contribs$sample, 'cohort')
  
  ## Calc medians
  pd_sig_contrib_med <- with(pd_contribs,{
    aggregate(contrib+1, list(signature=signature), median)
  })
  colnames(pd_sig_contrib_med)[2] <- 'contrib'
  med_pos <- with(pd_contribs,{
    agg <- aggregate(index, list(signature), max)
    out <- structure(agg$x, names=as.character(agg[,1]))
    out * 0.5
  })
  
  p_contribs <- ggplot(pd_contribs, aes(x=index, y=contrib+1)) +
    facet_grid(signature~., scales='free_y', space='free_x') +
    geom_point(aes(color=cohort), size=0.2) +
    scale_color_manual(values=c(Hartwig='#6B479A', PCAWG='#F48134')) +
    
    geom_point(data=pd_sig_contrib_med, aes(y=contrib), x=med_pos, color='black', shape=95, size=10) +
    scale_y_continuous(trans='log10', name='Contribution + 1') +
    xlab('Sample rank') +
    theme_bw() +
    theme(
      #strip.text.y=element_text(angle=0),
      strip.text.y=element_blank(),
      panel.grid.major.x=element_blank(),
      panel.grid.minor.x=element_blank(),
      panel.grid.minor.y=element_blank(),
      legend.position='top'
      #axis.title.x=element_blank(),
      #axis.text.x=element_blank(),
      #axis.ticks.x=element_blank()
    )
  
  ## Profiles --------------------------------
  # trimText <- function(v, width=25){
  #   #v=rownames(sig_profiles)
  #   #width=20
  #   
  #   v_short <- substring(v, 1, width)
  #   which_long <- nchar(v) > width
  #   v_short[which_long] <- paste0(v_short[which_long],'...')
  #   
  #   return(v_short)
  # }
  
  sel_sigs <- summ$denovo_name_uniq
  profiles_ss <- profiles[,sel_sigs,drop=F]
  
  p_profiles <-
    mutSigExtractor::plotContexts(
      t(profiles_ss), group=colnames(profiles_ss), 
      y.axis.var.scale=T, force.group.labels=T
    ) +
    #ggtitle(paste0('Rank: ', rank.name)) +
    ylab('Probability') +
    theme(
      strip.text.y=element_text(angle=0)
    )
  
  if(!is.null(plot.title)){
    p_profiles <- p_profiles + ggtitle(plot.title)
  }
  
  ## Merge --------------------------------
  cowplot::plot_grid(
    p_contribs, p_profiles,
    nrow=1, align='h', axis='tblr', rel_widths=c(0.2, 1)
  )
}

clusterReportWrapper <- function(profiles, clust_out=NULL){
  #profiles=sig_profiles.merged[['ID']]
  if(is.null(clust_out)){
    clust_out <- clusterDenovoSigs(profiles)
  }
  uniq_clusters <- unique(clust_out$summ$cluster)
  plots <- lapply(uniq_clusters, function(i){
    #i=1
    #print(i)
    summ <- subset(clust_out$summ, cluster==i)
    plot_title <- paste0(
      'Cluster ', i,
      '. Mean within cluster cos sim=',round(summ$mean_within_clust_cossim[1],2))
    clusterReport(summ, profiles, plot.title=plot_title)
  })
}

## Preliminary clustering to determin silhouette scores --------------------------------
##
MIN_COS_SIM <- 0.85
low_cossim_sigs <- subset(sig_metadata, matched_ref_cos_sim<MIN_COS_SIM, denovo_name_uniq, drop=T)

MUT_TYPES <- c('SBS','DBS','ID')
l_profiles <- lapply(MUT_TYPES, function(i){
  profiles <- sig_profiles.merged[[i]]
  profiles <- profiles[,colnames(profiles) %in% low_cossim_sigs,drop=F]
  return(profiles)
})
names(l_profiles) <- MUT_TYPES

## Manually choose k
if(F){
  plotSilhouette(l_profiles$SBS)
  plotSilhouette(l_profiles$DBS)
  plotSilhouette(l_profiles$ID)
}

## Clustering
cluster_out <- list(
  SBS=clusterDenovoSigs(l_profiles$SBS, k.override=27),
  DBS=clusterDenovoSigs(l_profiles$DBS, k.override=7),
  ID=clusterDenovoSigs(l_profiles$ID, k.override=8)
)

cluster_sil <- lapply(cluster_out,`[[`,'silhouette')
# plotSilhouette(cluster_sil$SBS)
# plotSilhouette(cluster_sil$DBS)
# plotSilhouette(cluster_sil$ID)

cluster_summ <- lapply(cluster_out,`[[`,'summ')
cluster_summ <- do.call(rbind, cluster_summ); rownames(cluster_summ) <- NULL

cluster_summ$cluster_name <- with(cluster_summ, paste0(
  mut_type,'_denovo_',cluster
))

## Plot sig profiles --------------------------------
## Denovo clusters
if(WRITE_OUTPUT){
  dir.create(paste0(wd,'/scripts/07_clustering/'), showWarnings=F)
  
  l_cluster_reports <- lapply(MUT_TYPES, function(i){
    clusterReportWrapper(profiles=l_profiles[[i]], clust_out=cluster_out[[i]])
  })
  names(l_cluster_reports) <- MUT_TYPES
  
  for(i in MUT_TYPES){
    #i='SBS'
    message('Plotting novel sigs for: ',i)
    p_silhouette <- plotSilhouette(cluster_sil[[i]])
    l_plots <- c(
      list(p_silhouette),
      l_cluster_reports[[i]]
    )
    
    out_path <- paste0(wd,'/scripts/07_clustering/novel_sigs.',i,'.pdf')
    pdf(out_path, 12, 7)
    for(j in l_plots){ plot(j) }
    dev.off()
  }
}

## All signatures
if(WRITE_OUTPUT){
  plotContextsWrapper <- function(m_sig_profiles, out.path=NULL){
    #m_sig_profiles=sig_profiles.merged$ID
    m_sig_profiles <- t(m_sig_profiles)
    l_plots <- lapply(rownames(m_sig_profiles), function(i){
      #i='Biliary.ID_A'
      i_m <- m_sig_profiles[i,,drop=F]
      plotContexts(i_m) + 
        ggtitle(i) +
        theme(panel.grid=element_blank())
    })
    
    if(is.null(out.path)){ return(l_plots) }
    
    #out.path='/home/lnguyen/hpc/cuppen/projects/P0025_PCAWG_HMF/passengers/processed/sigs_denovo/extractions/10_rm_no_consent/scripts/07_clustering/all_sigs.ID.pdf'
    pdf(out.path, 10, 2.5)
    for(i in l_plots){ plot(i) }
    dev.off()
  }
  
  plotContextsWrapper(sig_profiles.merged$SBS, paste0(wd,'/scripts/07_clustering/all_sigs.SBS.pdf'))
  plotContextsWrapper(sig_profiles.merged$DBS, paste0(wd,'/scripts/07_clustering/all_sigs.DBS.pdf'))
  plotContextsWrapper(sig_profiles.merged$ID, paste0(wd,'/scripts/07_clustering/all_sigs.ID.pdf'))
}

##
if(WRITE_OUTPUT){
  write.table(
    cluster_summ,
    paste0(wd,'/scripts/07_clustering/novel_denovo_clusters.txt'),
    sep='\t',quote=F,row.names=F
  )
}

## Reannotate sig metadata --------------------------------
sig_metadata$sig_name <- cluster_summ$cluster_name[ match(sig_metadata$denovo_name_uniq, cluster_summ$denovo_name_uniq) ]
denovo_clust_indexes <- grepl('denovo',sig_metadata$sig_name)
sig_metadata$is_denovo_clust <- denovo_clust_indexes

## For the remaining (high cos sim sigs), assign ref sig name
sig_metadata$sig_name[!denovo_clust_indexes] <- sig_metadata$matched_ref_name[!denovo_clust_indexes]

## Add cos sim to old signatures --------------------------------
current_vs_old_sig_cos_sims <- (function(){
   sig_profiles.current <- mergeSigProfileMatrices(sig_profiles, by='tissue.type')
   
   sig_profiles.old <- loadSigProfiles(
      paste0(base_dir,'/passengers/processed/sigs_denovo/extractions/09_fixed_smnvs/sig_profiles/')
   )
   sig_profiles.old <- mergeSigProfileMatrices(sig_profiles.old, by='tissue.type')
   
   l <- lapply(names(sig_profiles.current), function(i){
      #i='Biliary.SBS'
      #i='Skin.DBS'
      #print(i)
      m_current <- t(sig_profiles.current[[i]])
      if(any(grepl('0',rownames(m_current)))){ return(NULL) }
      
      m_old <- t(sig_profiles.old[[i]])
      #rownames(m_old) <- paste0(rownames(m_old), '.old')
      
      cos_sims <- mutSigExtractor::cosSim(m_current, m_old, all.vs.all=T)
      cos_sims <- as.matrix(cos_sims)
      old_sig.name <- rownames(cos_sims)[ apply(cos_sims, 2, which.max) ]
      old_sig.max_cos_sim <- apply(cos_sims, 2, max)
      
      data.frame(
         sig_name.current=rownames(m_current),
         sig_name.old=old_sig.name,
         cos_sim=old_sig.max_cos_sim,
         row.names=NULL
      )
   })
   do.call(rbind, l)
})()

sig_metadata$preprint_matched_sig_name <- current_vs_old_sig_cos_sims$sig_name.old[match(sig_metadata$denovo_name_uniq, current_vs_old_sig_cos_sims$sig_name.current)]
sig_metadata$preprint_sig_cos_sim <- current_vs_old_sig_cos_sims$cos_sim[match(sig_metadata$denovo_name_uniq, current_vs_old_sig_cos_sims$sig_name.current)]
sig_metadata$preprint_is_same_sig <- sig_metadata$denovo_name_uniq==sig_metadata$preprint_matched_sig_name

## Save original `sig_metadata` --------------------------------
sig_metadata.pre <- sig_metadata
if(WRITE_OUTPUT){ 
  write.table(
    sig_metadata, paste0(wd,'/scripts/07_clustering/sig_metadata.pre.txt'),
    sep='\t', quote=F, row.names=F
  )
}

## For denovo clusts, set non-applicable column values to NA
sig_metadata$etiology[denovo_clust_indexes] <- 'Unknown'
sig_metadata$matched_ref_name[denovo_clust_indexes] <- NA
sig_metadata$matched_ref_cos_sim[denovo_clust_indexes] <- NA

## Manually assign signatures  --------------------------------
if(F){
  cancerTypesWithMissingSigs <- function(sig.name='SBS1', sig.metadata=sig_metadata, simplify=TRUE){
    #sig.metadata=sig_metadata.pre
    #sig.name=c('SBS5','SBS40')
    
    sig.metadata$sig_exists <- sig.metadata$sig_name %in% sig.name
    
    l <- split(sig.metadata$sig_exists, sig.metadata$tissue_type)
    v <- unlist(lapply(l, function(i){ any(i) }))
    v <- !v
    if(!simplify){ return(v) }
    names(v)[v]
  }
  
  ##SBS1
  cancerTypesWithMissingSigs('SBS1')
  
  ## SBS5/40
  cancerTypesWithMissingSigs(c('SBS5','SBS40'), sig_metadata.pre)
  df <- subset(sig_metadata.pre, matched_ref_name %in% c('SBS5','SBS40') & matched_ref_cos_sim<MIN_COS_SIM)
  df[c('denovo_name_uniq','matched_ref_name','matched_ref_cos_sim')]
  subset(sig_metadata.pre, tissue_type=='Kidney')
  
  ## Misc
  subset(sig_metadata.pre, tissue_type=='Breast')
  subset(sig_metadata.pre, tissue_type=='Prostate')
  subset(sig_metadata.pre, tissue_type=='Thyroid')
  
  ## ----
}

if(WRITE_OUTPUT){
  cosSimWrapper <- function(profiles.denovo=sig_profiles.merged$SBS, profiles.ref=SBS_SIGNATURE_PROFILES_V3, top.n=5){
    #profiles.denovo=sig_profiles.merged$SBS
    #profiles.ref=SBS_SIGNATURE_PROFILES_V3
    
    l_cos_sims <- lapply(1:ncol(profiles.denovo), function(i){
      #i=1
      i_m <- profiles.denovo[,i,drop=F]
      cos_sims <- cosSim(i_m, profiles.ref, all.vs.all=T, by.row=F)
      
      out <- data.frame(
        denovo_name_uniq=colnames(profiles.denovo)[i],
        ref_name=rownames(cos_sims),
        cos_sim=cos_sims[,1],
        row.names=NULL
      )
      out <- out[order(out$cos_sim, decreasing=T),]
      out$index <- 1:nrow(out)
      return(out)
    })
    
    out <- do.call(rbind, l_cos_sims)
    out[out$index<=top.n,]
  }
  
  write.table(
    cosSimWrapper(), paste0(wd,'/scripts/07_clustering/cos_sim_to_cosmic.txt'),
    sep='\t', quote=F, row.names=F
  )
}

sig_manual_assignments <- c(
  ## SBS1 is expected in all tissue types
  Cervix.SBS_D='SBS1', ## Visually checked sig profile
  Liver.SBS_H='SBS1', ## Visually checked sig profile
  Lymphoid.SBS_B='SBS1', ## Visually checked sig profile
  Thyroid.SBS_C='SBS1', ## Visually checked sig profile
  
  ## SBS5/40 is expected in all tissue types
  Biliary.SBS_C='SBS5', ## SBS5: 0.74
  Breast.SBS_G='SBS5', ## Assignment from Arne
  Breast.SBS_K='SBS40', ## Assignment from Arne
  Cervix.SBS_C='SBS6', ## Assignment from Arne
  Cervix.SBS_F='SBS5', ## Assignment from Arne
  Cervix.SBS_E='SBS8', ## Assignment from Arne
  CNS.SBS_H='SBS40', ## SBS40: 0.71
  CNS.SBS_F='SBS5', ## Assignment from Arne
  Colorectum.SBS_E='SBS40', ## SBS5: 0.72
  Colorectum.SBS_I='SBS5', ## Assignment from Arne
  Gastric.SBS_H='SBS5', ## Assignment from Arne
  HeadAndNeck.SBS_K='SBS40', ## Assignment from Arne
  #Liver.SBS_G='SBS40', ## SBS3: 0.86, SBS40: 0.79, Visually checked sig profile
  Lung.SBS_H='SBS40', ## SBS40: 0.81
  Lymphoid.SBS_C='SBS5', ## SBS5: 0.73 ## No visual match
  Ovary.SBS_C='SBS5', ## SBS5: 0.73
  #Pancreas ## No match
  Pancreas_NET.SBS_C='SBS5', ## SBS12: 0.87, SBS5: 0.85
  Prostate.SBS_C='SBS5', ## SBS5: 0.81
  Prostate.SBS_F='SBS40', ## SBS40: 0.80
  Urothelial.SBS_F='SBS5', ## SBS40: 0.76
  Uterus.SBS_F='SBS40', ## SBS40: 0.83
  
  ## Misc
  Breast.SBS_J='SBS44', ## Assignment from Arne
  Colorectum.SBS_D='SBS30', ## Assignment from Arne
  Prostate.SBS_G='SBS26' ## SBS12: 0.9, SBS26: 0.9
)

sig_manual_assignments.cos_sims <- sapply(names(sig_manual_assignments), function(i){ ## Recalculate matched_ref_cos_sim
  #i='Pancreas.NET_C'
  #print(i)
  denovo_sig <- sig_profiles.merged$SBS[[i]]
  target_ref_sig <- sig_manual_assignments[[i]]
  
  if(!(target_ref_sig %in% colnames(mutSigExtractor::SBS_SIGNATURE_PROFILES_V3))){
    warning(i,': ',target_ref_sig, ' not in `SBS_SIGNATURE_PROFILES_V3`')
    return(NA)
  }
  
  ref_sig <- mutSigExtractor::SBS_SIGNATURE_PROFILES_V3[, target_ref_sig ]
  round(cosSim(denovo_sig, ref_sig),2)
})

sig_metadata$is_manually_annotated <- sig_metadata$denovo_name_uniq %in% names(sig_manual_assignments)

sig_manual_assignments.indexes <- match(names(sig_manual_assignments), sig_metadata$denovo_name_uniq)
sig_metadata$sig_name[sig_manual_assignments.indexes] <- sig_manual_assignments

#subset(sig_metadata, is_denovo_clust & sig_name %in% c('SBS5','SBS40'))
sig_metadata$is_denovo_clust[sig_manual_assignments.indexes] <- FALSE
sig_metadata$matched_ref_cos_sim[sig_manual_assignments.indexes] <- sig_manual_assignments.cos_sims
sig_metadata$matched_ref_cos_sim[sig_manual_assignments.indexes] <- sig_manual_assignments.cos_sims

## Retrieve or calculate mean_within_clust_cossim
sig_metadata$mean_within_clust_cossim <- cluster_summ$mean_within_clust_cossim[ match(sig_metadata$denovo_name_uniq, cluster_summ$denovo_name_uniq) ]
sig_metadata$mean_within_clust_cossim <- round(sig_metadata$mean_within_clust_cossim, 2)
sig_metadata$mean_within_clust_cossim[sig_manual_assignments.indexes] <- NA

## Edit colnames of sig_contribs matrix --------------------------------
sig_contribs.renamed <- lapply(sig_contribs, function(i){
  #i=sig_contribs$Thyroid
  colnames(i) <- sig_metadata$sig_name[ match(colnames(i), sig_metadata$denovo_name_uniq) ]
  #colnames(i) <- getSigMetadata(colnames(i),'sig_name')
  #i[,naturalsort::naturalorder(colnames(i)),drop=F]
  return(i)
})

##
sig_contribs.renamed <- do.call(plyr::rbind.fill, sig_contribs.renamed)
rownames(sig_contribs.renamed) <- unlist(lapply(sig_contribs, rownames), use.names=F)
sig_contribs.renamed[is.na(sig_contribs.renamed)] <- 0

## Remove empty placeholder signatures
sig_contribs.renamed <- sig_contribs.renamed[,colnames(sig_contribs.renamed)!='NA']

## Sort signatures alphabetically and by mut type
sig_contribs.renamed <- sig_contribs.renamed[,naturalsort::naturalorder(colnames(sig_contribs.renamed)),drop=F]
sig_contribs.renamed <- lapply(MUT_TYPES, function(i){ 
  sig_contribs.renamed[,grep(i,colnames(sig_contribs.renamed)),drop=F] 
})
sig_contribs.renamed <- do.call(cbind, sig_contribs.renamed)

## Calculate mean_within_clust_cossim for cosmic sigs --------------------------------
cosmic_matched_sigs.mean_within_clust_cossim <- (function(){
  cosmic_matched_sigs <- subset(sig_metadata, !is_denovo_clust, c(denovo_name_uniq, sig_name, mut_type))
  cosmic_matched_sigs <- unique(cosmic_matched_sigs)
  cosmic_matched_sigs <- cosmic_matched_sigs[cosmic_matched_sigs$sig_name!='NA',]

  cosmic_matched_sigs_uniq <- sort(unique(cosmic_matched_sigs$sig_name))
  cosmic_matched_sigs_uniq_type <- cosmic_matched_sigs$mut_type[ match(cosmic_matched_sigs_uniq, cosmic_matched_sigs$sig_name) ]
  cosmic_matched_sigs_uniq <- split(cosmic_matched_sigs_uniq, cosmic_matched_sigs_uniq_type)

  l <- lapply(names(cosmic_matched_sigs_uniq), function(i){
    #i='DBS'
    i.cosmic_matched_sigs_uniq <- cosmic_matched_sigs_uniq[[i]]
    i.sig_profiles.merged <- sig_profiles.merged[[i]]

    i.l <- lapply(i.cosmic_matched_sigs_uniq, function(j){
      #j='DBS1'
      #print(j)
      j.denovo_name_uniq <- cosmic_matched_sigs$denovo_name_uniq[ cosmic_matched_sigs$sig_name==j ]
      m <- i.sig_profiles.merged[,colnames(i.sig_profiles.merged) %in% j.denovo_name_uniq, drop=F]

      if(ncol(m)==1){ return(NA) }

      distances <- 1-mutSigExtractor:::cosSim.matrix(m, m, all.vs.all=T)
      diag(distances) <- NA
      #1 - median(rowMeans(distances, na.rm=T))
      1 - mean(rowMeans(distances, na.rm=T))
    })
    names(i.l) <- i.cosmic_matched_sigs_uniq
    i.v <- unlist(i.l)

    return(i.v)
  })
  unlist(l)
})()

sig_metadata$mean_within_clust_cossim <- (function(){
  m <- cbind(
    cosmic_matched_sigs.mean_within_clust_cossim[sig_metadata$sig_name],
    sig_metadata$mean_within_clust_cossim
  )
  m <- t(apply(m, 1, sort, na.last=T))
  m[,1]
})()
sig_metadata$mean_within_clust_cossim <- round(sig_metadata$mean_within_clust_cossim, 2)

## Make one cos sim measure --------------------------------
sig_metadata$cos_sim_merged <- (function(){
  df <- sig_metadata[,c('sig_name','matched_ref_cos_sim','mean_within_clust_cossim')]
  #df <- subset(df, !is.na(sig_name))
  #subset(df, grepl('SBS18',sig_name))
  
  df$merged <- df$mean_within_clust_cossim
  df$merged[!is.na(df$matched_ref_cos_sim)] <- df$matched_ref_cos_sim[!is.na(df$matched_ref_cos_sim)]
  df$merged[is.na(df$merged)] <- df$mean_within_clust_cossim[is.na(df$merged)]
  
  return(df$merged)
})()

## Export ================================
if(WRITE_OUTPUT){
  post_processed_dir <- paste0(wd,'/sig_contrib/fit_lsq.post_processed/')
  dir.create(post_processed_dir, showWarnings=F)
  
  write.table(
    sig_contribs.renamed,
    gzfile(paste0(post_processed_dir,'/denovo_contribs.lsq.post_processed.txt.gz')),
    sep='\t',quote=F
  )
  
  write.table(
    sig_metadata,
    gzfile(paste0(post_processed_dir,'/sig_metadata.post_processed.txt.gz')),
    sep='\t',quote=F,row.names=F
  )
  
  saveRDS(sig_profiles.merged, paste0(post_processed_dir,'/denovo_profiles.rds'))
}



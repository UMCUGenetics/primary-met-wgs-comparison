options(stringsAsFactors=F)

## Paths ================================
path_prefix <- c(
   hpc  ='/hpc/',
   local=path.expand('~/hpc/')
)
path_prefix <- path_prefix[ min(which(dir.exists(path_prefix))) ]
path_prefix <- sub('/hpc/','/',path_prefix)
base_dir <- paste0(path_prefix,'/hpc/cuppen/projects/P0025_PCAWG_HMF/')

devtools::load_all(paste0(path_prefix,'/hpc/cuppen/projects/P0013_WGS_patterns_Diagn/CHORD/processed/scripts_main/mutSigExtractor'))

## Paths ================================
args <- commandArgs(trailingOnly=T)
wd <- args[1]
#wd <- paste0(base_dir,'/passengers/processed/sigs_denovo/extractions/06_DR104-update4_pcawgSoftFilt/output/Biliary.dbs/')

context_translations <- read.delim(paste0(base_dir,'/passengers/processed/sigs_denovo/dep/mut_context_translations.txt'))

##
message('Getting output dir from SigProfiler...')
dir.mut_type <- grep(
   '(SBS|DBS|ID)\\d+$',
   list.dirs(wd, recursive=F),
   value=T
)
if(length(dir.mut_type)>1){
   dir.mut_type <- dir.mut_type[1]
   warning('Multiple mut type dirs were found selecting the first: ', dir.mut_type)
}

getSigProfilerMatrices <- function(dir.mut_type, what){
   #dir.mut_type
   #what='signatures'
   
   ## Get paths to matrices
   dir.solutions <- paste0(dir.mut_type,'/All_Solutions/')
   if(what=='signatures'){
      paths <- list.files(dir.solutions,'*_Signatures.txt$', recursive=T, full.names=T)
   } else if(what=='activities'){
      paths <- list.files(dir.solutions,'*_Activities.txt$', recursive=T, full.names=T)
   } else {
      stop("`what` must be 'signatures' or 'activities'")
   }
   
   ## Read txt files as matrices
   matrices <- lapply(paths, function(i){
      df <- read.delim(i)
      rownames(df) <- df[,1]
      as.matrix(df[,-1,drop=F])
   })
   
   ## Use mutSigExtractor context names
   if(what=='signatures'){
      matrices <- lapply(matrices, function(i){
         rownames(i) <- 
            context_translations$name.mutSigExtractor[ 
               match(rownames(i), context_translations$name.sigProfiler) 
            ]
         return(i)
      })
   }
   
   ## Get NMF rank
   rank_name <- sapply( strsplit(basename(paths),'_'), `[[`, 2 )
   rank_name <- sub('^S','',rank_name)
   names(matrices) <- rank_name
   matrices <- matrices[ naturalsort::naturalorder(rank_name) ]
   
   ## Get mut type
   mut_type <- stringr::str_extract(basename(dir.mut_type), '^[A-Z]+')
   class(matrices) <- c(class(matrices), mut_type)
   
   return(matrices)
}

## ================================
message('Matching denovo sigs to reference sigs...')

profiles.denovo.raw <- getSigProfilerMatrices(dir.mut_type, 'signatures')

sig_metadata <- read.delim(mutSigExtractor:::SIG_METADATA_PATH)
matchDenovoToRefSigs <- function(profiles.denovo, output='profiles'){
   #profiles.denovo=profiles.denovo.raw
   #output='cossims'
   
   ## init
   output <- match.arg(output, c('cossims','profiles'))
   
   ## Get relevant ref signature profiles
   mut_type <- grep(c('SBS|DBS|ID'), class(profiles.denovo), value=T)
   profiles.ref <- list(
      SBS=mutSigExtractor::SBS_SIGNATURE_PROFILES_V3,
      DBS=mutSigExtractor::DBS_SIGNATURE_PROFILES,
      ID=mutSigExtractor::INDEL_SIGNATURE_PROFILES
   )[[mut_type]]
   
   ## Cos sim function
   calcProfileCosSim <- function(profiles.denovo, profiles.ref){
      cosSim <- function(x, y) { sum(x*y)/sqrt(sum(x^2)*sum(y^2)) }
      
      out <- apply(profiles.denovo, 2, function(i){
         #i=profiles.de_novo[,1]
         apply(profiles.ref, 2, function(j){
            cosSim(i, j)
         })
      })
      
      t(out)
   }
   
   if(output=='cossims'){
      return(
         lapply(profiles.denovo, function(i){
            calcProfileCosSim(i, profiles.ref)
         })
      )
   }
   
   lapply(profiles.denovo, function(i){
      #i=profiles.denovo[[2]]
      ## Match each denovo sigs to a ref sig
      cos_sims <- calcProfileCosSim(i, profiles.ref)
      
      max_sim <- round(apply(cos_sims,1,max),2) * 100
      max_sim_ref_sig <- colnames(cos_sims)[ max.col(cos_sims) ]
      
      sig_etiology <- sig_metadata$sig_etiology[ match(max_sim_ref_sig, sig_metadata$sig_name) ]
      
      ## Add info to denovo sig names
      sig_names <- paste0(
         max_sim_ref_sig,'.',
         sig_etiology,'.',
         max_sim,'.',
         colnames(i)
      )
      colnames(i) <- sig_names
      
      return(i)
   })
}

## ================================
message('Prepping data for plotting...')
cossims.denovo <- matchDenovoToRefSigs(profiles.denovo.raw, output='cossims')
profiles.denovo <- matchDenovoToRefSigs(profiles.denovo.raw, output='profiles')
contribs.denovo <- getSigProfilerMatrices(dir.mut_type, 'activities')

denovoSigReport <- function(profiles, contribs, cossims, rank.name){
   if(F){
      rank.name='1'
      profiles=profiles.denovo[[rank.name]]
      contribs=contribs.denovo[[rank.name]]
      cossims=cossims.denovo[[rank.name]]
   }
   
   require(ggplot2)
   
   ## Reorder matrix based on ref sig name --------------------------------
   denovo_sig_order <- naturalsort::naturalorder(colnames(profiles))
   profiles <- profiles[,denovo_sig_order,drop=F]
   contribs <- contribs[,denovo_sig_order,drop=F]
   cossims <- cossims[denovo_sig_order,,drop=F]
   
   ## Contribs --------------------------------
   pd_contribs <- reshape2::melt(contribs)
   colnames(pd_contribs) <- c('sample','signature','contrib')
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
   
   ## Calc medians
   pd_sig_contrib_med <- with(pd_contribs,{
      aggregate(contrib+1, list(signature=signature), median)
   })
   colnames(pd_sig_contrib_med)[2] <- 'contrib'
   med_pos <- nrow(contribs)*0.5
   
   p_contribs <- ggplot(pd_contribs, aes(x=index, y=contrib+1)) +
      facet_grid(signature~., scales='free_y', space='free_x') +
      geom_point(size=0.2) +
      geom_point(data=pd_sig_contrib_med, aes(y=contrib), x=med_pos, color='red', shape=95, size=10) +
      scale_y_continuous(trans='log10', name='Contribution + 1') +
      xlab('Sample rank') +
      theme_bw() +
      theme(
         #strip.text.y=element_text(angle=0),
         strip.text.y=element_blank(),
         panel.grid.major.x=element_blank(),
         panel.grid.minor.x=element_blank(),
         panel.grid.minor.y=element_blank()
         #axis.title.x=element_blank(),
         #axis.text.x=element_blank(),
         #axis.ticks.x=element_blank()
      )
   
   ## Profiles --------------------------------
   trimText <- function(v, width=25){
      #v=rownames(sig_profiles)
      #width=20
      
      v_short <- substring(v, 1, width)
      which_long <- nchar(v) > width
      v_short[which_long] <- paste0(v_short[which_long],'...')
      
      return(v_short)
   }
   
   profiles_t <- t(profiles)
   rownames(profiles_t) <- trimText(rownames(profiles_t))
   
   p_profiles <-
      mutSigExtractor::plotContexts(
         profiles_t, group=rownames(profiles_t), 
         y.axis.var.scale=T, force.group.labels=T
      ) +
      ggtitle(paste0('Rank: ', rank.name)) +
      ylab('Probability') +
      theme(
         strip.text.y=element_text(angle=0)
      )
   
   ## Cossims --------------------------------
   ## Sort and melt
   top.n <- 5
   pd_cossims <- do.call(rbind, lapply(rownames(cossims), function(i){
      v <- sort(cossims[i,], decreasing=T)[1:top.n]
      data.frame(sig=i, feature=names(v), value=v, index=1:top.n, row.names=NULL)
   }))
   pd_cossims <- as.data.frame(lapply(pd_cossims, function(i){
      if(!is.numeric(i)){ i <- factor(i, unique(i)) }
      return(i)
   }))
   
   ## Bar labels
   pd_cossims$label <- paste0(as.character(pd_cossims$feature),' (',round(pd_cossims$value,2),')')
   #pd_cossims$label <- as.character(pd_cossims$feature)

   ## Plot
   p_cossims <- ggplot(pd_cossims, aes(x=index, y=value)) +
      facet_grid(sig~.) +
      geom_bar(stat='identity', width=1, size=0.25, fill='lightsteelblue', color='black') +
      geom_text(aes(label=label), y=0.02, angle=90, hjust=0, vjust=0.5, size=2.7) +
      labs(y='Cos sim', x='Rank') +
      theme_bw() +
      theme(
         #strip.text.y=element_blank(),
         panel.grid.minor=element_blank()
      )
   
   cowplot::plot_grid(
      p_contribs, p_profiles, p_cossims,
      nrow=1, align='h', axis='tblr', rel_widths=c(0.2, 1, 0.2)
   )
}

##
message('Making plots for:')
rank_names <- names(contribs.denovo)
reports <- lapply(rank_names, function(i){
   message('> Rank: ',i)
   denovoSigReport(
      profiles=profiles.denovo[[i]], 
      contribs=contribs.denovo[[i]], 
      cossims=cossims.denovo[[i]], 
      rank.name=i
   )
})
names(reports) <- stringr::str_pad(rank_names,2,pad='0')

##
message('Exporting plots for:')
report_dir <- paste0(wd,'/reports/')
dir.create(report_dir, showWarnings=F)

for(i in names(reports)){
   #i=names(reports)[1]
   
   message('> Rank: ',i)
   
   rank <- as.integer(i)
   height <- 100*rank + 150
   
   out_path <- paste0(report_dir,'/rank_',i,'.png')
   png(out_path, width=1600, height=height, res=125)
   plot(reports[[i]])
   dev.off()
}


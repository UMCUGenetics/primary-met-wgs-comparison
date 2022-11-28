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
wd <- paste0(base_dir,'/passengers/analysis/sig_comparison/08_20221111/')
plots_dir <- paste0(wd,'/plots/'); dir.create(plots_dir, showWarnings=F)
tables_dir <- paste0(wd,'/tables/'); dir.create(tables_dir, showWarnings=F)

nmf_dir <- paste0(base_dir,'/passengers/processed/sigs_denovo/extractions/12_fixed_seeds/')

## Dependencies --------------------------------
source(paste0(base_dir,'/passengers/analysis/common_functions/common_functions.R'))
library(mutSigExtractor)
library(ggplot2)

## Contexts ########################################################################################
## Load data ================================
contexts <- readRDS(paste0(base_dir, '/passengers/processed/mut_contexts/04_fixed_smnvs/matrices/contexts_merged.rds'))
contexts <- contexts[c('snv','indel','dbs')]
names(contexts) <- c('SBS','ID','DBS')

## Mut load S-curve ================================
plotEcdfMutLoad <- function(
   contexts, sample.groups, qvalue.thres=0.05,
   legend.justification=c(0, 0.5), output='plot', pcawg.sens.loss=NULL,
   show.div.lines=F, signif.label.style='show.signif'
){
   if(F){
      sample.groups=sample_groups$cancer_stage.cancer_type_code
      qvalue.thres=0.01
      legend.justification=c(0, 0.5)
      output='plot'
      pcawg.sens.loss=NULL
      show.div.lines=T
      signif.label.style='star.signif'
   }
   
   ## Init --------------------------------
   require(ggplot2)
   require(scales)
   
   if(!(output %in% c('plot','enr'))){ 
      stop("`output` must be 'plot' or 'enr'") 
   }
   
   if(!(signif.label.style %in% c('show.all','show.signif','star.signif'))){ 
      stop("`output` must be 'show.all','show.signif','star.signif'") 
   }
    
   ## Matrix to long dataframe --------------------------------
   mut_type_load <- do.call(cbind, lapply(contexts, rowSums))
   pd <- getSampleData(sample.groups, mut_type_load)
   pd$feature <- factor(pd$feature, c('SBS','ID','DBS'))
   
   ## Subset for selected samples
   mut_type_load <- mut_type_load[rownames(mut_type_load) %in% pd$sample_id,]
   
   ##
   log10p1 <- function(x){ log10(x+1) }
   pd$value.log10p1 <- log10p1(pd$value)
   
   ## Calculate x positions --------------------------------
   pd_split <- split(pd, paste0(pd$comparison_group   ,':',pd$label_group ,':',pd$feature))
   range01 <- function(x){(x-min(x))/(max(x)-min(x))}
   pd <- do.call(rbind, lapply(pd_split, function(i){
      #i=pd_split$`Hartwig:Breast:BRCA`
      i <- i[order(i$value),]
      i$rank <- 1:nrow(i)
      i$rank_norm <- range01(i$rank)
      return(i)
   }))
   rownames(pd) <- NULL
   
   pd <- pd[order(pd$label_group , pd$feature, pd$comparison_group),]
   
   ## Group
   pd$group <- paste0(pd$feature,' ',pd$comparison_group)
   pd$group <- factor(pd$group, unique(pd$group))
   
   ## Enrichment --------------------------------
   ## Enrichment functions require wide format data
   enr_input <- getSampleData(sample.groups, mut_type_load, out.format='wide')
   
   ## Convert comparison_group to TRUE/FALSE
   enr_input$comparison_bool <- ifelse(
      enr_input$comparison_group==levels(enr_input$comparison_group)[1],
      TRUE, FALSE
   )
   
   enr_input.split <- split(enr_input, enr_input$label_group)
   feature_names <- colnames(mut_type_load)
   
   ## Main
   #subset(pd, label_group=='BRCA, Local vs Primary')
   #subset(pd, label_group=='BRCA, Local vs Primary')
   enr <- lapply(names(enr_input.split), function(i){
      #i=names(m_split)[1]
      #i='BRCA, Local vs Primary'
      #print(i)
      enr_input.ss <- enr_input.split[[i]]
      out <- univarFeatSel.default(
         x=enr_input.ss[,feature_names,drop=F],
         y=enr_input.ss$comparison_bool,
         alternative='two.sided', avg.numeric.func='median', order.by.pvalue=F
      )
      cbind(label_group=i, out)
   })
   enr <- do.call(rbind, enr)
   rownames(enr) <- NULL
   
   enr$feature <- factor(enr$feature, levels(pd$feature))
   enr$label_group <- factor(enr$label_group, levels(pd$label_group))
   NULL -> enr$is_pass_feature -> enr$is_keep_feature
   
   ## p-adjust by SV type
   enr$index <- 1:nrow(enr)
   enr.split <- split(enr, enr$feature)
   enr <- do.call(rbind, lapply(enr.split, function(i){
      i$qvalue <- p.adjust(i$pvalue, method='bonferroni')
      return(i)
   }))
   enr <- enr[order(enr$index),]
   enr$index <- NULL
   
   ## Calculate fold change
   #reverseLog10p <- function(x){ 10^x - 1 }
   enr$avg_case.log10p1 <- log10p1(enr$avg_case)
   enr$avg_ctrl.log10p1 <- log10p1(enr$avg_ctrl)
   
   enr$fold_change <- enr$avg_case / enr$avg_ctrl
   enr$fold_change.filt <- enr$fold_change
   enr$fold_change.filt[enr$qvalue>=qvalue.thres] <- 0
   
   if(output=='enr'){ return(enr) }
   
   ## Medians --------------------------------
   medians <- reshape2::melt(
      enr[c('label_group', 'feature', 'avg_case.log10p1', 'avg_ctrl.log10p1')],
      measure.vars=c('avg_case.log10p1','avg_ctrl.log10p1')
   )
   medians$comparison_group <- c(
      'avg_case.log10p1'=levels(pd$comparison_group)[1], 
      'avg_ctrl.log10p1'=levels(pd$comparison_group)[2]
   )[as.character(medians$variable)]
   
   medians$group <- paste0(medians$feature,' ',medians$comparison_group)
   medians$group <- factor(medians$group, unique(medians$group))
   
   ## Eff sizes --------------------------------
   eff_sizes <- enr[c('label_group','feature','fold_change.filt','fold_change')]
   
   ## Select numbers with significant difference
   if(signif.label.style=='show.all'){
      eff_sizes$label <- paste0(round(eff_sizes$fold_change, 1),'x')
      
   } else if(signif.label.style=='show.signif'){
      eff_sizes$label <- paste0(round(eff_sizes$fold_change.filt, 1),'x')
      eff_sizes$label[eff_sizes$fold_change.filt==0] <- ''
      eff_sizes <- eff_sizes[nchar(eff_sizes$label)!=0,]
      
   } else if(signif.label.style=='star.signif'){
      eff_sizes$label <- paste0(round(eff_sizes$fold_change, 1),'x')
      eff_sizes$star <- ifelse(eff_sizes$fold_change.filt==0,'','*')
      eff_sizes$label <- paste0(eff_sizes$star, eff_sizes$label)
      eff_sizes$star <- NULL
      
   }
   
   if(signif.label.style=='show'){
      eff_sizes$label <- paste0(round(eff_sizes$fold_change.filt, 1),'x')
      eff_sizes$label[enr$fold_change.filt==0] <- ''
      eff_sizes <- eff_sizes[nchar(eff_sizes$label)!=0,]
   } else if(signif.label.style=='all'){
      #eff_sizes$fold_change.filt <- eff_sizes$fold_change
      eff_sizes$label <- paste0(round(eff_sizes$fold_change, 1),'x')
   }
   
   ## Add new lines
   eff_sizes$n_newlines <- unlist(lapply(
      split(eff_sizes$label_group, eff_sizes$label_group),
      function(i){ seq_along(i) }
   ))
   eff_sizes$n_newlines <- eff_sizes$n_newlines - 1
   
   eff_sizes$label_newlines <- sapply(eff_sizes$n_newlines, function(i){
      if(i==0){ return('') }
      paste(rep('\n',i), collapse='')
   })
   
   eff_sizes$label <- paste0(eff_sizes$label_newlines, eff_sizes$label)
   
   ## Use Hartwig colors
   eff_sizes$group <- paste0(eff_sizes$feature,' ',levels(pd$comparison_group)[1])
   eff_sizes$group <- factor(eff_sizes$group, unique(eff_sizes$group))
   
   ## Plot params --------------------------------
   ## y-axis
   y_labels <- c(0, 10^(1:10))
   y_breaks <- y_labels+1
   y_breaks <- log10(y_breaks)
   
   powerOf10Format <- function(x) {
      #x=y_labels
      x_string <- scales::scientific(x)
      exponent <- log10(x)
      out <- paste0('10^',exponent)
      out[x==0] <- 0
      parse(text=out)
   }
   
   ## Stripes
   v_stripes <- data.frame( 
      label_group=unique(pd$label_group)
   )
   
   tmp <- rep(c(F,T), nrow(v_stripes))
   v_stripes$has_stripe <- tmp[1:nrow(v_stripes)]
   
   ## Colors
   #RColorBrewer::brewer.pal(Inf, 'Paired')
   # colors <- c(
   #    "SBS Metastatic"="#1F78B4",  
   #    "SBS Primary"="#A6CEE3",
   #    "ID Metastatic"="#33A02C",
   #    "ID Primary"="#B2DF8A",
   #    "DBS Metastatic"="#E31A1C",
   #    "DBS Primary"="#FB9A99"
   # )
   colors <- structure(
      c("#1F78B4","#A6CEE3","#33A02C","#B2DF8A","#E31A1C","#FB9A99"),
      names=levels(pd$group)
   )
   
   ## Plot --------------------------------
   p <- ggplot(pd) + 
      facet_grid(~label_group) +
      
      ## Background stripes
      geom_rect(
         data=v_stripes, 
         mapping=aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=has_stripe), 
         show.legend=F
      ) +
      scale_fill_manual(values=c('FALSE'='white','TRUE'='grey95')) +
      
      ## Axis lines
      geom_hline(yintercept=c(-Inf, Inf), size=0.5) +
      #geom_vline(xintercept=c(-Inf, Inf), size=0.5) +
      
      ## Sorted scatter plots
      geom_point(aes(x=rank_norm, y=value.log10p1, color=group), size=0.2) +
      geom_segment(
         data=medians, 
         mapping=aes(x=0.2, xend=0.8, y=value, yend=value, color=group), 
         show.legend=F
      ) +
      scale_color_manual(values=colors) +
      
      ## Eff size text
      geom_text(
         data=eff_sizes, mapping=aes(label=label, color=group, y=max(pd$value.log10p1), x=0.5), 
         vjust=1, size=2, show.legend=F
      ) +
      
      ## Axes
      scale_x_continuous(
         sec.axis=dup_axis()
      ) +
      scale_y_continuous(
         breaks=y_breaks,
         labels=powerOf10Format(y_labels),
         sec.axis=dup_axis()
      ) +
      guides(
         colour=guide_legend(order=1, override.aes=list(size=2.5))
      ) +
      labs(
         color=paste0('Numbers: change in\nmedian number of\nmutations when q<',qvalue.thres,'\n'),
         y='No. of mutations'
      ) +
      
      ## Theme
      #theme_bw() +
      theme(
         panel.background=element_rect(fill='transparent', color='transparent'),
         panel.spacing.x=unit(0,'pt'),
         panel.spacing.y=unit(0,'pt'),
         panel.grid=element_blank(),
         axis.line.y=element_line(size=0.2),
         #axis.line.x=element_line(size=0.2),
         strip.text.x=element_text(angle=90, hjust=0, vjust=0.5, size=8),
         strip.text.y=element_text(angle=0, hjust=0, vjust=0.5),
         strip.background.x=element_blank(),
         #strip.placement.x='outside',
         strip.switch.pad.grid=unit(-2,'pt'),
         
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         
         axis.title.y=element_text(size=8),
         axis.text.y=element_text(size=7),
         axis.title.y.right=element_blank(),
         axis.text.y.right=element_blank(),
         axis.ticks.y.right=element_blank(),
         
         legend.key.size=unit(10,'pt'),
         legend.spacing.y=unit(2,'pt'),
         legend.key=element_blank(),
         legend.title=element_text(size=8),
         legend.text=element_text(size=7),
         legend.justification=legend.justification,
         plot.margin=unit(c(1,1,1,1),'pt')
      )
   
   
   ## Dividing lines
   if(show.div.lines){
      v_lines <- data.frame( 
         label_group=unique(pd$label_group)
      )
      v_lines$super_group <- sub(',.*$','',v_lines$label_group)
      v_lines$has_line <- !duplicated(v_lines$super_group, fromLast=T)
      v_lines$has_line[nrow(v_lines)] <- FALSE ## Last panel doesnt need line
      
      
      p <- p + 
         geom_vline(
            data=v_lines,
            mapping=aes(xintercept=Inf, alpha=has_line),
            size=0.5
         ) +
         scale_alpha_manual(values=c('TRUE'=1, 'FALSE'=0), guide="none")
   }
   
   return(p)
}

# ## Testing
# if(F){
#    plotEcdfMutLoad(
#       contexts=contexts,
#       sample.groups=sample_groups$cancer_stage.cancer_type_code,
#       output='plot', signif.label.style='star.signif'
#    )
# }

## Stacked bar plots ================================
plotStackBarsContexts <- function(
   contexts, sample.groups,
   mut.type='SBS', sel.comparison.group=1,
   hide.x.text=T, hide.legend=F, legend.justification=c(0.5, 0.5), legend.ncol=2
){
   
   if(F){
      sample.groups=sample_groups$cancer_type.cancer_stage
      mut.type='SBS'
      sel.comparison.group=1
      #show.compact.plot=F
      hide.x.text=F
      hide.legend=F
      legend.justification=c(0, 1)
      legend.ncol=1
   }
   
   ## Constants --------------------------------
   mut_subtype_colors <- list(
      SBS=c(
         "C>A"="#06BAEB",
         "C>G"="#737373",
         "C>T"="#E12825",
         "T>A"="grey",
         "T>C"="#A1CD63",
         "T>G"="#EDC4C5"
      ),
      
      DBS=c(
         "AC>NN"="#03BCEC",
         "AT>NN"="#0165CB",
         "CC>NN"="#9FCC62",
         "CG>NN"="#006401",
         "CT>NN"="#FE9796",
         "GC>NN"="#E22823",
         "TA>NN"="#FCB066",
         "TC>NN"="#FC8100",
         "TG>NN"="#CC98FD",
         "TT>NN"="#4C0199"
      ),
      
      ID=c(
         "1bp del. C"="#EEC498",
         "1bp del. T"="#EF8633",
         "1bp ins. C"="#B0E9AF",
         "1bp ins. T"="#3F8831",
         ">1bp del. at repeats"="#BC261A",
         ">1bp ins. at repeats"="#4E9AF7",
         "Del. microhomology"="#43137D"
      )
   )
   
   mut_subtype_parsers <- list(
      SBS=function(chr){ gsub('(^\\w\\[)|(\\]\\w$)','',chr) },
      DBS=function(chr){ sub('\\w{2}$','NN',chr) },
      ID=function(chr){
         sapply(colnames(m), function(i){
            if(grepl('del.1.C', i, fixed=T)) return("1bp del. C")
            if(grepl('del.1.T', i, fixed=T)) return("1bp del. T")
            if(grepl('ins.1.C', i, fixed=T)) return("1bp ins. C")
            if(grepl('ins.1.T', i, fixed=T)) return("1bp ins. T")
            if(grepl('^del.+rep', i, fixed=F)) return(">1bp del. at repeats")
            if(grepl('^ins.+rep', i, fixed=F)) return(">1bp ins. at repeats")
            return('Del. microhomology')
         })
      }
   )
   
   ## Transform contexts --------------------------------
   m <- contexts[[mut.type]]
   
   ## Order samples by mut load
   m <- m[order(rowSums(m)),]
   
   ## Merge mut subtypes
   mut_subtypes <- mut_subtype_parsers[[mut.type]](colnames(m))
   m <- do.call(cbind, lapply(unique(mut_subtypes), function(i){
      #i='1bp del. C'
      rowSums(m[,mut_subtypes==i])
   }))
   colnames(m) <- unique(mut_subtypes)
   
   ## Rel contrib
   m <- m/rowSums(m)
   m[is.na(m)] <- 0
   
   ## Plot data --------------------------------
   pd <- getSampleData(sample.groups, m)
   
   ## Subset for cohort
   if(is.numeric(sel.comparison.group)){
      sel_comparison_group <- levels(pd$comparison_group)[ sel.comparison.group ]
   } else {
      sel_comparison_group <- sel.comparison.group
   }
   
   pd <- pd[pd$comparison_group==sel_comparison_group,]
   
   ## Force samples to order by mut load
   m <- m[rownames(m) %in% pd$sample_id,]
   pd$sample_id <- factor(pd$sample_id, rownames(m))
   
   ## Plot --------------------------------
   ggplot(pd, aes(x=sample_id, y=value, fill=feature)) +
      facet_grid(~label_group, scales='free_x', switch='y') +
      geom_bar(stat='identity', width=1) +
      
      scale_fill_manual(values=mut_subtype_colors[[mut.type]]) +
      scale_y_continuous(
         expand=c(0,0),
         breaks=c(0.25, 0.5, 0.75, 1),
         labels=function(x){ paste0(x*100,'%') }
      ) +
      guides(fill=guide_legend(ncol=legend.ncol)) +
      labs(y=pd$comparison_group[1]) +
      
      theme_bw() +
      theme(
         rect=element_rect(fill='transparent'),
         plot.background=element_rect(color='transparent'),
         panel.grid=element_blank(),
         panel.spacing.x=unit(1, 'pt'),
         panel.spacing.y=unit(0, 'pt'),
         strip.text.x=if(hide.x.text){ element_blank() } else { element_text(angle=90, hjust=0, vjust=0.5) },
         strip.background.x=element_blank(),
         strip.text.y=element_text(angle=0, hjust=0, vjust=0.5),
         strip.text.y.left=element_text(angle=0, hjust=1, vjust=0.5),
         strip.background.y=element_blank(),
         #strip.placement='outside',
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.x=element_blank(),
         axis.text.y=element_text(size=7),
         axis.title.y=element_text(size=8, angle=0, vjust=0.5, hjust=1),
         legend.title=element_blank(),
         legend.key.size=unit(10,'pt'),
         legend.text=element_text(size=7),
         legend.spacing.y=unit(2,'pt'),
         legend.position=if(hide.legend){ 'none' } else { 'right' },
         legend.justification=legend.justification,
         plot.margin=unit(c(1,1,1,1),'pt')
      )
}

## Testing
# plotStackBarsContexts(
#    contexts, sample_groups$cancer_type.cancer_stage,
#    mut.type='DBS', sel.comparison.group=1,
#    hide.legend=F, hide.x.text=F, legend.ncol=2
# )

plotStackBarsContextsAllMutTypes <- function(
   contexts, sample.groups,
   export.path=NULL, export.width=13, export.height=7
){
   
   if(F){
      sample.groups=sample_groups$cancer_stage.cancer_type_code
      export.path=NULL
      export.width=13
      export.height=7
   }
   
   ## Make plots --------------------------------
   plots <- list()
   
   plots$SBS.case <- plotStackBarsContexts(
      contexts, sample.groups,
      mut.type='SBS', sel.comparison.group=1, hide.x.text=F
   )
   
   plots$SBS.ctrl <- plotStackBarsContexts(
      contexts, sample.groups,
      mut.type='SBS', sel.comparison.group=2
   )
   
   plots$DBS.case <- plotStackBarsContexts(
      contexts, sample.groups,
      mut.type='DBS', sel.comparison.group=1
   )
   
   plots$DBS.ctrl <- plotStackBarsContexts(
      contexts, sample.groups,
      mut.type='DBS', sel.comparison.group=2
   )
   
   plots$ID.case <- plotStackBarsContexts(
      contexts, sample.groups,
      mut.type='ID', sel.comparison.group=1
   )
   
   plots$ID.ctrl <- plotStackBarsContexts(
      contexts, sample.groups,
      mut.type='ID', sel.comparison.group=2
   )
   
   ## Combine plots --------------------------------
   p <- patchwork::wrap_plots(plots, ncol=1)
   
   # pdf(paste0(plots_dir,'/context_profiles.pdf'), 13, 7)
   # plot(p)
   # dev.off()
   
   ## Export --------------------------------
   if(is.null(export.path)){
      return(p)
   }
   
   #if(verbose){ message('Exporting combined plots with W x H: ', export.width,' x ',export.height) }
   if(!grepl('[.](png|pdf)$',export.path)){
      stop("`export.path` must end with pdf or png")
   }
   
   if(grepl('[.]png$',export.path)){
      png(export.path, export.width, export.height, units='in', res=300)
      plot(p)
      dev.off()
   }
   
   if(grepl('[.]pdf$',export.path)){
      pdf(export.path, export.width, export.height)
      plot(p)
      dev.off()
   }
}

# plotStackBarsContextsAllMutTypes(
#    contexts, sample_groups$cancer_stage.cancer_type_code,
#    export.path=paste0(plots_dir,'/contexts_rel_contribs.pdf'),
#    export.width=13, 
#    export.height=7
# )

## Signatures ######################################################################################
## Load data ================================
sig_contribs <- read.delim(paste0(nmf_dir,'/sig_contrib/fit_lsq.post_processed/denovo_contribs.lsq.post_processed.txt.gz'), check.names=F)
sig_etiologies <- read.delim(paste0(nmf_dir,'/scripts/08_etiology/sig_etiologies.txt'))

## Signature contributions by etiology ================================
getSigEtiology <- function(v=NULL, show.sig.names=F, output='string'){
  #v <- colnames(sig_contribs)
  
  if(!(output %in% c('string','df'))){ stop("`output` must be 'string' or 'df'") }
  
  etiology_lookup <- structure(
    sig_etiologies$etiology,
    names=sig_etiologies$sig_name
  )
  
  if(is.null(v)){ return(etiology_lookup) }
  
  ## Match etiology to sig --------------------------------
  df <- data.frame(
    input_string=v,
    etiology=etiology_lookup[v]
  )
  
  ## Set etiology as 'Unknown' if etiology is empty
  df$is_unknown <- is.na(df$etiology) | df$etiology=='Unknown'
  df$etiology[df$is_unknown] <- 'Unknown'
  
  if(output=='df'){ return(df) }
  if(!show.sig.names){ return(df$etiology) }
  
  ## Show sig names next to etiology --------------------------------
  require(stringr)
  
  ## Only keep sigs that are in the input vector
  etiology_lookup.ss <- etiology_lookup[names(etiology_lookup) %in% v]
  etiology_lookup.ss.split <- split(names(etiology_lookup.ss), etiology_lookup.ss)
  
  ## Assign "SBS7a/7b, DBS1, ID13" notation to each etiology
  sig_lookup <- unlist(lapply(etiology_lookup.ss.split, function(i){
    #i=etiology_lookup.ss.split[['MMRd']]
    
    ## Remove denovo sigs when there are matching COSMIC sigs
    ## ...unless all sigs are denovo sigs
    if(length(i)>1 & !all(grepl('_denovo_',i))){ 
      i <- i[!grepl('_denovo_',i)] 
    }
    
    ## Split signatures by mut type
    mut_type <- str_extract(i,'[A-Z]+')
    mut_type <- factor(mut_type, unique(mut_type))
    i_split <- split(i, mut_type)
    
    ## Collapse signature names by mut type
    i_l <- lapply(names(i_split), function(j){
      #j=names(i_split)[[1]]
      sig_nums <- str_replace(i_split[[j]],'^[A-Z]+','')
      out <- paste0(j, paste(sig_nums, collapse='/'))
      out <- gsub('/_denovo_','/',out) ## Hack to simplify denovo sig names
      return(out)
    })
    
    paste(unlist(i_l), collapse=',')
  }))
  
  ## Match sig names with etiology
  df$sig_names_pre <- sig_lookup[ match(df$etiology, names(sig_lookup)) ]
  df$sig_names_pre[df$is_unknown] <- df$input_string[df$is_unknown]
  
  ## Select the correct sig type for the sig name
  df$mut_type <- str_extract(df$input_string,'[A-Z]+')
  df$mut_type <- factor(df$mut_type, unique(df$mut_type))
  
  sig_names_pre.split <- strsplit(df$sig_names_pre,',')
  df$sig_names <- sapply(1:length(sig_names_pre.split), function(i){
    #i=73
    sel_mut_type <- df$mut_type[[i]]
    strings <- sig_names_pre.split[[i]]
    mut_type <- str_extract(strings,'[A-Z]+')
    strings[mut_type==sel_mut_type]
  })
  
  ##
  df$output_string <- paste0(df$sig_names,' | ',df$etiology)
  df$index <- 1:nrow(df)
  df <- df[order(df$mut_type, df$etiology),]
  df$output_string <- factor(df$output_string, unique(df$output_string))
  df <- df[order(df$index),]
  
  return(df$output_string)
}

mergeSigsByEtiology <- function(m){
  #m=sig_contribs
  
  ## Get etiology and mut type
  sig_info <- data.frame(
    sig_name=colnames(m),
    sig_names.etiology=getSigEtiology(colnames(m), show.sig.names=T),
    mut_type=stringr::str_extract(colnames(m), '^[A-Z]+')
  )
  
  #sig_info$sig_names <- sub(' \\|.+$','',sig_info$sig_names.etiology)
  #sig_info$etiology <- sub('^.+ \\| ','',sig_info$sig_names.etiology)
  
  ## Sort groups by mut type, then etiology
  sig_info$mut_type <- factor(sig_info$mut_type, unique(sig_info$mut_type))
  sig_info <- sig_info[order(sig_info$mut_type, sig_info$sig_names.etiology),]
  #sig_info <- sig_info[order(sig_info$mut_type, sig_info$etiology, sig_info$sig_names),]
  sig_info$sig_names.etiology <- factor(sig_info$sig_names.etiology, unique(sig_info$sig_names.etiology))
  rownames(sig_info) <- NULL
  
  ## Sum sig contribs
  out <- do.call(cbind, lapply(split(sig_info, sig_info$sig_names.etiology), function(i){
    i_m <- m[,i$sig_name,drop=F]
    if(ncol(i_m)==1){ return(i_m) }
    rowSums(i_m)
  }))
  colnames(out) <- levels(sig_info$sig_names.etiology)
  
  return(out)
}

etiology_contribs <- mergeSigsByEtiology(sig_contribs)

if(WRITE_OUTPUT){
  write.table(
    etiology_contribs,
    gzfile(paste0(tables_dir,'/etiology_contribs.txt.gz')),
    sep='\t', quote=F
  )
}

## Signature contributions by etiology ================================
## yarrr::piratepal("info")
etiology_group_colors <- c(
  Age="#F6BACE", 
  APOBEC="#6B8993FF", 
  DDRD="#F6F0D4FF", 
  Environmental="#95CE8AFF", 
  ROS="#94D4D4FF", 
  Treatment="#E7695DFF",
  Unknown='lightgrey'
) 

plotSigEnrichment <- function(
   sig.contribs, sample.groups, 
   rel.contrib=F, mut.type='SBS', output='plot', 
   
   ## Thresholds
   thres.qvalue=0.05, 
   thres.log2fc=0.5, ## only for abs contrib
   thres.median.diff=0.01, ## only for rel contrib
   thres.mut.load=0,
   
   ## Plot args
   show.unknown.etiology.sigs=T, show.neg.median.diff=T, show.x.axis=T
){
   if(F){
      sig.contribs=etiology_contribs
      sample.groups=sample_groups$cancer_stage.cancer_type_code
      
      thres.qvalue=1
      thres.log2fc=0.5
      thres.median.diff=0.01
      thres.mut.load=5
      
      rel.contrib=F
      mut.type='DBS'
      output='plot'
      show.unknown.etiology.sigs=F
      show.neg.median.diff=F
      show.x.axis=T
   }
   
   ## Init --------------------------------
   require(ggplot2)
   require(gggibbous)
   require(ggpattern)
   
   if(!(output %in% c('enr','plot','list.plots', 'list.plot_data'))){
      stop("`output` must be: 'enr','plot','list.plots', 'list.plot_data'")
   }
   
   if(!(mut.type %in% c('SBS','DBS','ID'))){
      stop("`mut.type` must be: 'SBS', 'DBS', 'ID'")
   }
   
   ## Prep input --------------------------------
   m <- sig.contribs
   
   ## Subset for mut type
   mut_types <- stringr::str_extract(colnames(m), '^[A-Z]+')
   m <- m[,mut_types==mut.type]
   
   ## Rel contrib
   if(rel.contrib){
      m <- m/rowSums(m)
      m[is.na(m)] <- 0
   }

   ## Min mut load filter
   if(thres.mut.load=='auto'){
      if(rel.contrib){ 
         thres.mut.load <- 0.05 
      } else {
         thres.mut.load <- switch(
            mut.type,
            SBS=25,
            DBS=5,
            ID=10
         )
      }
   }
   if(is.numeric(thres.mut.load)){
      m[m<thres.mut.load] <- NA
   }
   
   ## Enrichment functions require wide format data
   enr_input <- getSampleData(sample.groups, m, out.format='wide')
   
   ## Convert comparison_group to TRUE/FALSE
   enr_input$comparison_bool <- ifelse(
      enr_input$comparison_group==levels(enr_input$comparison_group)[1],
      TRUE, FALSE
   )
   
   ## Enrichment --------------------------------
   enr_input.split <- split(enr_input, enr_input$label_group)
   feature_names <- colnames(m)
   
   enr <- suppressWarnings({
      lapply(names(enr_input.split), function(i){
         #i=names(m_split)[1]
         #i='BRCA'
         #print(i)
         df_ss <- enr_input.split[[i]]
         out <- univarFeatSel.default(
            x=df_ss[,feature_names,drop=F],
            y=df_ss[,'comparison_bool'],
            alternative='two.sided',
            #show.conting=T, show.sample.size=T
            order.by.pvalue=F, avg.numeric.func='median'
         )
         out$qvalue <- p.adjust(out$pvalue, method='bonferroni')
         cbind(label_group=i, out)
      })
   }) 
   enr <- do.call(rbind, enr)
   enr$feature <- factor(enr$feature, unique(enr$feature))
   enr$label_group  <- factor(enr$label_group, unique(enr$label_group))
   
   ## Mut type
   enr$mut_type <- stringr::str_extract(enr$feature, '^[A-Z]+')
   enr$mut_type <- factor(enr$mut_type, unique(enr$mut_type))
   
   ## Remove unneeded columns
   NULL -> 
      rownames(enr) -> 
      #enr$mut_type -> 
      #enr$mut_type.cancer_type_label -> 
      enr$is_pass_feature ->
      enr$is_keep_feature
   
   ## Calculate difference in median
   enr$median_diff <- enr$avg_case - enr$avg_ctrl
   enr$fold_change <- (enr$avg_case+1) / (enr$avg_ctrl+1) ## Add pseudo count to prevent x/0 or 0/x cases
   enr$fold_change.log2 <- log2(enr$fold_change)
   
   ## Misc
   if(!rel.contrib){
      enr$is_signif <- enr$qvalue<thres.qvalue & abs(enr$fold_change.log2)>=thres.log2fc
   } else {
      enr$is_signif <- enr$qvalue<thres.qvalue & abs(enr$eff_size)>=0.1 ## eff_size is cliff_delta
   }
  
   ## 
   if(!show.unknown.etiology.sigs){
      enr <- enr[!grepl('_denovo_\\d+ \\| Unknown$', enr$feature),]
   }
   
   enr$has_contrib <- enr$avg_case>0 | enr$avg_ctrl>0
   
   ## Fix NA values
   enr <- within(enr,{
      avg_case[is.na(avg_case)] <- 0
      avg_ctrl[is.na(avg_ctrl)] <- 0
      median_diff[is.na(median_diff)] <- 0
      fold_change[is.na(fold_change)] <- 0
      fold_change.log2[is.na(fold_change.log2)] <- 0
      has_contrib[is.na(has_contrib)] <- FALSE
   })
   
   if(output=='enr'){ return(enr) }
   
   ## Subset for signif
   no_enriched_features <- sum(enr$is_signif)==0
   if(!no_enriched_features){
      sel_features <- enr$feature[enr$is_signif]
   } else {
      warning('No significantly enriched signatures. Using top 3 COSMIC signatures')
      sel_features <- levels(enr$feature)[ grep('denovo',levels(enr$feature), invert=T) ][ 1:3 ]
      enr <- subset(enr, feature %in% sel_features)
   }
   enr <- subset(enr, feature %in% sel_features)
   
   ## Force factor levels
   enr$feature <- factor(enr$feature, unique(enr$feature))
   enr$label_group <- as.factor(enr$label_group)
   
   ## Main plot data --------------------------------
   ## Melt
   pd <- reshape2::melt(enr, measure.vars=c('avg_case','avg_ctrl'), value.name='median_contrib')
   
   ## Relabel comparison_group, e.g. to Metastatic and Primary
   colnames(pd)[colnames(pd)=='variable'] <- 'comparison_group'
   
   comparison_group_1 <- levels(enr_input$comparison_group)[1]
   comparison_group_2 <- levels(enr_input$comparison_group)[2]
   
   pd$comparison_group <- ifelse(
      pd$comparison_group=='avg_case',
      comparison_group_1,
      comparison_group_2
   )
   
   ## Assign right moon to e.g. Metastatic vs [Primary], where [] is the chosen factor level
   pd$is_right_moon <- pd$comparison_group==comparison_group_2
   
   ##
   if(!rel.contrib){
      
      ## Median contribution
      pd$median_contrib.trans <- log10(pd$median_contrib + 1)
      
      ## Point size
      if(mut.type=='SBS'){
         point_sizes <- list(
            break_value=c(0, 100, 500, 1000, Inf),
            point_size= c(3, 4,   5,   6),
            label=c('<100','100-500','500-1000','>1000')
         )
      } else if(mut.type=='DBS'){
         point_sizes <- list(
            break_value=c(0, 10, 20, 40, Inf),
            point_size= c(3, 4,  5,  6),
            label=c('<10','10-20','20-40','>40')
         )
      } else if(mut.type=='ID'){
         point_sizes <- list(
            break_value=c(0, 50, 100, 200, Inf),
            point_size= c(3, 4,  5,   6),
            label=c('<50','50-100','100-200','>200')
         )
      }
      
   } else {
      
      ## Median contribution
      pd$median_contrib.trans <- pd$median_contrib
      
      if(mut.type=='SBS'){
         point_sizes <- list(
            break_value=c(0, 0.01, 0.02, 0.04, Inf),
            point_size= c(3, 4,   5,   6),
            label=c('<0.01','0.02','0.04','>0.04')
         )
      } else if(mut.type=='DBS'){
         point_sizes <- list(
            break_value=c(0, 0.05, 0.1, 0.2, Inf),
            point_size= c(3, 4,  5,  6),
            label=c('<0.05','0.05-0.1','0.1-0.2','>0.2')
         )
      } else if(mut.type=='ID'){
         point_sizes <- list(
            break_value=c(0, 0.05, 0.1, 0.2, Inf),
            point_size= c(3, 4,  5,   6),
            label=c('<0.05','0.05-0.1','0.1-0.2','>0.2')
         )
      }
      
   }
   
   pd$point_size_bin <- cut(
      x=abs(pd$median_diff),
      breaks=point_sizes$break_value,
      labels=point_sizes$label,
      right=F
   )
   
   ## Add not signif and sig absent levels
   point_size_bin_levels <- levels(pd$point_size_bin)
   pd$point_size_bin <- as.character(pd$point_size_bin)
   
   pd$point_size_bin[!pd$is_signif] <- 'Diff. not significant'
   pd$point_size_bin[!pd$has_contrib] <- 'Signature absent'
   pd$point_size_bin <- factor(pd$point_size_bin, c('Signature absent','Diff. not significant',point_size_bin_levels))
   
   point_sizes_final <- structure(
      c(0.5, 2, point_sizes$point_size),
      names=levels(pd$point_size_bin)
   )
   
   ## Outer ring to indicate significant increase
   pd$enriched_in <- ifelse(pd$fold_change.log2>0, comparison_group_1, comparison_group_2)
   pd$enriched_in[ !pd$is_signif ] <- 'Neither'
   pd$outer_circle_stroke <- ifelse(pd$is_signif, 0.9, 0.1)
   
   ## Factor as integer
   pd$xpos <- as.integer(pd$label_group)
   pd$ypos <- as.integer(pd$feature)
   
   
   ## Misc plot params --------------------------------
   ## For guide lines
   max_xpos <- max(pd$xpos)
   max_ypos <- max(pd$ypos)
   
   ## For overlay circles
   pd_dedup <- subset(pd, comparison_group==comparison_group_1)
   
   ## rev(RColorBrewer::brewer.pal(Inf,'Spectral'))
   fill_colors <- c("#3288BD","#66C2A5","#ABDDA4","#E6F598","#FFFFBF","#FEE08B","#FDAE61","#F46D43","#D53E4F","#9E0142")
   
   ## Moon plot --------------------------------
   p_moon <- ggplot(pd, aes(x=label_group, y=ypos)) +
      
      ## Guide stripes and lines
      GGally::geom_stripped_cols(odd='grey96', even='grey93') +
      #geom_hline(data=h_guides, mapping=aes(yintercept=ypos+0.5), size=0.8, color='white') +
      geom_hline(yintercept=(2:max_ypos)-0.5, size=0.8, color='white') +
      
      ## Moons
      geom_moon(
         mapping=aes(ratio=0.5, right=is_right_moon, fill=median_contrib.trans, size=point_size_bin),
         stroke=0.07, key_glyph=draw_key_full_moon
      ) +
      scale_size_manual(values=point_sizes_final)
   
   ## Fill scale
   if(!rel.contrib){
      p_moon <- p_moon +
         scale_fill_gradientn(
            colors=fill_colors, na.value=fill_colors[1],
            breaks=0:5, ## 0, 10, 100, ... in log10 scale
            limits=if(mut.type=='DBS' & rel.contrib){ c(NA,NA) } else { c(1, NA) },
            labels=function(x){ 
               x <- 10^x-1 ## Convert to normal scale
               x <- 10^round(log10(x)) ## Round to nearest 10, 100, etc...
               is_min <- which.min(x)
               if(x[is_min]!=0){ x[is_min] <- paste0('<',x[is_min]) } ## If min break is not 0, set to e.g. '<10'
               return(x)
            }
         )
   } else {
      p_moon <- p_moon +
         scale_fill_gradientn(
            colors=fill_colors, na.value=fill_colors[1],
            limits=c(0,1)
         )
   }
   
   p_moon <- p_moon +
      
      ## Overlay circles to mark significance
      geom_point(
         data=pd_dedup,
         mapping=aes(size=point_size_bin, color=enriched_in),
         shape=21, fill=NA, stroke=pd_dedup$outer_circle_stroke,
      ) +
      scale_color_manual(
         breaks=c(comparison_group_1, comparison_group_2),
         values=structure(
            c('#6B469A', '#F58134', 'black'),
            names=c(comparison_group_1, comparison_group_2, 'Neither')
         )
      ) +
      
      ## Axes
      scale_x_discrete(expand=expansion(add=0.5)) +
      scale_y_reverse(
         breaks=1:max_ypos,
         labels=levels(pd$feature),
         expand=expansion(add=0.5)
      ) +
      
      ## Legend
      guides(
         colour=guide_legend(
            order=1, title.position='top', direction='vertical',
            override.aes=list(size=4, stroke=0.9)
         ),
         size=guide_legend(
            order=2, title.position='top', direction='vertical', reverse=T,
            # label.position='bottom', nrow=1, keywidth=unit(5,'pt'), 
            # label.theme=element_text(angle=90, size=7, hjust=1, vjust=0.5)
         ),
         fill=guide_colorbar(
            order=3, title.position='top', direction='horizontal',
            frame.colour='black', ticks.colour='black', 
            label.theme=element_text(size=7), barwidth=unit(70,'pt')
         )
      ) +
      
      labs(
         colour="Significant increase\nin contribution in:", 
         size=sprintf("Median diff. in\n%s contribution", mut.type), 
         fill=sprintf('Median %s contribution', mut.type)
      ) +
      
      ## Theme
      theme_bw() +
      theme(
         panel.grid=element_blank(),
         #panel.spacing.y=unit(-0.5,'pt'),
         axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8),
         axis.text.y.left=element_text(size=7),
         axis.title=element_blank(),
         legend.text=element_text(size=7),
         legend.title=element_text(size=8),
         legend.key.height=unit(10,'pt'),
         #legend.key.width=unit(5,'pt'),
         legend.spacing.y=unit(2,'pt'),
         strip.text.y=element_text(size=8),
         plot.margin=unit(c(1,1,1,1),'pt')
      )
   
   ## Plot summed median diffs --------------------------------
   ## Select relevant data
   pd_diffs <- subset(
      pd, 
      comparison_group==comparison_group_1, 
      c(label_group, mut_type, feature, median_diff, is_signif)
   )
   
   #pd_diffs <- pd_diffs[pd_diffs$is_signif,]
   
   ## Get etiologies
   pd_diffs$etiology <- sub('.+\\| ','',as.character(pd_diffs$feature))
   etiology_groups <- sig_etiologies[ match(pd_diffs$etiology, sig_etiologies$etiology), c('etiology_group_1','etiology_group_2') ]
   etiology_groups$etiology_group_1[is.na(etiology_groups$etiology_group_1)] <- 'Unknown'
   pd_diffs <- cbind(pd_diffs, etiology_groups)
   pd_diffs$feature <- NULL
   
   ## Fill etiology group 2 with etiology group 1
   pd_diffs$etiology_group_2[nchar(pd_diffs$etiology_group_2)==0] <- NA
   
   # index <- nchar(pd_diffs$etiology_group_2)==0
   # pd_diffs$etiology_group_2[index] <- pd_diffs$etiology_group_1[index]
   # rm(index)
   
   ## Aggregate diffs with the same etiology 1 and 2
   pd_diffs$median_diff.sign <- ifelse(pd_diffs$median_diff>=0,'+','-')
   pd_diffs$tmp_group <- paste(
      pd_diffs$label_group,
      pd_diffs$etiology_group_1,
      pd_diffs$etiology_group_2,
      pd_diffs$median_diff.sign,
      pd_diffs$is_signif,
      sep='::'
   )
   pd_diffs$tmp_group <- as.integer(factor(pd_diffs$tmp_group, unique(pd_diffs$tmp_group)))
   
   pd_diffs_split <- split(pd_diffs, pd_diffs$tmp_group)
   pd_diffs <- do.call(rbind, lapply(pd_diffs_split, function(i){
      #i=pd_diffs_split$
      out <- i[1,]
      out$median_diff <- sum(i$median_diff)
      #out$is_signif <- any(i$is_signif)
      return(out)
   }))
   NULL -> rownames(pd_diffs) -> pd_diffs$tmp_group
   rm(pd_diffs_split)
   
   ## Calculate fraction for supp table
   pd_diffs$tmp_group <- paste(
      pd_diffs$cancer_type_label,
      pd_diffs$median_diff.sign,
      sep='::'
   )
   pd_diffs$tmp_group <- as.integer(factor(pd_diffs$tmp_group, unique(pd_diffs$tmp_group)))
   
   pd_diffs_split <- split(pd_diffs, pd_diffs$tmp_group)
   pd_diffs <- do.call(rbind, lapply(split(pd_diffs, pd_diffs$tmp_group), function(i){
      i$median_diff.total <- sum(i$median_diff)
      i$median_diff.frac <- i$median_diff/i$median_diff.total
      return(i)
   }))
   rm(pd_diffs_split)
   
   ## Force stack order
   pd_diffs_split <- split(pd_diffs, pd_diffs$tmp_group)
   pd_diffs <- do.call(rbind, lapply(pd_diffs_split, function(i){
      #i=pd_diffs_split$`5`
      i <- i[order(-i$is_signif, i$etiology, decreasing=T),]
      i$stack_index <- 1:nrow(i)
      return(i)
   }))
   rm(pd_diffs_split)
   NULL -> rownames(pd_diffs) -> pd_diffs$tmp_group
   
   ## Remove primary cancer enrichment
   if(!show.neg.median.diff){
      pd_diffs <- pd_diffs[pd_diffs$median_diff>=0,]
   }
   
   ## Plot
   if(nrow(pd_diffs)>0){
      p_diff <- ggplot(pd_diffs, aes(x=label_group, y=median_diff)) +
         
         geom_col_pattern(
            aes(group=stack_index, fill=etiology_group_1, pattern_fill=etiology_group_2, size=is_signif),
            position='stack', width=0.5,
            color='black', pattern_color=NA, pattern_angle=45, pattern_density=0.5, pattern_spacing=0.1,
            #size=ifelse(pd_diffs$is_signif, 0.3, 0)
         ) +
         scale_size_manual(values=c('TRUE'=0.3,'FALSE'=0), guide='none') +
         
         { if(show.neg.median.diff) geom_hline(yintercept=0, size=0.25) } +
         
         scale_fill_manual(values=etiology_group_colors) +
         scale_pattern_fill_manual(values=etiology_group_colors, guide='none') +
         scale_x_discrete(drop=FALSE, expand=expansion(add=0.5)) +
         
         guides(
            fill=guide_legend(override.aes=list(size=1))
         ) +
         labs(
            #x='Cancer type', 
            y=sprintf('Median diff.'), 
            fill='Etiology group'
         ) +
         
         theme_bw() +
         theme(
            panel.grid=element_blank(),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8),
            axis.text.y=element_text(size=7),
            axis.title=element_text(size=8),
            legend.title=element_text(size=8),
            legend.text=element_text(size=7),
            legend.key.size=unit(10,'pt'),
            legend.spacing.y=unit(2,'pt'),
            plot.margin=unit(c(3,1,1,1),'pt')
         )
   } else{
      pd_diffs <- data.frame(
         label_group=factor(levels(pd$label_group),levels(pd$label_group)),
         median_diff=0
      )
      
      p_diff<- ggplot(pd_diffs, aes(x=label_group, y=median_diff)) +
         geom_col() +
         scale_x_discrete(drop=FALSE, expand=expansion(add=0.5)) +
         labs(
            #x='Cancer type', 
            y=sprintf('Median diff.'), 
            fill='Etiology group'
         ) +
         theme_bw() +
         theme(
            panel.grid=element_blank(),
            axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=8),
            axis.text.y=element_text(size=7),
            axis.title=element_text(size=8),
            legend.title=element_text(size=8),
            legend.text=element_text(size=7),
            legend.key.size=unit(10,'pt'),
            legend.spacing.y=unit(2,'pt'),
            plot.margin=unit(c(3,1,1,1),'pt')
         )
   }
   
   
   ## Plot, no. of cancer types with signif metastatic enrichment --------------------------------
   if(!no_enriched_features){
      pd_signif <- subset(
         pd, 
         comparison_group==comparison_group_1 & is_signif, 
         c(label_group, mut_type, feature, enriched_in)
      )
      
      if(!no_enriched_features){
         pd_signif <- with(
            pd_signif, 
            aggregate(enriched_in==comparison_group_1, list(feature=feature, mut_type=mut_type), sum)
         )
         colnames(pd_signif)[ncol(pd_signif)] <- 'count'
      } else {
         pd_signif$count <- numeric()
      }
      
      
      ## Get etiology groups
      pd_signif$etiology <- sub('.+\\| ','',as.character(pd_signif$feature))
      etiology_groups <- sig_etiologies[ match(pd_signif$etiology, sig_etiologies$etiology), c('etiology_group_1','etiology_group_2') ]
      etiology_groups$etiology_group_1[is.na(etiology_groups$etiology_group_1)] <- 'Unknown'
      pd_signif <- cbind(pd_signif,etiology_groups)
      rownames(pd_signif) <- NULL
      
      ## Fill etiology group 2 with etiology group 1
      pd_signif$etiology_group_2[nchar(pd_signif$etiology_group_2)==0] <- NA
      
      ##
      max_signif_count <- max(pd_signif$count)
      p_signif <- ggplot(pd_signif, aes(x=count, y=feature)) +
         geom_col_pattern(
            aes(fill=etiology_group_1, pattern_fill=etiology_group_2),
            color='black', size=0.3, width=0.7,
            pattern_color=NA, pattern_angle=45, pattern_spacing=0.2, pattern_density=0.5
         ) +
         geom_text(aes(x=max_signif_count*0.02, label=count), hjust=0, size=2.7) +
         
         scale_fill_manual(values=etiology_group_colors) +
         scale_pattern_fill_manual(values=etiology_group_colors, guide='none') +
         scale_y_discrete(limits=rev, expand=expansion(add=0.5)) +
         
         guides(
            fill=guide_legend(override.aes=list(size=1))
         ) +
         labs(
            x=paste0('Cancer types with\nenrichment in\n', comparison_group_1,' vs ',comparison_group_2),
            y='Median\ndifference', 
            fill='Etiology group'
         ) +
         
         theme_bw() +
         theme(
            panel.grid=element_blank(),
            axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=7),
            axis.title=element_text(size=8),
            legend.title=element_text(size=8),
            legend.text=element_text(size=7),
            legend.key.size=unit(10,'pt'),
            legend.spacing.y=unit(2,'pt'),
            plot.margin=unit(c(1,1,1,1),'pt')
         )
   } else {
      
      pd_signif <- data.frame(
         feature=factor(sel_features, levels(pd$feature)),
         count=0
      )
      p_signif <- ggplot(pd_signif, aes(x=count, y=feature)) +
         geom_col() +
         scale_y_discrete(limits=rev, expand=expansion(add=0.5)) +
         labs(
            x=paste0('Cancer types with\nenrichment in\n', comparison_group_1,' vs ',comparison_group_2),
            y='Median\ndifference', 
            fill='Etiology group'
         ) +
         theme_bw() +
         theme(
            panel.grid=element_blank(),
            axis.text.x=element_text(size=8),
            axis.text.y=element_text(size=7),
            axis.title=element_text(size=8),
            legend.title=element_text(size=8),
            legend.text=element_text(size=7),
            legend.key.size=unit(10,'pt'),
            legend.spacing.y=unit(2,'pt'),
            plot.margin=unit(c(1,1,1,1),'pt')
         )
   }
   
   
   
   ## Plot data--------------------------------
   if(output=='list.plot_data'){
      return(list(
         main=pd,
         median_diff=pd_diffs,
         signif_count=pd_signif
      ))
   }
   
   ## Combine plots --------------------------------
   require(patchwork)
   
   plots <- list()
   
   plots$diff <- p_diff + 
      scale_y_continuous(expand=c(0,0)) +
      theme(
         axis.title.x=element_blank(), 
         axis.text.x=element_blank(), 
         axis.ticks.x=element_blank(),
         axis.title.y=element_text(angle=0, hjust=1, vjust=0.5),
         panel.border=element_blank(),
         axis.line.y=element_line(size=0.25)
      )
   
   plots$spacer <- plot_spacer() + 
      theme(
         plot.margin=unit(c(1,1,1,1),'pt')
      )
   
   plots$moon <- p_moon +
      theme(
         panel.border=element_blank(),
         axis.line=element_line(size=0.25)
      )
   
   #if(is.null(plot.signif.max.y)){ plot.signif.max.y <- max_signif_count }
   #if(plot.signif.max.y<max_signif_count){ stop("`plot.signif.max.y` must be >", max_signif_count) }
   plots$signif <- p_signif + 
      guides(fill='none') +
      #scale_x_continuous(expand=c(0,0), limits=c(0, max(max_signif_count, plot.signif.max.y))) +
      theme(
         axis.title.y=element_blank(), 
         axis.text.y=element_blank(), 
         axis.ticks.y=element_blank(),
         panel.border=element_blank(),
         axis.line.x=element_line(size=0.25)
      )
   
   if(!show.x.axis){
      plots$moon <- plots$moon + 
         theme(
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank()
         )
      
      plots$signif <- plots$signif +
         theme(
            axis.text.x=element_blank(),
            #axis.ticks.x=element_blank(),
            axis.title.x=element_blank()
         )
   }
   
   if(output=='list.plots'){ return(plots) }
   
   wrap_plots(
      plots,
      guides='collect',
      heights=c(1,2), widths=c(5,1)
   )
}

# ## Test
# if(F){
#    plotSigEnrichment(
#       etiology_contribs, sample_groups$cancer_stage.cancer_type_code,
#       rel.contrib=F, mut.type='DBS', output='plot', thres.mut.load='auto',
#       show.neg.median.diff=F, show.unknown.etiology.sigs=F, show.x.axis=T
#    )
   
#    enr <- plotSigEnrichment(
#       etiology_contribs, sample_groups$progression_status_code.cancer_type_code,
#       rel.contrib=F, mut.type='ID', output='enr', thres.mut.load='auto',
#       show.neg.median.diff=F, show.unknown.etiology.sigs=F, show.x.axis=T
#    )
#    enr <- enr[order(enr$feature),]
#    subset(enr, !(avg_case==0 & avg_ctrl==0))
# }

## Main figure ######################################################################################
plotSmnvEnrichment <- function(
   sample.groups, 
   mut.contexts=contexts,
   etiology.contribs=etiology_contribs,
   sig.rel.contrib=F,
   sig.thres.mut.load='auto',
   sigs.show.unknown.etiology.sigs=F,
   sigs.show.neg.median.diff=F,
   patchwork.widths=NULL,
   patchwork.heights=NULL,
   export.width=NULL,
   export.height=NULL,
   export.path=NULL,
   verbose=T,
   
   ecdf.signif.label.style='show.signif',
   sig.thres.qvalue=0.05,
   
   x.label.height=7,
   ecdf.plot.height.scale=1,
   context.plot.height.scale=1,
   diff.plot.height.scale=1
){
   
   if(F){
      mut.contexts=contexts
      etiology.contribs=etiology_contribs
      sample.groups=sample_groups$cancer_stage.cancer_type_code
      sig.rel.contrib=F
      sig.thres.mut.load='auto'
      sigs.show.unknown.etiology.sigs=F
      sigs.show.neg.median.diff=F
      export.path=NULL
      patchwork.widths=NULL
      patchwork.heights=NULL
      export.width=NULL
      export.height=NULL
      verbose=T
      
      ecdf.signif.label.style='star.signif'
      sig.thres.qvalue=0.05
      
      x.label.height=NULL
      ecdf.plot.height.scale=1
      context.plot.height.scale=1
      diff.plot.height.scale=1
      
      
   }
   
   require(patchwork)

   ## Make a list of all necessary plots --------------------------------
   if(verbose){ message('Generating individual plots...') }
   plots <- c(
      list(
         ecdf.main = plotEcdfMutLoad(contexts=mut.contexts, sample.groups, signif.label.style=ecdf.signif.label.style),
         ecdf.spacer = plot_spacer() + theme(plot.margin=unit(c(1,1,1,1),'pt'))
      ),
      
      list(
         contexts.hartwig.main = plotStackBarsContexts(contexts=mut.contexts, sample.groups, mut.type='SBS', sel.comparison.group=1),
         contexts.hartwig.spacer = plot_spacer() + theme(plot.margin=unit(c(1,1,1,1),'pt')),
         contexts.pcawg.main = plotStackBarsContexts(contexts=mut.contexts, sample.groups, mut.type='SBS', sel.comparison.group=2),
         contexts.pcawg.spacer = plot_spacer() + theme(plot.margin=unit(c(1,1,1,1),'pt'))
      ),
      
      unlist(recursive=F, x=list(
         SBS=plotSigEnrichment(
            sig.contribs=etiology.contribs, 
            sample.groups=sample.groups,
            mut.type='SBS', output='list.plots', show.x.axis=F,
            show.neg.median.diff=sigs.show.neg.median.diff,
            show.unknown.etiology.sigs=sigs.show.unknown.etiology.sigs,
            rel.contrib=sig.rel.contrib,
            thres.qvalue=sig.thres.qvalue,
            thres.mut.load=sig.thres.mut.load
         ),
         
         DBS=plotSigEnrichment(
            sig.contribs=etiology.contribs, 
            sample.groups=sample.groups,
            mut.type='DBS', output='list.plots', show.x.axis=F,
            show.neg.median.diff=sigs.show.neg.median.diff, 
            show.unknown.etiology.sigs=sigs.show.unknown.etiology.sigs, 
            rel.contrib=sig.rel.contrib,
            thres.qvalue=sig.thres.qvalue,
            thres.mut.load=sig.thres.mut.load
         ),
         
         ID=plotSigEnrichment(
            sig.contribs=etiology.contribs, 
            sample.groups=sample.groups,
            mut.type='ID', output='list.plots', show.x.axis=T,
            show.neg.median.diff=sigs.show.neg.median.diff, 
            show.unknown.etiology.sigs=sigs.show.unknown.etiology.sigs, 
            rel.contrib=sig.rel.contrib,
            thres.qvalue=sig.thres.qvalue,
            thres.mut.load=sig.thres.mut.load
         )
      ))
   )
   
   # plotSigEnrichment(
   #    sig.contribs=etiology.contribs, 
   #    sample.groups=sample.groups,
   #    mut.type='SBS', show.x.axis=F,
   #    show.neg.median.diff=F,
   #    show.unknown.etiology.sigs=F,
   #    rel.contrib=sig.rel.contrib
   # )
   
   ## Fix max x-axis value for right plot (no. of signif cancer types) --------------------------------
   if(verbose){ message('Calculating signif plot x-axis limits...') }
   
   getAxisLimits <- function(ggobject, axis='x', side='both'){
      
      if(!(axis %in% c('x','y'))){ stop("`axis` must be 'x' or 'y'") }
      range_name <- if(axis=='x'){ 'x.range' } else { 'y.range' }
      
      axis_range <- ggplot_build(ggobject)$layout$panel_params[[1]][[range_name]]
      
      if(!(side %in% c('both','min','max'))){ stop("`side` must be 'both', 'min', or 'max'") }
      if(side=='both'){ return(axis_range) }
      if(side=='min'){ return(axis_range[1]) }
      if(side=='max'){ return(axis_range[2]) }
   }
   
   plots.signif.max_x <- max(
      getAxisLimits(plots$SBS.signif, side='max'),
      getAxisLimits(plots$DBS.signif, side='max'),
      getAxisLimits(plots$ID.signif, side='max')
   )
   
   suppressMessages({
      scale_x_template <- scale_x_continuous(
         expand=c(0,0), 
         limits=c(0, plots.signif.max_x),
         breaks=scales::pretty_breaks()
      )
      plots$SBS.signif <- plots$SBS.signif + scale_x_template
      plots$DBS.signif <- plots$DBS.signif + scale_x_template
      plots$ID.signif <- plots$ID.signif + scale_x_template
   })
   
   
   ## Auto-calculate plot dimensions --------------------------------
   if(verbose & (is.null(patchwork.widths) | is.null(patchwork.heights))){
      message('Calculating plot heights...') 
   }
   
   #' ECDF:
   #' E- E0
   #' 
   #' Contexts:
   #' C- C0
   #' C- C0
   #' 
   #' Signatures:
   #' SD S0
   #' SM SS
   #' DD D0
   #' DM DS
   #' ID I0
   #' IM IS 
   #'
   #' Total rows = 9
   
   ## Auto width
   if(is.null(patchwork.widths)){
      
      #patchwork.widths <- c(12,1)
      
      x_range <- getAxisLimits(plots$SBS.moon, axis='x', side='both')
      x_range <- x_range[2] - x_range[1]
      
      patchwork.widths <- c(
         x_range,
         x_range/10
      )
   }
   
   ## Auto heights
   if(is.null(patchwork.heights)){
      
      heights.moon_plots <- lapply(
         list(SBS=plots$SBS.moon, DBS=plots$DBS.moon, ID=plots$ID.moon),
         function(i){
            #i=plots$SBS.moon
            y_range <- getAxisLimits(i, axis='y')
            abs(y_range[2] - y_range[1])
         }
      )
      
      heights.ecdf_plot <- 8 * ecdf.plot.height.scale
      heights.context_plots <- 3 * context.plot.height.scale ## single context plot height
      heights.diff_plots <- 4 * diff.plot.height.scale
      
      patchwork.heights <- c(
         heights.ecdf_plot, ## ECDF
         
         heights.context_plots, ## Contexts, Hartwig
         heights.context_plots, ## Contexts, PCAWG
         
         heights.diff_plots, ## SBS, diff
         heights.moon_plots$SBS, ## SBS, main
         heights.diff_plots, ## DBS, diff
         heights.moon_plots$DBS, ## DBS, main
         heights.diff_plots, ## ID, diff
         heights.moon_plots$ID  ## ID, main
      )
   }
      
   ## Combine --------------------------------
   if(verbose){ message('Combining plots...') }
   
   p <- patchwork::wrap_plots(
      plots,
      guides='collect', 
      ncol=2,
      widths=patchwork.widths,
      heights=patchwork.heights
   )
   
   ## Auto-calculate export dimensions --------------------------------
   # if(is.null(x.label.height)){
   #    max_x_label_nchar <- max(nchar(unique(as.character(plots$ecdf.main$data$label_group))))
   #    char_unit_height <- 7/17
   #    x_label_height <- max_x_label_nchar * char_unit_height
   # } else {
   #    x_label_height <- x.label.height
   # }
   
   if(is.null(export.width)){
      export.width <- 0.35 * (sum(patchwork.widths) + 7 + 5) ## patchwork widths + axis label width + legend widths
   }
   
   if(is.null(export.height)){
      export.height <- 0.2 * (sum(patchwork.heights) + 2*x.label.height) ## patchwork heights + (2 x axis label height)
   }
   
   # if(F){
   #    pdf(paste0(plots_dir,'/smnv_main.pdf'), export.width, export.height)
   #    plot(p)
   #    dev.off()
   # }

   
   ## Export --------------------------------
   if(is.null(export.path)){
      return(p)
   }
   
   if(verbose){ message('Exporting combined plots with W x H: ', export.width,' x ',export.height) }
   if(!grepl('[.](png|pdf)$',export.path)){
      stop("`export.path` must end with pdf or png")
   }
   
   if(grepl('[.]png$',export.path)){
      png(export.path, export.width, export.height, units='in', res=300)
      plot(p)
      dev.off()
   }
   
   if(grepl('[.]pdf$',export.path)){
      pdf(export.path, export.width, export.height)
      plot(p)
      dev.off()
   }
}

if(WRITE_OUTPUT){
   ## Absolute contribution --------------------------------
   plotSmnvEnrichment(
      sample_groups$cancer_stage.cancer_type,
      export.path=paste0(plots_dir,'/smnv--abs--cancer_stage.cancer_type.pdf'),
      ecdf.signif.label.style='star.signif',
      x.label.height=14
   )
   
   plotSmnvEnrichment(
      sample_groups$cancer_stage.cancer_type_code,
      export.path=paste0(plots_dir,'/smnv--abs--cancer_stage.cancer_type_code.pdf'),
      ecdf.signif.label.style='star.signif'
   )
   
   plotSmnvEnrichment(
      sample_groups$cancer_stage.cancer_type_code_subtype_subset,
      export.path=paste0(plots_dir,'/smnv--abs--cancer_stage.cancer_type_code_subtype_subset.pdf'),
      ecdf.signif.label.style='star.signif'
   )
   
   plotSmnvEnrichment(
      sample_groups$progression_status_code.cancer_type_code,
      export.path=paste0(plots_dir,'/smnv--abs--progression_status_code.cancer_type_code.pdf'),
      ecdf.signif.label.style='star.signif',
      patchwork.widths=c(3,1),
      export.height=12, export.width=8
   )
   
   plotSmnvEnrichment(
      sample_groups$metastatic_location.cancer_type_code_subtype_subset,
      export.path=paste0(plots_dir,'/smnv--abs--metastatic_location.cancer_type_code_subtype_subset.pdf'),
      ecdf.signif.label.style='star.signif',
      export.height=14
   )
   
   plotSmnvEnrichment(
      sample_groups$metastatic_location.cancer_type_code,
      export.path=paste0(plots_dir,'/smnv--abs--metastatic_location.cancer_type_code.pdf'),
      ecdf.signif.label.style='star.signif'
   )
   
   ## Relative contribution --------------------------------
   plotSmnvEnrichment(
      sample_groups$cancer_stage.cancer_type, sig.rel.contrib=T,
      export.path=paste0(plots_dir,'/smnv--rel--cancer_stage.cancer_type.pdf'),
      ecdf.signif.label.style='star.signif',
      x.label.height=14
   )
   
   plotSmnvEnrichment(
      sample_groups$cancer_stage.cancer_type_code, sig.rel.contrib=T,
      export.path=paste0(plots_dir,'/smnv--rel--cancer_stage.cancer_type_code.pdf'),
      ecdf.signif.label.style='star.signif'
   )
   
   plotSmnvEnrichment(
      sample_groups$cancer_stage.cancer_type_code_subtype_subset, sig.rel.contrib=T,
      export.path=paste0(plots_dir,'/smnv--rel--cancer_stage.cancer_type_code_subtype_subset.pdf'),
      ecdf.signif.label.style='star.signif'
   )
   
   plotSmnvEnrichment(
      sample_groups$progression_status_code.cancer_type_code, sig.rel.contrib=T,
      export.path=paste0(plots_dir,'/smnv--rel--progression_status_code.cancer_type_code.pdf'),
      ecdf.signif.label.style='star.signif',
      patchwork.widths=c(3,1),
      export.height=12, export.width=8
   )
   
   plotSmnvEnrichment(
      sample_groups$metastatic_location.cancer_type_code_subtype_subset, sig.rel.contrib=T,
      export.path=paste0(plots_dir,'/smnv--rel--metastatic_location.cancer_type_code_subtype_subset.pdf'),
      ecdf.signif.label.style='star.signif'
   )
   
   plotSmnvEnrichment(
      sample_groups$metastatic_location.cancer_type_code, sig.rel.contrib=T,
      export.path=paste0(plots_dir,'/smnv--rel--metastatic_location.cancer_type_code.pdf'),
      ecdf.signif.label.style='star.signif'
   )
   
   ## All signatures --------------------------------
   plotSmnvEnrichment(
      sample_groups$cancer_stage.cancer_type_code,
      export.path=paste0(plots_dir,'/smnv--abs_full--cancer_stage.cancer_type_code--all_sigs.pdf'),
      ecdf.signif.label.style='star.signif',
      sigs.show.unknown.etiology.sigs=T,
      sigs.show.neg.median.diff=T
   )
}

## Denovo sigs analysis ######################################################################################
## Load data ================================
sig_metadata <- read.delim(paste0(nmf_dir,'/sig_contrib/fit_lsq.post_processed/sig_metadata.post_processed.txt.gz'), check.names=F)
denovo_profiles <- readRDS(paste0(nmf_dir,'/sig_contrib/fit_lsq.post_processed/denovo_profiles.rds'))

## Group denovo sig profiles ================================
getDenovoSigProfileBySigName <- function(
   sig.name, sig.metadata=sig_metadata, denovo.profiles=denovo_profiles,
   show.warnings=T
){
   if(F){
      sig.name='SBS25'
      sig.metadata=sig_metadata
      denovo.profiles=denovo_profiles
      show.warnings=T
   }
   
   main <- function(i){
      #i=sig.name
      mut_type <- stringr::str_extract(i, '[A-Z]{2,3}')
      target_denovo_sigs <- subset(sig.metadata, sig_name==i, denovo_name_uniq, drop=T)
      
      missing_denovo_sigs <- target_denovo_sigs[!(target_denovo_sigs %in% colnames(denovo.profiles[[mut_type]]))]
      if(show.warnings & length(missing_denovo_sigs)>0){
         warning(
            'Some denovo sigs matched to ',i,' are absent in `denovo.profiles`:\n',
            paste(missing_denovo_sigs, collapse=', ')
         )
      }
      
      target_denovo_sigs <- target_denovo_sigs[!(target_denovo_sigs %in% missing_denovo_sigs)]
      denovo.profiles[[mut_type]][,target_denovo_sigs,drop=F]
   }
   
   if(length(sig.name)==1){ return(main(sig.name)) }
   
   out <- lapply(sig.name, main)
   names(out) <- sig.name
   return(out)
}

cosSimToCosmic <- function(m, target.cosmic.sigs=NULL, output='plot'){
   if(F){
      m=getDenovoSigProfileBySigName('DBS_denovo_1')
      m=data_ss[[i]]
      output='plot'
      target.cosmic.sigs=c('SBS1')
   }
   
   ## Checks --------------------------------
   if(!(output %in% c('raw','plot'))){
      stop("`output` must be 'raw' or 'plot'")
   }
   
   if(!(is.data.frame(m) | is.matrix(m))){
      stop("`m` must be a data frame or matrix")
   }
   
   ## Parse mut type of input matrix --------------------------------
   mut_type <- NULL
   if(any(grepl('\\[', rownames(m)))){
      mut_type <- 'SBS'
   } else if(any(grepl('[A-Z]{2}>[A-Z]{2}', rownames(m)))){
      mut_type <- 'DBS'
   } else if(any(grepl('del|ins', rownames(m)))){
      mut_type <- 'ID'
   }
   if(is.null(mut_type)){ stop('Input contains invalid mut types') }
   
   ## Get relevant COSMIC sig profile matrix --------------------------------
   cosmic_profiles <- switch (
      mut_type,
      SBS=mutSigExtractor::SBS_SIGNATURE_PROFILES_V3,
      DBS=mutSigExtractor::DBS_SIGNATURE_PROFILES,
      ID=mutSigExtractor::INDEL_SIGNATURE_PROFILES
   )
   
   if(!is.null(target.cosmic.sigs)){
      cosmic_profiles <- cosmic_profiles[,colnames(cosmic_profiles) %in% target.cosmic.sigs, drop=F]
   }
   
   
   ## Calc cos sim --------------------------------
   cos_sims <- mutSigExtractor:::cosSim.matrix(
      t(m), 
      t(cosmic_profiles), 
      all.vs.all=T
   )
   
   if(is.vector(cos_sims)){
      cos_sims <- t(cos_sims)
      rownames(cos_sims) <- colnames(cosmic_profiles)
   }
   
   if(output=='raw'){ return(cos_sims) }
   
   ## Plot --------------------------------
   require(ggplot2)
   pd <- reshape2::melt(as.matrix(cos_sims))
   colnames(pd) <- c('ref_sig','input_sig','cos_sim')
   ggplot(pd, aes(x=ref_sig, y=input_sig)) +
      coord_cartesian(expand=0) +
      geom_tile(aes(fill=cos_sim), color='black', size=0.2) +
      geom_text(aes(label=round(cos_sim, 2))) +
      scale_fill_distiller(palette='YlGnBu', limits=c(0,1)) +
      theme_bw() +
      theme(
         panel.grid=element_blank()
      )
}

## Plot denovo sigs grouped by sig name ================================
## Plot data --------------------------------
denovo_profiles.by_sig_name <- (function(){
   l <- lapply(colnames(sig_contribs), function(i){
      #print(i)
      #i='SBS41'
      getDenovoSigProfileBySigName(i)
   })
   names(l) <- colnames(sig_contribs)
   return(l)
})()

## Main --------------------------------
plotDenovoProfilesBySigName <- function(
   data=denovo_profiles.by_sig_name, 
   mut.type='SBS', output='plots.list',
   sig.etiologies=sig_etiologies,
   export.path=NULL, export.width=NULL, export.height=NULL, export.height.scale=0.8,
   verbose=T
){
   if(F){
      data=denovo_profiles.by_sig_name
      mut.type='SBS'
      output='plots.combined'
      sig.etiologies=sig_etiologies
      export.path=NULL
      export.width=NULL
      export.height=NULL
      export.height.scale=0.8
      verbose=T
   }
   
   ## Init --------------------------------
   require(ggplot2)
   require(patchwork)
   
   if(!(mut.type %in% c('SBS','DBS','ID'))){
      stop("`mut.type` must be 'SBS', 'DBS' or 'ID'")
   }
   
   if(!(output %in% c('plots.list','plots.combined'))){
      stop("`output` must be 'plots.list' or 'plots.combined'")
   }
   
   ## Subset by mut type --------------------------------
   data_ss <- data[grepl(paste0('^',mut.type), names(data))]
   data_ss <- data_ss[sapply(data_ss, ncol)>0 ]
   
   ## Plot sig profiles --------------------------------
   if(verbose){ message('Generating ggplots...') }
   plots <- lapply(names(data_ss), function(i){
      #print(i)
      #i='SBS1'
      m <- data_ss[[i]]
      
      ## Add cos sim to COSMIC for non-denovo clusters
      if(!grepl('denovo',i)){
         cos_sims <- cosSimToCosmic(m, target.cosmic.sigs=i, output='raw')
         cos_sims <- cos_sims[1,]
         cos_sims <- round(cos_sims, 2)
         colnames(m) <- paste0(colnames(m), ' | ', cos_sims)
      }
      
      ## Add etiology
      etiology <- sig.etiologies$etiology[ match(i, sig.etiologies$sig_name) ]
      if(is.na(etiology)){ etiology <- 'Unknown' }
      
      plotContexts(t(m), group=colnames(m), force.group.labels=T) +
         #ggtitle(paste0('Assigned to: ',i)) +
         labs(
            title=paste0(i,' | Etiology: ', etiology),
            y='Probability'
         ) +
         theme(
            panel.grid.major.y=element_blank(),
            panel.spacing.y=unit(3,'pt'),
            strip.text.y=element_text(angle=0, hjust=0),
            #axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            plot.title=element_text(face='bold')
         )
   })
   
   if(is.null(export.path) & output=='plots.list'){
      return(plots)
   }
   
   ## Calculate plot dimensions --------------------------------
   if(verbose){ message('Combining plots...') }
   patchwork_heights <- sapply(data_ss, ncol)
   plots_combined <- patchwork::wrap_plots(
      plots,
      guides='collect', 
      ncol=1,
      heights=patchwork_heights
   )
   
   if(is.null(export.path) & output=='plots.combined'){ return(plots_combined) }
   
   if(is.null(export.width)){ export.width <- 14 }
   
   if(is.null(export.height)){ 
      export.height <- sum(patchwork_heights) * export.height.scale
      if(export.height>200){ 
         message('`export.height` (',export.height,') > 200 inches. Using 200 inches as export height') 
         export.height <- 200
      }
   }
   
   # pdf(paste0(plots_dir,'/denovo_sig_profiles.pdf'), export.width, 200)
   # plot(plots_combined)
   # dev.off()
   
   ## Export --------------------------------
   if(verbose){ message('Exporting combined plots with W x H: ', export.width,' x ',export.height) }
   if(!grepl('[.](png|pdf)$',export.path)){
      stop("`export.path` must end with pdf or png")
   }
   
   if(grepl('[.]png$',export.path)){
      png(export.path, export.width, export.height, units='in', res=300)
      plot(plots_combined)
      dev.off()
   }
   
   if(grepl('[.]pdf$',export.path)){
      pdf(export.path, export.width, export.height)
      plot(plots_combined)
      dev.off()
   }
}

## Export --------------------------------
if(WRITE_OUTPUT){
   plotDenovoProfilesBySigName(mut.type='SBS', export.path=paste0(plots_dir,'/sig_profiles--SBS.pdf'))
   plotDenovoProfilesBySigName(mut.type='DBS', export.path=paste0(plots_dir,'/sig_profiles--DBS.pdf'))
   plotDenovoProfilesBySigName(mut.type='ID', export.path=paste0(plots_dir,'/sig_profiles--ID.pdf'))
}

## Cos sim of denovo sigs to cosmic sigs ================================
## --------------------------------
cosSimToCosmicAllDenovo <- function(
   data=denovo_profiles.by_sig_name,
   export.path=NULL
){
   mut_types <- c('SBS','DBS','ID')
   
   ## Annotated dataframe of cos sims. Rows: denovo sigs, columns: cosmic sigs --------------------------------
   cos_sims <- lapply(mut_types, function(mut_type){
      #mut_type='SBS'
      data_ss <- data[grepl(paste0('^', mut_type),names(data))]
      
      do.call(rbind, lapply(names(data_ss), function(sig_name){
         #sig_name='SBS1'
         m <- data_ss[[sig_name]]
         
         if(nrow(m)==0 | ncol(m)==0){ return(NULL) }
         
         cos_sims_1 <- cosSimToCosmic(m, output='raw')
         cos_sims_1 <- t(cos_sims_1)
         
         data.frame(
            #mut_type=mut_type,
            assigned_sig_name=sig_name,
            denovo_sig_name=rownames(cos_sims_1),
            cos_sims_1,
            check.names=F, row.names=NULL
         )
      }))
   })
   names(cos_sims) <- mut_types
   ## Each list object is for a mut type
   
   ## Write excel doc --------------------------------
   if(is.null(export.path)){ return(cos_sims) }
   
   require(openxlsx)
   
   real_data_start_col <- 3
   wb <- createWorkbook()
   for(mut_type in mut_types){
      ## Create sheet for mut type
      #mut_type='SBS'
      addWorksheet(wb, sheetName=mut_type)
      df <- cos_sims[[mut_type]]
      writeData(wb, sheet=mut_type, x=df)
      
      ## Color cells by cos sim values
      conditionalFormatting(
         wb, 
         sheet=mut_type,
         cols = 1:ncol(df), rows = 1:nrow(df),
         style = c('#FCFCFF', '#63BE7B'),
         #rule = c(0, 255),
         type = "colourScale"
      )
      
      ## Make cells narrower
      setColWidths(
         wb,
         sheet=mut_type,
         cols=real_data_start_col:ncol(df),
         widths=5
      )
      
      ## Freeze top row
      freezePane(wb, sheet=mut_type, firstActiveRow=2)
   
      ## Add horizontal guide lines
      first_assigned_sig_name_instance <- which(!duplicated(df$assigned_sig_name))
      addStyle(
         wb, sheet=mut_type, 
         style=createStyle(border='Bottom'), 
         cols=rep(1:ncol(df), length(first_assigned_sig_name_instance)),
         rows=rep(first_assigned_sig_name_instance, each=ncol(df))
      )
   }
   
   saveWorkbook(
      wb, 
      #file=paste0(tables_dir,"/cos_sims--denovo_sigs_vs_cosmic_sigs.xlsx"), 
      file=export.path,
      overwrite=TRUE
   )
}

# ## Cos sim to SBS1 for denovo sigs not assigned to SBS1 --------------------------------
# if(WRITE_OUTPUT){
#    (function(){
#       df <- cosSimToCosmicAllDenovo()$SBS
#       df <- df[,c('assigned_sig_name','denovo_sig_name','SBS1')]
#       colnames(df)[ncol(df)] <- paste0('cos_sim.', colnames(df)[ncol(df)])
#       df <- df[order(df$cos_sim.SBS1, decreasing=T),]
#       df <- subset(
#          df, 
#          assigned_sig_name!='SBS1' & grepl('Breast|Prostate|Kidney|Thyroid', denovo_sig_name)
#       )
#       rownames(df) <- NULL
#       
#       pdf(paste0(tables_dir,'/cos_sim_to_SBS1--denovo_sigs_not_assigned_to_SBS1.pdf'), 5, 4)
#       gridExtra::grid.table(df[1:10,], rows=NULL)
#       dev.off()
#    })()
# }

## Sig hypermutator enrichment ================================
plotSigEnrichmentHypermutators <- function(
   m, min.mut.load=c(SBS=10000, DBS=500, ID=1000), qvalue.thres=0.05, output='plot'
){
   if(F){
      m=etiology_contribs #log10(sig_contribs+1)
      min.mut.load=c(SBS=10000, DBS=500, ID=1000)
      qvalue.thres=0.05
      output='plot'
   }
   
   #t(m['CPCT02010523T',])
   
   ## Init --------------------------------
   require(ggplot2)
   if(!(output %in% c('plot','enr','enr.formatted'))){
      stop("`output` must be 'plot', 'enr', or 'enr.formatted'")
   }
   
   ## Prep input --------------------------------
   df <- data.frame(m, check.names=F)
   
   ## Set low contribs to 0
   mut_types <- stringr::str_extract(colnames(df), '^[A-Z]+')
   min_mut_loads <- min.mut.load[mut_types]
   for(i in 1:ncol(df)){
      #i=1
      df[,i][ df[,i]<min_mut_loads[i] ] <- 0
   }
   
   ## Convert to is_hypermutator TRUE/FALSE
   df <- df>0
   
   ## Add metadata
   df <- cbind( df, getSampleMetadata(rownames(m), c('cancer_type_code','cohort')) )
   df$is_hmf_sample <- df$cohort=='Hartwig'
   
   ## Enrichment --------------------------------
   df_split <- split(df, df$cancer_type_code)
   feature_names <- colnames(m)
   
   enr <- lapply(names(df_split), function(i){
      #i=names(m_split)[1]
      #i='BRCA'
      #print(i)
      df_ss <- df_split[[i]]
      out <- univarFeatSel.default(
         x=df_ss[,feature_names,drop=F],
         y=df_ss[,'is_hmf_sample'],
         alternative='two.sided',
         show.conting=T, show.sample.size=T,
         order.by.pvalue=F
      )
      cbind(cancer_type_code=i, out)
   })
   enr <- do.call(rbind, enr)
   NULL -> enr$is_keep_feature -> enr$is_pass_feature
   
   #enr$feature <- factor(enr$feature, unique(enr$feature))
   #enr$cancer_type_code <- factor(enr$cancer_type_code, unique(enr$cancer_type_code))
   
   ## p-adjust by mut type and cancer type
   enr <- subset(enr, case.true>=5 | ctrl.true>=5) ## Rm these rows to have better qvalue
   enr$mut_type <- stringr::str_extract(enr$feature, '^[A-Z]+')
   enr$mut_type.cancer_type_code <- paste0(enr$mut_type,':',enr$cancer_type_code)
   enr <- do.call(rbind, lapply(split(enr, enr$mut_type.cancer_type_code), function(i){
      i$qvalue <- p.adjust(i$pvalue, method='bonferroni')
      return(i)
   }))
   #enr <- enr[order(enr$mut_type.cancer_type_code, enr$feature),]
   enr$mut_type <- factor(enr$mut_type, c('SBS','DBS','ID'))
   rownames(enr) <- NULL
   
   ## Labels
   #enr$feature <- sub('^\\w+ \\| ','',enr$feature)
   enr$label <- paste0(enr$cancer_type_code, ' | ', enr$feature)
   enr$label[enr$qvalue>qvalue.thres] <- ''
   
   if(output=='enr'){ return(enr) }
   
   if(output=='enr.formatted'){
      sel_cols <- c(
         cancer_type='cancer_type_code',
         sig_group='feature',
         pvalue='pvalue',
         qvalue='qvalue',
         cramer_v='eff_size',
         prop_hypermutator.metastatic='avg_case',
         prop_hypermutator.primary='avg_ctrl',
         n_hypermutator.metastatic='case.true',
         n_metastatic='case.false',
         n_hypermutator.primary='ctrl.true',
         n_primary='ctrl.false'
      )
      enr <- enr[sel_cols]
      colnames(enr) <- names(sel_cols)
      enr$prop_hypermutator.metastatic[is.na(enr$prop_hypermutator.metastatic)] <- 0
      enr$prop_hypermutator.primary[is.na(enr$prop_hypermutator.primary)] <- 0
      return(enr)
   }
   
   ## Plot --------------------------------
   ggplot(enr, aes(x=eff_size, y=-log10(qvalue), color=mut_type)) +
      
      ## Guide lines
      geom_hline(yintercept=0, color='grey') +
      geom_vline(xintercept=0, color='grey') +
      geom_hline(yintercept=-log10(qvalue.thres), color='grey', linetype='dotted') +
      
      ## Points
      geom_point() +
      ggrepel::geom_text_repel(aes(label=label), size=2.7, show.legend=F) +
      
      scale_color_brewer(palette='Set1') +
      
      ## Misc
      labs(x="Effect size (Cramer's V)", color='Mut. type') +
      
      theme_bw() +
      theme(
         panel.grid=element_blank()
      )
}

if(WRITE_OUTPUT){
   p <- plotSigEnrichmentHypermutators(etiology_contribs)
   
   pdf(paste0(plots_dir,'/sig_hypermutator_volcano.pdf'), 9, 7)
   plot(p)
   dev.off()
}

## Supplementary tables ######################################################################################
## Functions ================================
formatMatrix <- function(m, rowname.label='sample_id', anon.sample.ids=T){
   #contribs=sig_contribs
   
   ## Rownames as column 1 --------------------------------
   m <- data.frame(
      rownames = rownames(m),
      m, 
      check.names=F,
      row.names=NULL
   )
   
   colnames(m)[1] <- rowname.label
   
   if(anon.sample.ids){ 
      m[[rowname.label]] <- anonymizeSampleIds(m[[rowname.label]])
   }
   
   return(m)
}

formatSigMetadataTable <- function(sig_metadata){
   
   ## Select relevant columns
   sel_cols <- c(
      sig_name.denovo = 'denovo_name_uniq',
      sig_name.final = 'sig_name',
      is_denovo_clust = 'is_denovo_clust',
      etiology = 'etiology',
      cos_sim.matched_ref = 'matched_ref_cos_sim',
      cos_sim.mean_within_clust_cossim = 'mean_within_clust_cossim',
      cos_sim.merged = 'cos_sim_merged'
   )
   
   sig_metadata_ss <- sig_metadata[sel_cols]
   colnames(sig_metadata_ss) <- names(sel_cols)
   
   ## Add etiology
   sig_metadata_ss$etiology <- getSigEtiology(sig_metadata_ss$sig_name.final)
   
   ## Remove placeholder 'signatures'
   sig_metadata_ss <- subset(sig_metadata_ss, !is.na(sig_name.final))
   
   return(sig_metadata_ss)
}

formatEcdfEnrichmentTable <- function(sample.groups){

   enr <- plotEcdfMutLoad(
      contexts=contexts,
      sample.groups=sample.groups,
      output='enr'
   )
   
   sel_cols <- c(
      cancer_type = 'label_group',
      mut_type = 'feature',
      pvalue = 'pvalue',
      qvalue = 'qvalue',
      median_metastatic = 'avg_case',
      median_primary = 'avg_ctrl',
      fold_change = 'fold_change'
   )

   enr_ss <- enr[sel_cols]
   colnames(enr_ss) <- names(sel_cols)
   rownames(enr_ss) <- NULL

   return(enr_ss)
}

formatSigEnrichmentTable <- function(sample.groups){
   mut_types <- c('SBS','DBS','ID')
   enr <- lapply(mut_types, function(i){
      i_enr <- plotSigEnrichment(
         sample.groups=sample.groups,
         sig.contribs=etiology_contribs,
         mut.type=i,
         output='enr'
      )
      i_enr$mut_type <- i
      return(i_enr)
   })
   enr <- do.call(rbind, enr)
   
   sel_cols <- c(
      cancer_type = 'label_group',
      mut_type = 'mut_type',
      sig_group = 'feature',
      pvalue = 'pvalue',
      qvalue = 'qvalue',
      median_metastatic = 'avg_case',
      median_primary = 'avg_ctrl',
      median_diff = 'median_diff',
      fold_change.p1 = 'fold_change',
      fold_change.log2p1 = 'fold_change.log2',
      is_signif = 'is_signif',
      has_contrib_in_metastatic_or_primary = 'has_contrib'
   )

   enr_ss <- enr[sel_cols]
   colnames(enr_ss) <- names(sel_cols)
   rownames(enr_ss) <- NULL

   return(enr_ss)
}

formatEtiologyGroupDiffsTable <- function(sample.groups){
   #sample.groups=sample_groups$progression_status_code.cancer_type_code
   mut_types <- c('SBS','DBS','ID')
   diffs <- lapply(mut_types, function(i){
      #i='DBS'
      i_plot_data <- plotSigEnrichment(
         sample.groups=sample.groups,
         sig.contribs=etiology_contribs,
         mut.type=i,
         show.neg.median.diff=T,
         output='list.plot_data'
      )
      i_diffs <- i_plot_data$median_diff
      if(all(!i_diffs$is_signif)){ return(NULL) }
      return(i_diffs)
   })
   diffs <- do.call(rbind, diffs)
   
   sel_cols <- c(
      cancer_type = 'label_group',
      mut_type = 'mut_type',
      median_diff = 'median_diff',
      is_signif = 'is_signif',
      etiology = 'etiology',
      etiology_group_1 = 'etiology_group_1',
      etiology_group_2 = 'etiology_group_2',
      median_diff.sign = 'median_diff.sign',
      median_diff.total = 'median_diff.total',
      median_diff.frac = 'median_diff.frac'
   )
   
   diffs_ss <- diffs[sel_cols]
   colnames(diffs_ss) <- names(sel_cols)
   rownames(diffs_ss) <- NULL
   
   return(diffs_ss)
}

## Misc tables ================================
## Denovo sig profiles --------------------------------
if(WRITE_OUTPUT){
   denovo_profiles.export <- lapply(denovo_profiles, function(i){
      formatMatrix(t(i), rowname.label='sig_name.denovo')
   })
   
   openxlsx::write.xlsx(
      denovo_profiles.export,
      paste0(tables_dir,"/smnv--denovo_sig_profiles.xlsx")
   )
}

## Cos sim denovo vs cosmic --------------------------------
if(WRITE_OUTPUT){
   cosSimToCosmicAllDenovo(export.path=paste0(tables_dir,"/smnv--denovo_sigs_vs_cosmic_sigs_cos_sims.xlsx"))
}

## Main ================================
if(WRITE_OUTPUT){
   l <- list()
   
   l$sig_metadata <- formatSigMetadataTable(sig_metadata)
   
   l$sig_etiologies <- sig_etiologies
   
   l$sig_contribs <- formatMatrix(sig_contribs) 
   l$etiology_contribs <- formatMatrix(etiology_contribs) 
   
   l$smnv_burden_enrichment <- formatEcdfEnrichmentTable(sample_groups$cancer_stage.cancer_type)
   
   l$etio_enr.cancer_stage <- formatSigEnrichmentTable(sample_groups$cancer_stage.cancer_type_code)
   l$etio_enr.cancer_stage_subtypes <- formatSigEnrichmentTable(sample_groups$cancer_stage.cancer_type_code_subtype_subset)
   l$etio_enr.metastatic_location <- formatSigEnrichmentTable(sample_groups$metastatic_location.cancer_type_code_subtype_subset)
   l$etio_enr.progression_status <- formatSigEnrichmentTable(sample_groups$progression_status_code.cancer_type_code)
   
   l$etio_diff.cancer_stage <- formatEtiologyGroupDiffsTable(sample_groups$cancer_stage.cancer_type_code)
   l$etio_diff.cancer_stage_subtypes <- formatEtiologyGroupDiffsTable(sample_groups$cancer_stage.cancer_type_code_subtype_subset)
   l$etio_diff.metastatic_location <- formatEtiologyGroupDiffsTable(sample_groups$metastatic_location.cancer_type_code_subtype_subset)
   l$etio_diff.progression_status <- formatEtiologyGroupDiffsTable(sample_groups$progression_status_code.cancer_type_code)

   l$hypermutator_enr <- plotSigEnrichmentHypermutators(etiology_contribs, output='enr.formatted')
   
   openxlsx::write.xlsx(
      l, 
      paste0(tables_dir,"/smnv--mutational_signature_analysis.xlsx")
   )
   
   # ##
   # with(l$smnv_burden_enrichment, {
   #    lapply(split(fold_change, mut_type), function(i){
   #       c(mean=mean(i), sd=sd(i))
   #    })
   # })
   
}



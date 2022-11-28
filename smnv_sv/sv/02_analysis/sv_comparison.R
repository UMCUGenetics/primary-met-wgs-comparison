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
wd <- paste0(base_dir,'/passengers/analysis/sv_comparison/08_20221111')
dir.create(paste0(wd,'/plots/'), showWarnings=F)
dir.create(paste0(wd,'/tables/'), showWarnings=F)

## Dependencies --------------------------------
source(paste0(base_dir,'/passengers/analysis/common_functions/common_functions.R'))

library(ggplot2)
#library(patchwork)

## Load data ================================
## LINX SVs
linx_svs_raw <- cacheAndReadData(
   paste0(base_dir,'/passengers/processed/linx_tables/linx.vis_sv_data.merged.txt.gz')
)

linx_svs <- linx_svs_raw
linx_svs <- linx_svs[linx_svs$sample %in% sample_metadata$sample_id,]
linx_svs$sample <- factor(linx_svs$sample, sample_metadata$sample_id) ## Ensures that samples without SVs are also counted

## Aneuploidy score, LOH, WGD, and TP53 mut enrichment
genome_status <- read.delim(paste0(wd,'/ext_data/genome_status.txt'))
rownames(genome_status) <- genome_status$sample_id
genome_status$sample_id <- NULL

## SV contexts ================================
## DEL, DUP --------------------------------
countDelDupByLen <- function(linx_svs, bin.breaks=c(0,10^(3:7),Inf), resolved.types=c('DEL','DUP'), output='contexts'){
   
   if(!(output %in% c('contexts','raw'))){
      stop("`output` must be 'contexts' or 'raw'")
   }
   
   ## --------------------------------
   df <- subset(
      linx_svs, 
      ResolvedType %in% resolved.types, 
      c(sample, ClusterId, ResolvedType, ChrStart, ChrEnd, PosStart, PosEnd)
   )
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
   m <- table(df_flat$sample, df_flat$context)
   m <- unclass(m)
   return(m)
}

dels <- countDelDupByLen(linx_svs, resolved.types='DEL', bin.breaks=c(0,10000,Inf))
dups <- countDelDupByLen(linx_svs, resolved.types='DUP', bin.breaks=c(0,10000,Inf))

## Complex len --------------------------------
countComplexByLen <- function(linx_svs, bin.breaks=c(0,25,50,100,200,400,Inf), include.inv.trans=T, output='contexts'){
   
   if(!(output %in% c('contexts','raw'))){
      stop("`output` must be 'contexts' or 'raw'")
   }
   
   ## --------------------------------
   df <- linx_svs[,c('sample', 'ClusterId', 'ResolvedType')]
   
   if(!include.inv.trans){
      df <- df[df$ResolvedType=='COMPLEX',]
   } else {
      df <- df[df$ResolvedType=='COMPLEX' | grepl('INV|TRANS',df$ResolvedType),]
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
   
   #complex_len$context <- paste0('COMPLEX_',complex_len$len_bin)
   
   if(output=='raw'){ return(complex_len) }
   
   ## Make context matrix -------------------------------
   m <- table(complex_len$sample, complex_len$len_bin)
   m <- unclass(m)
   colnames(m) <- paste0('COMPLEX_',colnames(m))
   return(m)
}

complex_svs <- countComplexByLen(linx_svs, bin.breaks=c(0,20,Inf), include.inv.trans=T)

## Other SVs --------------------------------
countResolvedTypes <- function(linx_svs){
   
   ## --------------------------------
   df <- linx_svs[,c('sample','ClusterId','ResolvedType')]
   df <- unique(df)
   
   ## --------------------------------
   resolved_types <- names(sort(table(df$ResolvedType), decreasing=T))
   resolved_types_order <- c(
      'LINE',
      'COMPLEX',
      'DEL','DUP',
      grep('INV',resolved_types,value=T),
      grep('TRANS',resolved_types,value=T),
      'DOUBLE_MINUTE',
      'INS',
      resolved_types
   )
   resolved_types_order <- unique(resolved_types_order)
   
   df$ResolvedType <- factor(df$ResolvedType, resolved_types_order)
   
   ## --------------------------------
   m <- table(df$sample, df$ResolvedType)
   m <- unclass(m)
   return(m)
}

all_svs <- countResolvedTypes(linx_svs)

## Profiles for: indel length, SV DEL/DUP length, and COMPLEX cluster size ================================
plotSvLengthDistribution <- function(linx.svs){
   if(F){
      linx.svs = subset(linx_svs, sample %in% sample_groups$cancer_stage.cancer_type$sample_id)
   }
   
   ## Calc lengths --------------------------------
   del_dup_lengths <- countDelDupByLen(linx.svs, output='raw')
   
   complex_lengths <- with(
      subset(linx.svs, ResolvedType=='COMPLEX'),
      aggregate(sample, list(sample=sample, ClusterId=ClusterId), length)
   )
   
   ## Plot data --------------------------------
   prepDensityPlotData <- function(
      df, mut.type=NULL, sample.prop=NULL, seed=1,
      required.cols=c(sample='sample', value='length')
   ){
      if(F){
         df=subset(del_dup_lengths, ResolvedType=='DEL')
         required.cols=c(sample='sample', value='length')
         mut.type='test'
      }
      
      ## Select columns
      pd <- df[required.cols]
      colnames(pd) <- names(required.cols)
      
      ## Normalize density by no. of events per sample
      mut_load <- unclass(table(pd$sample))
      pd$mut_load <- mut_load[ as.character(pd$sample) ]
      pd$weight <- 1/pd$mut_load
      pd$mut_load <- NULL
      
      ## Add mut type
      if(!is.null(mut.type)){ pd$mut_type <- mut.type }
      
      ## Add metadata
      pd$cohort <- getSampleMetadata(pd$sample, 'cohort')
      pd$label_group <- getSampleMetadata(pd$sample, 'cancer_type_code')
      pd$sample <- NULL
      rownames(pd) <- NULL
      
      return(pd)
   }
   
   density_data <- rbind(
      prepDensityPlotData(subset(del_dup_lengths, ResolvedType=='DEL'), mut.type='DEL'),
      prepDensityPlotData(subset(del_dup_lengths, ResolvedType=='DUP'), mut.type='DUP'),
      prepDensityPlotData(complex_lengths, required.cols=c(sample='sample', value='x'), mut.type='COMPLEX')
   )
   density_data$mut_type <- factor(density_data$mut_type, unique(density_data$mut_type))
   
   ## Plot --------------------------------
   ##
   vlines <- rbind(
      c('DEL',10000),
      c('DUP',10000),
      c('COMPLEX',20)
   )
   vlines <- as.data.frame(vlines)
   colnames(vlines) <- c('mut_type', 'xintercept')
   vlines$xintercept <- log10( as.numeric(vlines$xintercept) )
   vlines$mut_type <- factor(vlines$mut_type, levels(density_data$mut_type))
   
   ##
   facet_labels <- c(
      DEL='SV deletion\n(length)',
      DUP='SV duplication\n(length)',
      COMPLEX='Complex SV\n(no. of SVs)'
   )
   
   ##
   ggplot(density_data, aes(x=log10(value), y=label_group)) +
      facet_wrap(~mut_type, nrow=1, labeller=labeller(.cols=facet_labels)) +
      ggridges::geom_density_ridges(aes(height=..density.., weight=weight, color=cohort), stat='density', fill=NA) +
      geom_vline(data=vlines, mapping=aes(xintercept=xintercept), linetype='dotted') +
      
      scale_y_discrete(limits=rev) +
      scale_x_continuous(breaks=0:10, labels=function(x){ 10^x }, position='top') +
      scale_color_manual(
         values=c(Hartwig='#6B469A', PCAWG='#F58134'),
         labels=c(Hartwig='Metastatic', PCAWG='Primary')
      ) +
      
      theme_bw() +
      theme(
         panel.grid.minor=element_blank(),
         axis.title.y=element_blank(),
         axis.text.x=element_text(angle=90, hjust=0, vjust=0.5),
         axis.title.x=element_blank(),
         strip.placement='outside',
         strip.background=element_blank(),
         legend.position='top',
         legend.title=element_blank()
      )
}

if(WRITE_OUTPUT){
   pdf(paste0(wd,'/plots/sv_size.density.pdf'), 8, 8)
   plot(plotSvLengthDistribution(
      subset(linx_svs, sample %in% sample_groups$cancer_stage.cancer_type$sample_id)
   ))
   dev.off()
}

## SV type load ################################################################################
## Prep features ================================
sv_type_load <- (function(){
   df <- data.frame(
      #SV_LOAD=sv_load,
      dels,
      dups,
      complex_svs,
      LINE=all_svs[,'LINE'],
      check.names=F
   )
   
   df <- cbind(rowSums(df), df)
   #df <- df[sample_whitelist,]
   
   return(df)
})()

# if(F){ 
#    write.table(sv_type_load, paste0(wd,'/sv_type_load.txt'), sep='\t', quote=F) 
# }

sv_type_names <- c(
   'Total no. of SVs', 
   'Deletion <10kb', 'Deletion >=10kb', 
   'Duplication <10kb', 'Duplication >=10kb', 
   'Complex <20 SVs', 'Complex >=20 SVs', 
   'LINE'
)
colnames(sv_type_load) <- sv_type_names

sv_type_colors <- c(
   '#D3D3D3',
   "#FB9A99","#ed6e6f",#"#d35456",
   "#B2DF8A","#5bb356",
   "#A6CEE3","#4b93c3",
   "#e1cbac"
)
sv_type_colors <- colorspace::desaturate(sv_type_colors, amount=0.5)
names(sv_type_colors) <- sv_type_names
#scales::show_col(sv_type_colors)
#colorspace::desaturate('#b4815b', amount=0.5) |> scales::show_col()

## SV type S curves ================================
plotEcdfSvTypeLoad <- function(
   sv.type.load, sample.groups,
   
   ## Stat testing
   pvalue.thres=0.01, qvalue.thres=0.05, 
   fc.thres=c(less=0.8, greater=1.2), 
   signif.label.style='show.signif',
   
   ## Misc
   hide.axis.text.x=F, panel.1.rel.height=2, pcawg.sens.loss=NULL, output='plot'
){
   if(F){
      sv.type.load=sv_type_load #sv_type_load[,1,drop=F]
      sample.groups=sample_groups$cancer_type_code.cancer_stage
      
      pvalue.thres=0.01
      qvalue.thres=0.05
      fc.thres=c(less=0.8, greater=1.2)
      signif.label.style='star.signif'
      
      hide.axis.text.x=T
      facet.strip.colors=sv_type_colors
      panel.1.rel.height=2
      output='plot'
      #pcawg.sens.loss=pcawg_sens_loss
   }
   
   ## Init --------------------------------
   require(ggplot2)
   #require(scales)
   
   if(!(output %in% c('plot','enr','averages'))){ stop("`output` must be 'plot', 'enr', or 'averages") }
   
   if(is.vector(sv.type.load)){ 
      sv.type.load <- as.matrix(data.frame(sv.type.load, check.names=F)) 
   }
   
   if(!(signif.label.style %in% c('show.all','show.signif','star.signif'))){ 
      stop("`output` must be 'show.all','show.signif','star.signif'") 
   }
   
   # ##
   # if(!is.null(pcawg.sens.loss)){
   #    sample_cancer_types <- getSampleMetadata(rownames(m),'cancer_type')
   #    sample_cohorts <- getSampleMetadata(rownames(m),'cohort')
   #    
   #    sample_sens_loss <- pcawg.sens.loss['SV']
   #    sample_sens_loss <- sample_sens_loss[sample_cancer_types,,drop=F]
   #    sample_sens_loss <- sample_sens_loss[,rep(1,ncol(m)),drop=F]
   #    rownames(sample_sens_loss) <- NULL
   #    
   #    sample_sens_gain <- 1+sample_sens_loss
   #    sample_sens_gain[sample_cohorts=='Hartwig',] <- 1
   #    
   #    dimnames_m <- dimnames(m)
   #    m <- m*sample_sens_gain
   #    dimnames(m) <- dimnames_m
   # }
   # 
   
   ## Prep plot data --------------------------------
   ## Get data based on sample groups
   pd <- getSampleData(sample.groups, sv.type.load)
   
   ## Force original feature order
   pd$feature <- factor(pd$feature, colnames(sv.type.load))
   
   ## Log transform
   log10p1 <- function(x){ log10(x+1) }
   pd$value.log10p1 <- log10p1(pd$value)
   
   ## Sort by SV load and calculate x-axis indexes --------------------------------
   pd_split <- split(pd, paste0(pd$comparison_group,':',pd$label_group,':',pd$feature))
   range01 <- function(x){(x-min(x))/(max(x)-min(x))}
   pd <- do.call(rbind, lapply(pd_split, function(i){
      i <- i[order(i$value),]
      i$rank <- 1:nrow(i)
      i$rank_norm <- range01(i$rank)
      return(i)
   }))
   rownames(pd) <- NULL
   
   pd <- pd[order(pd$label_group, pd$feature, pd$comparison_group),]
   
   ## Enrichment --------------------------------
   ## Enrichment functions require wide format data
   enr_input <- getSampleData(sample.groups, sv.type.load, out.format='wide')
   
   ## Convert comparison_group to TRUE/FALSE
   enr_input$comparison_bool <- ifelse(
      enr_input$comparison_group==levels(enr_input$comparison_group)[1],
      TRUE, FALSE
   )
   
   ## Calculate enrichment
   enr_input.split <- split(enr_input, enr_input$label_group)
   feature_names <- colnames(sv.type.load)
   
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
   
   ## Force factor level order
   enr$feature <- factor(enr$feature, feature_names)
   enr$label_group <- factor(enr$label_group, levels(pd$label_group))
   
   ## p-adjust by SV type
   enr.split <- split(enr, enr$feature)
   enr <- do.call(rbind, lapply(enr.split, function(i){
      i$qvalue <- p.adjust(i$pvalue, method='bonferroni')
      return(i)
   }))
   
   ## Remove unneeded columns
   NULL -> enr$is_pass_feature -> enr$is_keep_feature
   
   ## Calculate difference in median
   enr$avg_case.log10p1 <- log10p1(enr$avg_case)
   enr$avg_ctrl.log10p1 <- log10p1(enr$avg_ctrl)
   enr$median_diff.log10p1 <- enr$avg_case.log10p1 - enr$avg_ctrl.log10p1
   
   enr$fold_change <- (enr$avg_case+1) / (enr$avg_ctrl+1)
   enr$fold_change.error_direction <- (function(){
      x <- enr$avg_case/enr$avg_ctrl
      out <- rep(0, length(x))
      out[is.na(x)] <- -1
      out[x==Inf] <- 1
      out[enr$avg_case==0 & enr$avg_ctrl==0] <- 0
      return(out)
   })()
   
   ## Filter out not signif values
   enr$median_diff.log10p1.filt <- enr$median_diff.log10p1
   enr$fold_change.filt <- enr$fold_change
   
  
   if(!is.null(qvalue.thres)){
      signif_type <- 'q'
      signif_thres <- qvalue.thres
      which_signif <- enr$qvalue<qvalue.thres & (enr$fold_change <= fc.thres['less'] | enr$fold_change >= fc.thres['greater'])
   } else {
      signif_type <- 'p'
      signif_thres <- pvalue.thres
      which_signif <- enr$pvalue<pvalue.thres & (enr$fold_change <= fc.thres['less'] | enr$fold_change >= fc.thres['greater'])
   }
   
   enr$median_diff.log10p1.filt[!which_signif] <- 0
   enr$fold_change.filt[!which_signif] <- NA
   
   max_eff_size <- max(abs(enr$median_diff.log10p1.filt))
   
   ## For fold change label, add < or > to divide by zero cases. Also round label value
   enr$label <- (function(){
      if(signif.label.style=='show.signif'){
         out <- round(enr$fold_change.filt, 1)
      } else {
         out <- round(enr$fold_change, 1)
      }
      
      inequality <- rep('', nrow(enr))
      inequality[enr$fold_change.error_direction==1] <- '>'
      inequality[enr$fold_change.error_direction==-1] <- '<'
      
      out <- paste0(inequality, out,'x')
      
      if(signif.label.style=='show.signif'){
         out[!which_signif] <- ''
      }
      
      if(signif.label.style=='star.signif'){
         stars <- rep('*', length(out))
         stars[!which_signif] <- ''
         out <- paste0(stars, out)
      }
      
      return(out)
   })()
   
   if(output=='enr'){ 
      out <- enr
      out$label_group <- sub('\n',' ',out$label_group)
      return(out)
   }
   
   ## Averages --------------------------------
   case_name <- levels(pd$comparison_group)[1]
   ctrl_name <- levels(pd$comparison_group)[2]
   
   averages <- reshape2::melt(
      enr[c('label_group', 'feature', 'avg_case.log10p1', 'avg_ctrl.log10p1')],
      measure.vars=c('avg_case.log10p1','avg_ctrl.log10p1')
   )
   averages$comparison_group <- c(
      'avg_case.log10p1'=case_name, 
      'avg_ctrl.log10p1'=ctrl_name
   )[as.character(averages$variable)]
   
   if(output=='avg'){ return(averages) }
   
   ## Plot params --------------------------------
   ## y-axis
   y_labels <- c(0, 10^(1:10))
   y_breaks <- y_labels+1
   y_breaks <- log10(y_breaks)
   
   ## Stripes
   v_stripes <- data.frame( 
      label_group=unique(pd$label_group)
   )
   
   tmp <- rep(c(F,T), nrow(v_stripes))
   v_stripes$has_stripe <- tmp[1:nrow(v_stripes)]
   
   ## Colors
   ## c(Metastatic='#6B469A', Primary='#F58134')
   colors <- structure(
      c('#6B469A','#F58134'),
      name=c(case_name, ctrl_name)
   )
   
   fills <-  c(
      '#FEC565','#F1CA85', ## oranges
      'white',
      '#d5c6de','#a07fb6' ## purples
   )
   
   ## Facet heghts
   facet_heights <- rep(1, length(levels(pd$feature)))
   facet_heights[1] <- panel.1.rel.height
   
   ## Plot --------------------------------
   ggplot(pd) + 
      ## Facets
      facet_grid( 
         if(ncol(sv.type.load)==1){ '~label_group' } else { 'feature~label_group' },
         switch='y'
      ) +
      
      ggh4x::force_panelsizes(rows=facet_heights) +
      
      ## Eff size background
      geom_rect(data=enr, mapping=aes(xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, fill=median_diff.log10p1.filt)) +
      geom_text(data=enr, mapping=aes(y=max(pd$value.log10p1), x=0.5, label=label), size=2, vjust=1) +
      scale_fill_gradientn(colors=fills, limits=c(-max_eff_size, max_eff_size)) +
      
      ## Sorted scatter plots
      geom_point(aes(x=rank_norm, y=value.log10p1, color=comparison_group), size=0.3) +
      scale_color_manual(values=colors) +
      
      ## Average line
      geom_segment(
         data=averages,
         mapping=aes(x=0.2, xend=0.8, y=value, yend=value, color=comparison_group), size=0.8, color='black',
         show.legend=F
      ) +
      geom_segment(
         data=averages, 
         mapping=aes(x=0.2, xend=0.8, y=value, yend=value, color=comparison_group), size=0.5,
         show.legend=F
      ) +
      
      ## Axis lines
      geom_hline(yintercept=-Inf, size=0.8) +
      geom_vline(xintercept=-Inf, size=0.8) +
      
      ##
      scale_y_continuous(breaks=y_breaks, labels=y_labels) +
      guides(
         colour=guide_legend(
            order=1, override.aes=list(size=1.7),
            keyheight=0.7, keywidth=0.7
         ),
         fill=guide_colorbar(
            order=2, direction='horizontal', title.position='top',
            label.position="bottom", label.hjust=0.5, label.vjust=0.5, label.theme=element_text(angle=90, size=8),
            frame.colour='black', ticks.colour='black', barwidth=4, barheight=0.8,
         )
      ) +
      labs(
         color='Cohort',
         fill=paste0("Difference in\nlog10(median)\nif ",signif_type,"<", signif_thres),
         y='Count'
      ) +
      
      theme(
         panel.spacing.x=unit(2,'pt'),
         panel.spacing.y=unit(2,'pt'),
         panel.grid=element_blank(),
         panel.border=element_blank(),
         axis.line.y=element_line(size=0.2),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank(),
         strip.text.x=if(hide.axis.text.x){ element_blank() } else { element_text(angle=90, hjust=0, vjust=0.5, size=8) },
         strip.text.y.left=element_text(angle=0, hjust=1, vjust=0.5),
         strip.placement='outside',
         strip.background.x=element_blank(),
         strip.background.y=element_blank(),
         legend.title=element_text(size=9),
         legend.key=element_rect(fill=NA),
         plot.margin=margin(1,1,1,1)
      )
}

# ## Test
# plotEcdfSvTypeLoad(
#    sv_type_load,
#    sample_groups$cancer_stage.cancer_type_and_code
# )

## Genome status ================================
plotGenomeStatus <- function(
   genome.status, sample.groups,
   pvalue.thres=0.05, grid.ncol=2, hide.axis.text.x=F
){
   if(F){
      genome.status=genome_status
      sample.groups=sample_groups$cancer_type_code.cancer_stage
      pvalue.thres=0.05
      grid.ncol=2
      hide.axis.text.x=T
   }
   
   require(ggplot2)
   
   ## Enrichment --------------------------------
   ## Enrichment functions require wide format data
   enr_input <- getSampleData(sample.groups, genome.status, out.format='wide')
   
   ## Convert comparison_group to TRUE/FALSE
   enr_input$comparison_bool <- ifelse(
      enr_input$comparison_group==levels(enr_input$comparison_group)[1],
      TRUE, FALSE
   )
   
   ## Calculate enrichment
   enr_input.split <- split(enr_input, enr_input$label_group)
   feature_names <- colnames(genome.status)
   
   enr <- lapply(names(enr_input.split), function(i){
      #i=names(df_split)[1]
      #print(i)
      enr_input.ss <- as.data.frame(enr_input.split[[i]])
      out <- univarFeatSel.default(
         x=enr_input.ss[,feature_names,drop=F], 
         y=enr_input.ss$comparison_bool,
         alternative='two.sided', avg.numeric.func='median', order.by.pvalue=F, show.sample.size=T
      )
      cbind(label_group=i, out)
   })
   enr <- do.call(rbind, enr)
   rownames(enr) <- NULL
   
   enr$is_keep_feature <- NULL
   enr$is_pass_feature <- NULL
   enr$is_signif <- enr$pvalue<pvalue.thres
   
   ## Plot --------------------------------
   ## Cancer type order
   enr$label_group <- factor(enr$label_group, levels(enr_input$label_group))
   
   ## Significance direction
   enr$signif_direction <- 0
   enr$signif_direction[enr$is_signif & enr$eff_size>0] <- 1
   enr$signif_direction[enr$is_signif & enr$eff_size<0] <- -1
   enr$signif_direction <- as.factor(enr$signif_direction)
   
   ## Feature labels
   enr$feature_label <- c(
      aneuploidy_score='ANEU',
      loh_prop_in_diploid_samples='LOH',
      has_wgd='WGD',
      has_tp53_mut='TP53'
   )[ enr$feature ]
   
   enr$feature_label <- factor(enr$feature_label, unique(enr$feature_label))
   
   enr$feature_color_label <- enr$feature_label
   enr$feature_color_label[enr$signif_direction!=1] <- NA
   
   ## Grid positions
   grid_template <- matrix(unique(enr$feature_label), ncol=grid.ncol)
   grid_template <- reshape2::melt(grid_template)
   
   enr$xpos <- grid_template$Var1[match(enr$feature_label, grid_template$value)]
   enr$ypos <- -grid_template$Var2[match(enr$feature_label, grid_template$value)]
   
   ## Colors
   n_features <- length(levels(enr$feature_label))
   color_palette <- rep('red4', n_features)
   #color_palette <- c('indianred4', 'red2', 'orangered2','goldenrod4')
   #scales::show_col(color_palette)
   
   ##
   ggplot(enr, aes(x=xpos, y=ypos)) +
      facet_wrap(label_group~., nrow=1) +
      geom_tile(aes(fill=feature_color_label), color='white', size=0.25, alpha=0.2) +
      geom_text(aes(label=feature_label, color=feature_color_label), size=1.8) +
      scale_color_manual(values=color_palette, na.value='lightgrey') +
      scale_fill_manual(values=color_palette, na.value='lightgrey') +
      coord_cartesian(expand=F) +
      theme_bw() +
      theme(
         panel.grid=element_blank(),
         panel.spacing.x=unit(2,'pt'),
         panel.border=element_blank(),
         #panel.border=element_rect(color='lightgrey'),
         legend.position='none',
         axis.title=element_blank(),
         axis.text=element_blank(),
         axis.ticks=element_blank(),
         axis.line=element_line(color=NA),
         strip.text.x=if(hide.axis.text.x){ element_blank() } else { element_text(angle=90, hjust=0, vjust=0.5, size=8) },
         strip.background.x=element_rect(color=NA, fill=NA),
         plot.margin=margin(1,1,1,1)
      )
}

## Test
# plotGenomeStatus(
#    genome_status, sample_groups$cancer_type_code_subtype.cancer_stage
# )

## Rel contrib ================================
plotRelSvTypeLoad <- function(
   sv.type.load, sample.groups,
   sv.load.col.index=1, sv.type.colors=NULL, hide.axis.text.x=F,
   legend.justification=c(0.5, 0.5)
){
   if(F){
      sv.type.load=sv_type_load[,-1]
      sample.groups=sample_groups$cancer_type_code.cancer_stage
      
      sv.load.col.index=1
      sv.type.colors=sv_type_colors
      hide.axis.text.x=F
      legend.justification=c(0.5, 0)
   }
   
   ## Prep plot data --------------------------------
   if(is.vector(sv.type.load)){ 
      sv.type.load <- as.matrix(data.frame(sv.type.load, check.names=F)) 
   }
   
   ## Get data based on sample groups
   pd <- getSampleData(sample.groups, sv.type.load)
   
   ## Force original feature order
   pd$feature <- factor(pd$feature, colnames(sv.type.load))
   
   ##
   agg <- with(pd, {
      out <- aggregate(
         value,
         list(label_group=label_group, comparison_group=comparison_group, feature=feature),
         sum
      )
      colnames(out)[length(out)] <- 'value'
      return(out)
   })
   
   ## --------------------------------
   ggplot(agg, aes(x=comparison_group, y=value)) +
      facet_grid(~label_group) +
      
      geom_bar(aes(fill=feature), stat='identity', position='fill', width=1) +
      geom_vline(xintercept=1.5, size=0.3, color='white') +
      { if(!is.null(sv.type.colors)){ scale_fill_manual(values=sv.type.colors, limits=force) }} +
      
      ## Axis lines
      geom_hline(yintercept=-Inf, size=0.8) +
      geom_vline(xintercept=-Inf, size=0.8) +
      
      ## Scales
      scale_y_continuous(expand=c(0,0), labels=function(x){ paste0(x*100,'%') } ) +
      scale_x_discrete(expand=c(0,0)) +
      guides(
         fill=guide_legend(
            keyheight=0, keywidth=0.5, 
            label.theme=element_text(size=8), 
            override.aes=list(size=3)
         )
      ) +
      
      labs(
         y='% of SVs', 
         x='comparison_group', 
         fill=sprintf('Left bars: %s\nRight bars: %s\n', levels(pd$comparison_group)[1], levels(pd$comparison_group)[2])
      ) +
      
      theme_bw() +
      theme(
         panel.grid=element_blank(),
         panel.spacing.x=unit(3,'pt'),
         panel.border=element_blank(),
         axis.title.y=element_text(size=9, angle=0, hjust=1, vjust=0.5),
         axis.title.x=element_blank(),
         axis.text.x=element_blank(),
         axis.ticks.x=element_blank(),
         axis.line.x=element_blank(),
         axis.line.y=element_line(size=0.25),
         strip.text.x=if(hide.axis.text.x){ element_blank() } else { element_text(angle=90, hjust=0, vjust=0.5, size=8) },
         strip.background.x=element_rect(color=NA, fill=NA),
         strip.placement='outside',
         legend.title=element_text(size=9),
         legend.margin=margin(c(0,5,0,5)),
         legend.justification=legend.justification,
         plot.margin=margin(1,1,1,1)
      )
}

## Combine ================================
plotSummarySvTypeLoad <- function(
   sv.type.load=sv_type_load, 
   genome.status=genome_status, 
   sample.groups,
   export.path=NULL, export.width=NULL, export.height=NULL,
   plot.ecdf.signif.label.style='show.signif'
){
   if(F){
      sv.type.load=sv_type_load
      genome.status=genome_status
      sample.groups=sample_groups$cancer_type_code_subtype_subset.metastatic_location
      export.path=NULL
      export.width=NULL
      export.height=NULL
   }
   
   ## Plots --------------------------------
   suppressWarnings({
      p_genome_status <- plotGenomeStatus(
         genome_status, sample.groups, 
         hide.axis.text.x=F
      )
      
      p_ecdf <- plotEcdfSvTypeLoad(
         sv.type.load, sample.groups, 
         hide.axis.text.x=T,
         signif.label.style=plot.ecdf.signif.label.style
      )
      
      p_sv_frac <- plotRelSvTypeLoad(
         sv.type.load[,-1], sample.groups, 
         hide.axis.text.x=T, sv.type.colors=sv_type_colors, legend.justification=c(0.5, 0)
      )
   })
   
   p_combined <- patchwork::wrap_plots(
      p_genome_status, p_ecdf, p_sv_frac,
      ncol=1, heights=c(1, 12, 3), guides='collect'
   )
   
   # pdf(paste0(wd,'/plots/sv_type_load.ecdf.pdf'), 13, 7)
   # plot(p_combined)
   # dev.off()
   
   ## Auto export dimensions --------------------------------
   if(is.null(export.width)){
      n_label_groups <- length(unique(sample.groups$label_group))
      export.width <- (n_label_groups+6) * 0.5 ## Add 6 to account for legends
   }
   
   if(is.null(export.height)){
      export.height <- 7
   }

   ## Export --------------------------------
   if(is.null(export.path)){ return(p_combined) }
   
   if(!grepl('[.](png|pdf)$',export.path)){
      stop("`export.path` must end with pdf or png")
   }
   
   if(grepl('[.]png$',export.path)){
      png(export.path, export.width, export.height, units='in', res=300)
      plot(p_combined)
      dev.off()
   }
   
   if(grepl('[.]pdf$',export.path)){
      pdf(export.path, export.width, export.height)
      plot(p_combined)
      dev.off()
   }
}


## Linear model, SV load vs confounding factors ################################################################################
## Prep features ================================
sample_ids <- levels(linx_svs$sample)

m_treatment_filled <- (function(){
   
   m <- m_treatment[sample_ids,]
   rownames(m) <- sample_ids
   m[is.na(m)] <- FALSE
   return(m)
})()

variables <- cbind(
   getSampleMetadata(sample_ids, c('ploidy','is_hrd','is_msi')),
   n_treatments=rowSums(m_treatment_filled[sample_ids,]),
   had_radiotherapy=getSampleMetadata(sample_ids, 'had_radiotherapy'),
   m_treatment_filled[sample_ids,],
   m_gene_status[sample_ids,]
)
rownames(variables) <- sample_ids

## Fit linear model ================================
## Functions --------------------------------
lmFit <- function(
   y, x, sample.groups,
   y.log.trans=T,
   min.other.bool.true.frac=0.05, ## Min fraction of TRUE samples, other feature types
   min.treatment.true.frac=0.05, ## Min fraction of TRUE samples, treatment features
   min.driver.true.samples.all=15, ## Min no. of samples with driver mutation, in both Hartwig+PCAWG
   min.driver.true.samples.comparison.group=10, ## Min no. of samples with driver mutation, in Hartwig or PCAWG
   show.warnings=F, verbose=T
){
   # df <- data.frame(
   #    sample=names(sv_load),
   #    sv_load=sv_load,
   #    getSampleMetadata(names(sv_load), c('tissue_group','cancer_type_code','ploidy','is_hmf_sample','had_radiotherapy','is_hrd','is_msi')),
   #    row.names=NULL
   # )
   if(F){
      y=sv_type_load$`Complex >=20 SVs`
      x=variables
      sample.groups=sample_groups$cancer_type_code.progression_status
      # comparison.group.var='cohort'
      # label.group.var='cancer_type_code'
      
      #bool.names='auto'
      min.other.bool.true.frac=0.05
      min.treatment.true.frac=0.05
      #min.driver.true.samples.all=15
      #min.driver.true.samples.comparison.group=10
      
      min.driver.true.samples.all=5
      min.driver.true.samples.comparison.group=3

      
      y.log.trans=T
      show.warnings=F
      lm.type='lm'
      verbose=T
   }
   
   ## Init --------------------------------
   ## Checks
   if(!is.vector(y) & !is.numeric(y)){ stop('`y` must be a numeric vector') }
   if(!is.matrix(x) & !is.data.frame(x)){ stop('`x` must be a matrix or dataframe') }
   
   ## Get sample names
   if(length(names(y))>0){
      sample_names <- names(y)
   } else if(length(rownames(x))>0){
      sample_names <- rownames(x)
   } else {
      stop('`y` must have names or `x` must have rownames')
   }
   
   feature_names <- colnames(x)
   
   ## --------------------------------
   if(verbose){ message('Adding sample metadata...') }
   df <- data.frame(
      y,
      x,
      row.names=rownames(variables), check.names=F
   )
   
   df <- getSampleData(sample.groups, df, out.format='wide')
   
   ## Remove comparison_group sample sizes from label_group
   df$label_group <- sub(' \\(.+\\)$', '', df$label_group)
   df$label_group <- factor2(df$label_group)
   
   ## Fit linear model --------------------------------
   ## 'log10(sv_load+1) ~ ploidy + had_radiotherapy + is_hrd + is_msi'
   feature_names_escaped <- paste0('`',feature_names,'`')
   if(y.log.trans){
      fit_formula <- paste0('log10(y + 1) ~ ',paste(feature_names_escaped, collapse=' + '))
   } else {
      fit_formula <- paste0('y ~ ',paste(feature_names_escaped, collapse=' + '))
   }
   
   main <- function(features, min.treatment.true.frac=0.05, min.driver.true.samples=10, min.other.bool.true.frac=0.05){
      if(F){
         features=df
         min.treatment.true.frac=0.05
         min.driver.true.samples=10
      }
      
      
      fits <- list()
      summ <- lapply(unique(features$label_group), function(i){
         #i='Breast:BRCA'
         #i='Kidney:KIRC'
         
         
         ## Subset for cancer type --------------------------------
         i_df <- features[features$label_group==i,,drop=F]
         i_df <- lapply(i_df, function(column){
            if(all(is.na(column))){
               if(show.warnings){ warning(i,':  had all non-NA values') }
               column <-
                  if(is.logical(column)){
                     FALSE
                  } else if(is.numeric(column)){
                     0
                  } else {
                     ''
                  }
            }
            return(column)
         })
         i_df <- as.data.frame(i_df, check.names=F)
         
         # ##
         # i_df_ss <- subset(i_df, comparison_group=='PCAWG')
         # fit <- lm(formula ='log10(y+1) ~ SETD2', data=i_df_ss)
         # summary(fit)
         # ##
         
         ## Filter boolean features for those with enough TRUE samples --------------------------------
         ## Treatment
         treatment_features <- colnames(m_treatment_filled)
         #colSums(i_df[,treatment_features])
         for(j in treatment_features){
            #j='Platinum'
            j_v <- i_df[,j]
            j_true_frac <- sum(j_v, na.rm=T)/length(j_v)
            if(j_true_frac>=min.treatment.true.frac){ next } ## Dont do anything if enough TRUE samples
            i_df[,j] <- rep(FALSE, length(j_v)) ## Else set everything to FALSE
         }
         #colSums(i_df[,treatment_features])
         
         ## Drivers
         driver_features <- colnames(m_gene_status)
         #colSums(i_df[,driver_features])
         for(j in driver_features){
            #j='BRCA1'
            j_v <- i_df[,j]
            if(sum(j_v, na.rm=T)>=min.driver.true.samples){ next } 
            i_df[,j] <- rep(FALSE, length(j_v)) 
         }
         #colSums(i_df[,driver_features])
         
         ## Other boolean features
         bool_features <- sapply(i_df, is.logical)
         bool_features <- bool_features[names(bool_features) %in% feature_names]
         bool_features <- bool_features[ !(names(bool_features) %in% c(treatment_features, driver_features)) ]
         bool_features <- names(bool_features)
         
         for(j in bool_features){
            #j='had_radiotherapy'
            j_v <- i_df[,j]
            j_true_frac <- sum(j_v, na.rm=T)/length(j_v)
            if(j_true_frac>=min.other.bool.true.frac){ next } 
            i_df[,j] <- rep(FALSE, length(j_v)) 
         }
         
         ## Linear model --------------------------------
         fit <- lm(data=i_df, formula=fit_formula, na.action=na.exclude)
         fits[[i]] <<- fit
         
         ## Extract coef summary stats --------------------------------
         fit_summ <- summary(fit)$coefficients
         colnames(fit_summ) <- c('coef.value','coef.se','coef.tvalue','coef.pvalue')
         
         ## Force all features (i.e. also those with all NAs) to be in the output
         fit_summ <- as.data.frame(fit_summ)
         fit_summ <- cbind(feature=rownames(fit_summ), fit_summ)
         rownames(fit_summ) <- NULL
         fit_summ$feature <- sub('TRUE$','',fit_summ$feature)
         
         ## Add back missing features
         missing_features <- feature_names[!(feature_names %in% fit_summ$feature)]
         fit_summ_missing <- fit_summ[0,]
         fit_summ_missing <- fit_summ_missing[rep(1,length(missing_features)),]
         rownames(fit_summ_missing) <- NULL
         fit_summ_missing$feature <- missing_features
         
         fit_summ_missing$coef.pvalue[is.na(fit_summ_missing$coef.pvalue)] <- 1
         fit_summ_missing[is.na(fit_summ_missing)] <- 0
         
         fit_summ <- rbind(fit_summ, fit_summ_missing)
         fit_summ <- fit_summ[match(feature_names, fit_summ$feature),]
         #fit_summ$coef.qvalue <- p.adjust(fit_summ$coef.pvalue, 'BH')
         
         ## Output --------------------------------
         out <- data.frame(
            label_group=i,
            fit_summ[,c('feature','coef.value','coef.pvalue')],
            row.names=NULL
         )
         
         #out$lm.r_squared_adj=summary(fit)$adj.r.squared
         out$lm.pvalue <- with(
            summary(fit), 
            pf(
               fstatistic[1L], 
               fstatistic[2L],
               fstatistic[3L],
               lower.tail = FALSE
            )
         )
         
         return(out)
      })
      summ <- do.call(rbind, summ)
      summ$feature <- factor(summ$feature, unique(summ$feature))
      
      return(list(
         fits=fits,
         summ=summ
      ))
   }
   
   ## Fit LM on whole dataset and by Hartwig/PCAWG cohorts
   if(verbose){ message('Fitting LM on all samples...') }
   lm_out.all <- main(
      features = df, 
      min.treatment.true.frac = min.treatment.true.frac,
      min.driver.true.samples = min.driver.true.samples.all,
      min.other.bool.true.frac = min.other.bool.true.frac
   )
   
   case_name <- levels(df$comparison_group)[1]
   if(verbose){ message(sprintf('Fitting LM on %s samples...', case_name)) }
   lm_out.case <- main(
      features = subset(df, comparison_group==case_name),
      min.treatment.true.frac = min.treatment.true.frac,
      min.driver.true.samples = min.driver.true.samples.comparison.group,
      min.other.bool.true.frac = min.other.bool.true.frac
   )
   
   ctrl_name <- levels(df$comparison_group)[2]
   if(verbose){ message(sprintf('Fitting LM on %s samples...', ctrl_name)) }
   lm_out.ctrl <- main(
      features = subset(df, comparison_group==ctrl_name),
      min.treatment.true.frac = min.treatment.true.frac,
      min.driver.true.samples = min.driver.true.samples.comparison.group,
      min.other.bool.true.frac = min.other.bool.true.frac
   )
   
   ## --------------------------------
   if(verbose){ message('Formatting output...') }
   summ_merged <- cbind(
      lm_out.all$summ,
      case = lm_out.case$summ[,c('coef.value','coef.pvalue','lm.pvalue')],
      ctrl = lm_out.ctrl$summ[,c('coef.value','coef.pvalue','lm.pvalue')]
   )
   
   lm_out <- list(
      summ=summ_merged,
      fits=lm_out.all$fits,
      variables=list(
         x=df[,feature_names], 
         y=df$y, 
         metadata=df[c('comparison_group','label_group')]
      ),
      params=list(
         y.log.trans=y.log.trans,
         min.treatment.true.frac=min.treatment.true.frac,
         min.driver.true.samples.all=min.driver.true.samples.all,
         min.driver.true.samples.comparison.group=min.driver.true.samples.comparison.group
      )
   )
   class(lm_out) <- c('lm.out',class(lm_out))
   
   return(lm_out)
}

print.lm.out <- function(lm.out){
   #lm.out=lm_out
   cat('Objects in `lm.out`:\n')
   cat( paste0('$',names(lm.out)) )
   cat( '\n\n' )
   cat( '> head(lm.out$summ, n=5)\n' )
   print(head(lm.out$summ, n=5))
}

## Fit LM for all SV types --------------------------------
lmFitAllSvTypes <- function(
   sv.type.load=sv_type_load,
   x=variables,
   sample.groups,
   verbose=T,
   ...
){
   l_lm_out <- lapply(colnames(sv.type.load), function(i){
      #i='Deletion <10kb'
      if(verbose){ message('\n## ',i) }
      lmFit(
         y=sv.type.load[,i], 
         x=variables,
         sample.groups,
         verbose=verbose,
         ...
      )
   })
   names(l_lm_out) <- colnames(sv.type.load)
   return(l_lm_out)
}

#l_lm_out <- lmFitAllSvTypes(sample.groups=sample_groups$cancer_stage.cancer_type_code)

## Plot linear model ================================
## Constants --------------------------------
variable_groups <- (function(){
   out <- structure(rep('Misc', ncol(variables)), names=colnames(variables))
   out['ploidy'] <- 'Aneuploidy'
   out['is_hrd'] <- 'DNA repair defect'
   out['is_msi'] <- 'DNA repair defect'
   
   out[c('n_treatments','had_radiotherapy',colnames(m_treatment_filled))] <- 'Treatment'
   
   out[colnames(m_gene_status)] <- 'Gene mutation'
   out['TP53'] <- 'Aneuploidy'
   
   out <- factor(out, c('Aneuploidy','DNA repair defect','Treatment','Gene mutation','Misc'))
   
   return(out)
})()

# variable_group_colors <- c(
#    'Aneuploidy'="#1B9E77", 
#    'DNA repair defect'="#CD5C5C", 
#    'Treatment'="#7570B3", 
#    'Gene mutation'="#E7298A", 
#    'Misc'='grey'
# )

variable_group_colors <- c(
   'Aneuploidy'="#679a50", 
   'DNA repair defect'="#e93e95",
   'Treatment'="#00668a", 
   'Gene mutation'="#d79b26"#,
   #'Misc'='grey'
)

## Plotting for one SV type --------------------------------
plotLm <- function(
   lm.out, 
   output='volcano',
   lm.pvalue.thres=0.01, ## Linear model pvalue
   lm.coef.min.thres=0, ## Linear model coef
   lm.coef.pvalue.thres=0.01, ## Linear model coef pvalues
   y.pvalue.thres=0.01, ## Cancer types with SV load increase
   
   ## Volcano only
   label.max.overlaps=20, 
   axis.break.y=10, axis.break.y.heights=c(0.2, 1), 
   
   ## Lollipop only
   enr.eff.size.min.thres=0, ## Feature enrichment eff size 
   enr.pvalue.thres=Inf,

   plot.title=NULL,
   feature.group.colors=variable_group_colors,
   seed=1, verbose=T
){
   if(F){
      lm.out=l_lm_out$`Duplication >=10kb`
      #lm.out=l_lm_out$LINE
      
      output='lollipop';
      lm.pvalue.thres=0.01; ## Linear model pvalue
      lm.coef.min.thres=0; ## Linear model coef
      lm.coef.pvalue.thres=0.01; ## Linear model coef pvalues
      y.pvalue.thres=0.01; ## Cancer types with SV load increase
      
      ## Volcano only
      label.max.overlaps=20; 
      axis.break.y=10; axis.break.y.heights=c(0.2, 1); 
      
      ## Lollipop only
      enr.eff.size.min.thres=0; ## Feature enrichment eff size 
      enr.pvalue.thres=1;
      facet.limits=NULL;
      
      plot.title=NULL;
      feature.group.colors=variable_group_colors;
      seed=1; verbose=T
   }
   
   ## https://stats.stackexchange.com/questions/28938/why-do-linear-regression-and-anova-give-different-p-value-in-case-of-consideri
   ## The t-test statistic (and its p-value) is a test of whether b=0
   ## The F-test on the anova() printout is whether the added variable significantly reduces the residual sum of squares.
   
   ## Init --------------------------------
   require(ggplot2)
   require(ggrepel)
   require(ggh4x)
   
   if(!(output %in% c('raw','volcano','lollipop','lollipop.data'))){ 
      stop("`output` must be 'raw', 'volcano', 'lollipop', 'lollipop.data'")
   }
   
   ## Subset for cancer types with high SV load difference ------------------------------------
   all_label_groups <- unique(lm.out$variables$metadata$label_group)
   sel_label_groups <- all_label_groups
   case_name <- levels(lm.out$variables$metadata$comparison_group)[1]
   
   if(!is.null(y.pvalue.thres)){
      df_y <- cbind(
         lm.out$variables$metadata,
         y=lm.out$variables$y
      )
      df_y$is_hmf_sample <- lm.out$variables$metadata$comparison_group==case_name
      
      if(lm.out$params$y.log.trans){ df_y$y <- log10(df_y$y + 1) }
      
      df_y_split <- split(df_y, df_y$label_group)
      enr <- suppressWarnings({
         lapply(names(df_y_split), function(i){
            #i=names(features_split)[1]
            feature_ss <- df_y_split[[i]]
            out <- univarFeatSel.default(
               x=feature_ss['y'],
               y=feature_ss$is_hmf_sample,
               alternative='two.sided',
               avg.numeric.func='median',
               order.by.pvalue=F
            )
            cbind(label_group=i, out)
         })
      })
      enr <- do.call(rbind, enr)
      enr$eff_size <- enr$avg_case - enr$avg_ctrl
      sel_label_groups <- subset(
         enr, 
         pvalue<y.pvalue.thres, 
         label_group, drop=T
      )
      
      if(verbose & length(sel_label_groups)!=length(all_label_groups)){
         message(
            'The following label groups were kept with SV load pvalue<',y.pvalue.thres,':\n',
            paste(sel_label_groups, collapse=', ')
         )
      }
   }
   
   ## Calculate eff sizes, HMF vs PCAWG ------------------------------------
   if(verbose){ message('Calculating feature effect sizes...') }
   features <- data.frame(
      comparison_bool=lm.out$variables$metadata$comparison_group==case_name,
      label_group=lm.out$variables$metadata$label_group,
      lm.out$variables$x,
      check.names=F,
      row.names=NULL
   )
   feature_names <- colnames(lm.out$variables$x)
   
   features_split <- split(features, features$label_group)
   
   enr <- suppressWarnings({
      lapply(sel_label_groups, function(i){
         #i=names(features_split)[1]
         #i='BRCA'
         #print(i)
         feature_ss <- features_split[[i]]
         out <- univarFeatSel.default(
            x=feature_ss[,feature_names,drop=F],
            y=feature_ss$comparison_bool,
            alternative='two.sided',
            avg.numeric.func='median',
            order.by.pvalue=F
         )
         out <- cbind(label_group=i, out)
         
         if(all(is.na(out$avg_ctrl))){ out$avg_ctrl <- out$avg_case }
         
         return(out)
      })
   })
   enr <- do.call(rbind, enr)
   
   ## Linear model summary ------------------------------------
   summ <- subset(lm.out$summ, label_group %in% sel_label_groups)
   summ$label_group <- as.factor(summ$label_group)
   
   ## Remove bad regressions
   uniq_label_groups <- levels(summ$label_group)
   summ$lm.pvalue[is.na(summ$lm.pvalue)] <- 1
   summ <- summ[summ$lm.pvalue < lm.pvalue.thres,]
   
   removed_label_groups <- uniq_label_groups[!(uniq_label_groups %in% unique(summ$label_group))]
   if(verbose & length(removed_label_groups)>0){
      message(
         'The following cancer types were removed for having `lm.pvalue`<',lm.pvalue.thres,':\n  ',
         paste0(removed_label_groups, collapse=', ')
      )
   }
   
   ## 
   if(verbose){ message('Matching LM features with feature effect sizes...') }
   indexes <- match(
      paste0(summ$label_group,':',summ$feature),
      paste0(enr$label_group,':',enr$feature)
   )
   summ <- cbind(
      summ,
      enr[indexes,c('feature_type','pvalue','eff_size_metric','eff_size')]
   )
   summ$eff_size[is.na(summ$eff_size)] <- 0
   
   ## Add feature group for coloring
   summ$feature_group <- variable_groups[ as.character(summ$feature) ]
   
   if(output=='raw'){ return(summ) }
   
   ## Volcano plot ------------------------------------
   if(output=='volcano'){
      
      ## Make labels
      trimString <- function(string, width=17){
         string <- as.character(string)
         out <- strtrim(string, width=width)
         is_long_string <- nchar(string)>=width
         out[is_long_string] <- paste0(out[is_long_string],'...')
         return(out)
      }
      
      summ$label <- sprintf(
         '%s (%s)',
         trimString(summ$feature),
         sub('^.*:','',summ$label_group)
      )
      
      summ$show_label <- with(summ,{
         lm.pvalue < lm.pvalue.thres
         coef.pvalue < lm.coef.pvalue.thres & 
         coef.value > lm.coef.min.thres &
            (
               (#case.lm.pvalue < lm.pvalue.thres & 
                  case.coef.pvalue < lm.coef.pvalue.thres &
                  case.coef.pvalue > lm.coef.min.thres) |
                  
               (#pcawg.lm.pvalue < lm.pvalue.thres & 
                  ctrl.coef.pvalue < lm.coef.pvalue.thres &
                  ctrl.coef.pvalue > lm.coef.min.thres)
            )
      })
      summ$label[!summ$show_label] <- ''
      
      ## Filter eff size
      summ$eff_size[abs(summ$pvalue)>=enr.pvalue.thres] <- 0
      
      ## Force signif points to be plotted last
      summ <-  summ[order(nchar(summ$label), decreasing=F),]
      
      ## Max values
      max_eff_size <- max(abs(summ$eff_size))
      max_lm_coef <- max(abs(summ$coef.value))
      
      ## y-axis break
      if(!is.null(axis.break.y)){
         summ$y_group <- ifelse(-log10(summ$coef.pvalue) >= axis.break.y, 'high', 'low')
      } else {
         summ$y_group <- 'low'
      }
      n_y_groups <- length(unique(summ$y_group))
      
      ## Main
      p <- ggplot(summ, aes(x=coef.value, y=-log10(coef.pvalue))) +
         
         ## Points
         geom_point(
            aes(fill=eff_size, size=show_label), shape=21, stroke=0.4,
            color=ifelse(summ$show_label, 'black', 'grey')#,
            #size=ifelse(summ$show_label, 3, 1.5)
         ) +
         scale_fill_gradient2(
            low='#F58134', mid='white', high='#6B469A',
            limits=c(-max_eff_size, max_eff_size)
         ) +
         scale_size_manual(
            values=c('TRUE'=3, 'FALSE'=1.5), 
            labels=c('TRUE'='True','FALSE'='False'),
            name='Is significant?'
         ) +
         
         ## Labels
         geom_text_repel(
            data=subset(summ, nchar(label)>0),
            mapping=aes(label=label, color=feature_group),
            size=2.7, min.segment.length=0, segment.size=0.2, 
            max.iter=1e7, max.overlaps=label.max.overlaps,
            show.legend=T, seed=seed
         ) +
         scale_color_manual(values=feature.group.colors) +
         
         ## Axes
         coord_cartesian(clip='off') +
         scale_y_continuous(
            minor_breaks=1e-8 ## Use minor breaks to hack in 0 lines
         ) + 
         scale_x_continuous(
            minor_breaks=1e-8,
            limits=c(-max_lm_coef, max_lm_coef)
         ) +
         { if(n_y_groups>1){ facet_grid(y_group~., scales='free_y') } } +
         { if(n_y_groups>1){ force_panelsizes(rows=axis.break.y.heights) } } +
         
         ## Misc
         guides(
            size=guide_legend(order=1, reverse=T, override.aes=list(label="")),
            fill=guide_colorbar(order=2,ticks.colour='black', frame.colour='black', barheight=4.5),
            color=guide_legend(override.aes=list(size=3.5)) 
         ) +
         
         labs(
            x='LM coefficient',
            y='-log10(LM coefficient pvalue)',
            fill='Enrichment\nin metastatic',
            color='Feature group',
            title=if(!is.null(plot.title)){ plot.title } else { waiver() }
         ) +
         
         theme_bw() +
         theme(
            panel.grid=element_blank(),
            panel.grid.minor=element_line(color='grey'),
            panel.border=element_blank(),
            axis.line.x=element_line(size=0.25),
            axis.line.y=element_line(size=0.25),
            strip.background.y=element_blank(),
            strip.text.y=element_blank()
         )
      
      return(p)
   }
   
   ## Lollipop plot ------------------------------------
   if(output %in% c('lollipop', 'lollipop.data')){
      ## Subset for features contributing to increased SV load
      summ <- subset(
         summ,
         lm.pvalue < lm.pvalue.thres &
            coef.pvalue < lm.coef.pvalue.thres & 
            coef.value > lm.coef.min.thres &
            
            eff_size > enr.eff.size.min.thres &
            pvalue < enr.pvalue.thres &
            
            (
               (#case.lm.pvalue < lm.pvalue.thres & 
                  case.coef.pvalue < lm.coef.pvalue.thres &
                  case.coef.pvalue > lm.coef.min.thres) |
                  
               (#ctrl.lm.pvalue < lm.pvalue.thres & 
                  ctrl.coef.pvalue < lm.coef.pvalue.thres &
                  ctrl.coef.pvalue > lm.coef.min.thres)
            )
      )
      
      if(nrow(summ)==0){
         if(output=='lollipop'){
            return(ggplot() + theme_void())
         } else {
            return(list(
               plot = ggplot() + theme_void(),
               nrow = 0,
               coef.max = 0,
               eff_size.max = 0
            ))
         }
      }
      
      # summ <- subset(
      #    summ, 
      #    coef.pvalue < lm.coef.pvalue.thres 
      #    & coef.value > lm.coef.min.thres 
      #    & eff_size > enr.eff.size.min.thres
      #    & pvalue < enr.pvalue.thres
      # )
      
      ## x-axis ordering
      summ <- summ[
         order(
            summ$feature_group,
            summ$feature,
            -summ$coef.value
         )
      ,]
      
      summ$y_id <- 1:nrow(summ)
      #summ$cancer_type_code <- sub('^.+:', '', summ$label_group)
      summ$y_label <- paste0(summ$y_id,'::',summ$label_group)
      summ$y_label <- factor(summ$y_label,unique(summ$y_label))
      
      summ$feature <- factor(summ$feature, unique(summ$feature))
      
      ## Melt
      pd <- summ[,c('feature_group','feature','y_label','coef.value','eff_size')]
      pd <- reshape2::melt(pd, measure.vars=c('coef.value','eff_size'))
      pd$variable <- factor(pd$variable, unique(pd$variable))
      
      ## Colors
      if(is.null(feature.group.colors)){
         feature_group_colors <- RColorBrewer::brewer.pal(8,'Dark2')
      } 
      # else {
      #    facet_labels_y <- as.character(summ$feature_group)[ !duplicated(summ$feature) ]
      #    feature_group_colors <- feature.group.colors[ facet_labels_y ]
      # }
      
      pd$y_label_grouped <- interaction(pd$y_label, pd$feature, pd$feature_group)
      
      #require(facetscales)
      p <- ggplot(pd, aes(y=y_label_grouped, yend=y_label_grouped, x=value, fill=feature_group)) +
         ## Axes
         facet_grid(
            ~variable, scales='free_x',
            labeller=labeller(.cols=c(coef.value='\nLM coefficient', eff_size='Enrichment\nin metastatic'))
         ) +
         guides(y=guide_axis_nested(extend=0.7)) +
         scale_y_discrete(
            limits=rev,
            labels=function(x){ x <- sub('^.+::','',x) }
         ) +
         scale_x_continuous(position='top') +
         
         ## Main
         geom_vline(xintercept=0, size=0.3) +
         geom_segment(
            aes(x=0, xend=value, color=feature_group),
            size=1.5, show.legend=F
         ) +
         geom_point(shape=21, size=3) +
         scale_fill_manual(values=feature.group.colors, drop=FALSE) +
         scale_color_manual(values=feature.group.colors, drop=FALSE) +
         labs(
            fill='Feature group', 
            title=if(is.null(plot.title)){ waiver() } else { plot.title }
         ) +
         
         ## Theme
         theme_bw() +
         theme(
            panel.border=element_blank(),
            panel.grid.major.x=element_blank(),
            panel.grid.minor.x=element_blank(),
            panel.grid.minor.y=element_blank(),
            panel.grid.major.y=element_line(),
            axis.line.y.left=element_line(size=0.3),
            axis.line.x.top=element_line(size=0.3),
            axis.title=element_blank(),
            axis.text.x.top=element_text(angle=30, hjust=0, vjust=0),
            strip.placement='outside',
            strip.background.x=element_blank()
         )
      
      if(output=='lollipop.data'){
         return(list(
            plot = p,
            nrow = nrow(subset(pd, variable=='coef.value')),
            coef.max = max(subset(pd, variable=='coef.value', value)),
            eff_size.max = max(subset(pd, variable=='eff_size', value))
         ))
      }
      
      return(p)
   }
}

#plotLm(l_lm_out$`Complex <20 SVs`, output='lollipop')

## Plotting for all SV types --------------------------------
plotLmAllSvTypes <- function(
   l.lm.out,
   plot.type='volcano',
   patchwork.ncol=2, ## Only for volcano
   export.path=NULL,
   export.width=NULL,
   export.height=NULL,
   verbose=F
){
   if(F){
      l.lm.out <- l_lm_out
      plot.type='lollipop'
      patchwork.ncol=2
      export.path=NULL
      export.width=NULL
      export.height=NULL
      verbose=F
   }

   if(!(plot.type %in% c('volcano','lollipop'))){
      stop("`plot.type` must be 'volcano' or 'lollipop'")
   }

   ## ----------------------------
   if(plot.type=='volcano'){
      ## Main
      lp <- lapply(names(l.lm.out), function(i){
         #i='LINE'
         if(verbose){ message('\n## ', i) }
         plotLm(l.lm.out[[i]], output='volcano', plot.title=i, verbose=F)
      })
      names(lp) <- names(l.lm.out)
      p <- patchwork::wrap_plots(lp, ncol=patchwork.ncol, guides='collect')

      ## Auto export dimensions
      if(is.null(export.width)){
         unit_width <- 5
         export.width <- patchwork.ncol * unit_width
      }

      if(is.null(export.height)){
         unit_height <- 4.25
         patchwork_nrow <- ceiling(length(lp) / patchwork.ncol)
         export.height <- patchwork_nrow * unit_height
      }

      # pdf(paste0(wd,'/plots/volcano.coef_pvalue.coef_value.pdf'), 10, 17)
      # plot(p)
      # dev.off()
   }

   ## ----------------------------
   if(plot.type=='lollipop'){
      ## Main
      lp <- lapply(names(l.lm.out), function(i){
         #i='Deletion >1kb'
         if(verbose){ message('\n## ', i) }
         plotLm(lm.out=l.lm.out[[i]], output='lollipop.data', plot.title=i, verbose=F)
      })
      names(lp) <- names(l.lm.out)

      ## Auto facet limits by adding an invisible point
      coef_max <- max(sapply(lp, function(i) i$coef.max ))
      eff_size_max <- max(sapply(lp, function(i) i$eff_size.max ))

      lp <- lapply(lp, function(i){
         #i=lp[[1]]
         if(i$nrow>0){
            
            facet_limits <- i$plot$data
            facet_limits <- facet_limits[facet_limits$y_label_grouped==facet_limits$y_label_grouped[1],]
            facet_limits$value[facet_limits$variable=='coef.value'] <- coef_max
            facet_limits$value[facet_limits$variable=='eff_size'] <- eff_size_max
            
            # facet_limits <- data.frame(
            #    variable=c('coef.value','eff_size'),
            #    value=c(coef_max, eff_size_max),
            #    y_label_grouped=i$plot$data$y_label[1],
            #    feature_group=i$plot$data$feature_group[1]
            # )
            
            i$plot <- i$plot + geom_point(data=facet_limits, mapping=aes(y=-Inf), show.legend=F, alpha=0)
         }
         return(i)
      })

      ## Combine
      plot_heights <- sapply(lp,`[[`,'nrow')
      p <- patchwork::wrap_plots(
         lapply(lp,`[[`,'plot'),
         heights=plot_heights,
         guides='collect'
      )

      ## Auto export dimensions
      if(is.null(export.width)){
         export.width <- 7
      }

      if(is.null(export.height)){
         #sum(plot_heights + 4) * 0.17
         #export.height <- sum(plot_heights) * 0.25
         export.height <- sum(plot_heights + 4) * 0.17 ## Add 4 to pad for x axis labels and plot title
      }

      # pdf(paste0(wd,'/plots/lollipop.lm_coef.eff_size.pdf'), export.width, export.height)
      # plot(p)
      # dev.off()
   }

   ## Export --------------------------------
   if(is.null(export.path)){
      return(p)
   }

   if(verbose){ message('Exporting combined plots...') }
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

# plotLmAllSvTypes(
#    l_lm_out,
#    plot.type='volcano',
#    export.path=paste0(wd,'/plots/volcano.coef_pvalue.coef_value.pdf')
# )
#
# plotLmAllSvTypes(
#    l_lm_out,
#    plot.type='lollipop',
#    export.path=paste0(wd,'/plots/lollipop.lm_coef.eff_size.pdf')
# )

## SV main figure ################################################################################
formatMatrix <- function(m, rowname.label='sample_id', anon.sample.ids=T){
   #contribs=sig_contribs
   
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

plotSvEnrichment <- function(
   sample.groups,
   export.dir, verbose=T,
   
   lm.min.thres.mode=NULL,
   lm.min.other.bool.true.frac=0.05,
   lm.min.treatment.true.frac=0.05,
   lm.min.driver.true.samples.all=15,
   lm.min.driver.true.samples.comparison.group=10,
   
   plot.ecdf.signif.label.style='star.signif',
   plot.ecdf.export.height=NULL,
   plot.lollipop.height=NULL
){
   
   if(F){
      sample.groups=sample_groups$cancer_stage.cancer_type_code_supertype_subtype
      #export.dir=paste0(wd,'/plots/cancer_type_code_subtype_subset.cancer_stage/')
      verbose=F
      
      lm.min.thres.mode='cancer_type_subsets'
      lm.min.other.bool.true.frac=0.05 ## Min fraction of TRUE samples, other feature types
      lm.min.treatment.true.frac=0.05 ## Min fraction of TRUE samples, treatment features
      lm.min.driver.true.samples.all=15 ## Min no. of samples with driver mutation, in both Hartwig+PCAWG
      lm.min.driver.true.samples.comparison.group=10 ## Min no. of samples with driver mutation, in Hartwig or PCAWG
      
      plot.ecdf.signif.label.style='star.signif'
      plot.lollipop.height=NULL
   }
   
   ## Init ----------------------------
   if(!dir.exists(export.dir)){
      dir.create(export.dir)
   }
   
   ## SV type load summary ----------------------------
   if(verbose){ message('Plotting SV type load summary...') }
   plotSummarySvTypeLoad(
      sample.groups=sample.groups,
      export.path=paste0(export.dir,'/01_ecdf.pdf'),
      export.height=plot.ecdf.export.height,
      plot.ecdf.signif.label.style=plot.ecdf.signif.label.style
   )
   
   ## Fit linear model ----------------------------
   ##
   if(verbose){ message('Fitting linear models...') }
   
   ## Default LM min sample thresholds
   if(!is.null(lm.min.thres.mode)){
      if(lm.min.thres.mode=='cancer_type'){
         lm.min.other.bool.true.frac=0.05
         lm.min.treatment.true.frac=0.05
         lm.min.driver.true.samples.all=15
         lm.min.driver.true.samples.comparison.group=10
         
      } else if(lm.min.thres.mode=='cancer_type_subsets'){
         lm.min.other.bool.true.frac=0.05
         lm.min.treatment.true.frac=0.05
         lm.min.driver.true.samples.all=5
         lm.min.driver.true.samples.comparison.group=3
         
      } else {
         stop("`lm.min.thres.mode` must be 'cancer_type' or 'cancer_type_subsets'")
      }
   }
   
   ## Main
   l_lm_out <- lmFitAllSvTypes(
      sample.groups=sample.groups, verbose=F,
      
      min.other.bool.true.frac=lm.min.other.bool.true.frac, 
      min.treatment.true.frac=lm.min.treatment.true.frac, 
      min.driver.true.samples.all=lm.min.driver.true.samples.all, 
      min.driver.true.samples.comparison.group=lm.min.driver.true.samples.comparison.group
   )
   
   ## Plot linear model output ----------------------------
   if(verbose){ message('Plotting volcano...') }
   plotLmAllSvTypes(
      l_lm_out,
      plot.type='volcano',
      export.path=paste0(export.dir,'/02_volcano.pdf')
   )

   if(verbose){ message('Plotting lollipop...') }
   plotLmAllSvTypes(
      l_lm_out,
      plot.type='lollipop',
      export.path=paste0(export.dir,'/03_lollipop.pdf'),
      export.height=plot.lollipop.height
   )
   
   ## SV load enr table ----------------------------
   if(verbose){ message('Exporting tables...') }
   enr_sv <- plotEcdfSvTypeLoad(sv.type.load=sv_type_load, sample.groups=sample.groups, output='enr')
   
   sel_cols <- c(
      cancer_type = 'label_group',
      sv_type = 'feature',
      pvalue = 'pvalue',
      median_metastatic = 'avg_case',
      median_primary = 'avg_ctrl',
      fold_change.p1 = 'fold_change.filt',
      median_metastatic.log10p1 = 'avg_case.log10p1',
      median_primary.log10p1 = 'avg_ctrl.log10p1',
      median_diff.log10p1 = 'median_diff.log10p1'
   )
   
   enr_sv <- enr_sv[sel_cols]
   colnames(enr_sv) <- names(sel_cols)
   rownames(enr_sv) <- NULL
   
   ## LM table ----------------------------
   ##
   lm_summ <- lapply(names(l_lm_out), function(i){
      #i="Total no. of SVs"
      i_summ <- l_lm_out[[i]]$summ
      cbind(sv_type=i, i_summ)
   })
   lm_summ <- do.call(rbind, lm_summ)
   colnames(lm_summ)[colnames(lm_summ)=='label_group'] <- 'cancer_type'
   
   # ## Assign LM number
   # lm_summ$lm_name <- paste0(lm_summ$sv_type,' -- ',lm_summ$label_group)
   # lm_summ$lm_name <- factor(lm_summ$lm_name, unique(lm_summ$lm_name))
   # lm_summ$lm_index <- as.integer(lm_summ$lm_name)
   # lm_summ$lm_name <- NULL
   
   tables <- list()
   tables$sv_type_burden <- formatMatrix(sv_type_load)
   tables$sv_type_burden_enrichment <- enr_sv
   tables$lm_sv_type_burden_vs_features <- lm_summ
   
   ## Merged tables ----------------------------
   openxlsx::write.xlsx(
      tables, 
      paste0(export.dir,'/analysis.xlsx')
   )
}

## Exec --------------------------------
if(WRITE_OUTPUT){
   ## Cancer type
   plotSvEnrichment(
      sample.groups=sample_groups$cancer_stage.cancer_type_code,
      export.dir=paste0(wd,'/plots/cancer_stage.cancer_type_code/'),
      lm.min.thres.mode='cancer_type'
   )
   
   plotSvEnrichment(
      sample.groups=sample_groups$cancer_stage.cancer_type_code,
      export.dir=paste0(wd,'/plots/cancer_stage.cancer_type_code/'),
      lm.min.thres.mode='cancer_type'
   )
   
   plotSvEnrichment(
      sample.groups=sample_groups$cancer_stage.cancer_type_and_code,
      export.dir=paste0(wd,'/plots/cancer_stage.cancer_type_and_code'),
      lm.min.thres.mode='cancer_type'
   )
   
   ## Cancer subtypes
   plotSvEnrichment(
      sample.groups=sample_groups$cancer_stage.cancer_type_code_subtype_subset,
      export.dir=paste0(wd,'/plots/cancer_stage.cancer_type_code_subtype_subset/'),
      lm.min.thres.mode='cancer_type_subsets',
      plot.lollipop.height=12
   )
   
   plotSvEnrichment(
      sample.groups=sample_groups$cancer_stage.cancer_type_code_supertype_subtype,
      export.dir=paste0(wd,'/plots/cancer_stage.cancer_type_code_supertype_subtype'),
      lm.min.thres.mode='cancer_type'
   )
   
   ## Progression status
   plotSvEnrichment(
      sample.groups=sample_groups$progression_status_code.cancer_type_code,
      export.dir=paste0(wd,'/plots/progression_status_code.cancer_type_code/'),
      lm.min.thres.mode='cancer_type_subsets'
   )
   
   ## Metastatic location
   plotSvEnrichment(
      sample.groups=sample_groups$metastatic_location.cancer_type_code_subtype_subset,
      export.dir=paste0(wd,'/plots/metastatic_location.cancer_type_code_subtype_subset/'),
      lm.min.thres.mode='cancer_type_subsets'
   )
   
   plotSvEnrichment(
      sample.groups=sample_groups$metastatic_location.cancer_type_code,
      export.dir=paste0(wd,'/plots/metastatic_location.cancer_type_code/'),
      lm.min.thres.mode='cancer_type_subsets'
   )
}


## Copy supp tables ################################################################################
if(WRITE_OUTPUT){
   out_dirs <- list.dirs(paste0(wd,'/plots/'), recursive=F)
   tables_dir <- paste0(wd,'/tables/')
   
   for(i in out_dirs){
      #i=out_dirs[[1]]
      origin_path <- paste0(i,'/analysis.xlsx')
      dest_path <- paste0(wd,'/tables/sv--', basename(i), '.xlsx')
      file.copy(origin_path, dest_path)
   }
}













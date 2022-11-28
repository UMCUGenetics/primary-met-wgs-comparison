#' Performs pairwise testing and selects significant features
#'
#' @description For numerical variables, wilcoxon tests are performed. For categorical variables,
#' fisher exact tests are performed. The first factor level is assumed to be the negative outcome,
#' while the other levels are grouped together as the positive outcome. For example,
#' with the factor `as.factor(c('none','loh+pathogenic','deep_deletion'))`, 'none' is considered the
#' negative outcome.
#'
#' When y is a factor (multiclass classification), multiple one-vs-rest pairwise tests (i.e. one for
#' each class label) are performed for each feature. A feature is kept if any of the pairwise tests
#' give a significant pvalue/qvalue.
#'
#' @param x A dataframe of features where rows are samples
#' @param y A vector of class labels for each sample in `x`. For binary classification a logical
#' vector. For multiclass classification a character or factor vector.
#' @param alternative A vector containing 'two.sided','greater' or 'less', corresponding to each
#' feature
#' @param max.pvalue pvalue threshold for keeping features. Default: 0.01
#' @param min.cliff.delta 0 to 1. Cliff delta threshold for keeping continuous features. Applies to
#' both +ve and -ve values
#' @param min.cramer.v 0 to 1. Cliff delta threshold for keeping categorical features. Applies
#' to both +ve and -ve values
#' @param sel.top.n.features Limit the total number of features that are selected
#' @param min.features Minimum number of features to keep. Prevents outputting no features
#' @param show.conting Show contingency matrix for each feature?
#' @param show.avg Show case and control group averages per feature?
#' @param avg.numeric.func Can be 'iqm' (interquartile mean), 'mean', or 'median'
#' @param show.sample.size Show the sample size of case and control groups?
#' @param whitelist A character vector of feature names in x to keep, regardless of statistical
#' enrichment
#' @param order.by.pvalue If FALSE, rows of the output dataframe will be ordered by the same
#' order of the features (i.e. colnames) in `x`. If TRUE, it will be sorted by pvalue
#' @param verbose Show progress messages?
#' @param ... Arguments that can be passed to `univarFeatSel.default()`
#'
#' @return A dataframe containing the pvalue, effect size, and other summary statistics for each
#' feature
#' @export
#'
univarFeatSel <- function(x, y, order.by.pvalue=TRUE, verbose=F, ...){

   if(is.logical(y)){
      univarFeatSel.default(x=x, y=y, order.by.pvalue=order.by.pvalue, verbose=verbose, ...)
   } else if(is.factor(y) | is.character(y)) {
      ## Multiclass one vs rest enrichment
      enr <- do.call(rbind, lapply(unique(y), function(i){
         #i='Breast'
         if(verbose){ message('## ', i) }
         out <- univarFeatSel(x=x, y=y==i, order.by.pvalue=FALSE, verbose>=2, ...)
         cbind(class=i, out)
      }))
      if(order.by.pvalue){ enr <- enr[order(enr$pvalue),] }
      return(enr)
   } else {
      stop('`y` must be a logical, factor, or character vector')
   }
}

#' @rdname univarFeatSel
#' @export
#'
univarFeatSel.default <- function(
   x, y,
   alternative=NULL,
   max.pvalue=0.01, min.cliff.delta=0.1, min.cramer.v=0.1, sel.top.n.features=NULL, min.features=2,
   show.conting=FALSE,
   show.avg=TRUE, avg.numeric.func='iqm',
   show.sample.size=FALSE, whitelist=NULL, order.by.pvalue=TRUE,
   verbose=FALSE
){
   if(F){
      x=features[,-1]
      x=x[,grep('^rmd',colnames(x),invert=T)]
      y=features[,1]=='Skin'

      max.pvalue=0.01
      min.cliff.delta=0.1
      min.cramer.v=0.1
      min.features=2
      verbose=T

      alternative <- 'two.sided'
      #alternative[ grep('(^purple)|(^rmd)',colnames(x)) ] <- 'two.sided'

      #output.type='new.x'
      sel.top.n.features=NULL
   }

   ## Checks --------------------------------
   if(!is.data.frame(x)){ stop('x must be a dataframe') }
   if(!is.logical(y)){ stop('y must be a logical vector') }
   if(length(y)!=nrow(x)){ stop('length(y) does not equal nrow(x)') }
   if(any(sapply(x, is.character))){ stop('characters must be converted to factors') }
   if(is.null(colnames(x))){ stop('x must have colnames') }
   if(!(avg.numeric.func %in% c('iqm','trimmed_mean','mean','median'))){ stop('avg.numeric.func must be iqm, trimmed_mean, mean, or median') }

   ## Init --------------------------------
   ## Alternative vector
   if(is.null(alternative)){
      alternative <- rep('two.sided',ncol(x))
   } else if(length(alternative)==1){
      alternative <- rep(alternative,ncol(x))
   } else {
      if( !all(alternative %in% c('two.sided','greater','less')) ){
         stop("`alternative` contains values other than 'two.sided','greater' or 'less'")
      }

      if(length(alternative)!=ncol(x)){
         stop("`alternative` must be the same length as the number of features")
      }
   }

   if(length(names(alternative))==0){
      names(alternative) <- colnames(x)
   }

   ## Convert data types
   if(verbose){ message('Converting factors to logicals (assuming 1st factor level as negative effect)') }
   x <- lapply(x, function(i){
      if(is.factor(i)){ return(as.integer(i)>1) }
      return(i)
   })

   x <- as.data.frame(x, check.names=F)

   ## Template contingency matrix
   conting_template <- matrix(
      NA, nrow=0, ncol=4,
      dimnames=list(NULL, c('case.true','case.false','ctrl.true','ctrl.false'))
   )

   ## Wilcox test on numeric features --------------------------------
   is_numeric <- sapply(x, is.numeric)
   x_numeric <- x[,is_numeric, drop=F]
   x_numeric <- as.matrix(x_numeric)

   #contingencyMatrix(x=x_numeric, y=y, avg.func='iqm')

   pvalues_numeric <- numeric()
   cliff_delta <- numeric()
   conting_numeric <- conting_template

   if(ncol(x_numeric)!=0){
      if(verbose){ message('Performing wilcox tests for numeric features...') }
      pvalues_numeric <- wilcoxTest.data.frame(
         x_numeric[y,,drop=F],
         x_numeric[!y,,drop=F],
         alternative = alternative[colnames(x_numeric)]
      )
      pvalues_numeric[is.na(pvalues_numeric)] <- 1

      if(verbose){ message("Calculating Cliff's delta for numeric features...") }
      cliff_delta <- cliffDelta.data.frame(x_numeric[y,,drop=F], x_numeric[!y,,drop=F])
   }

   ## Fisher test on logical features --------------------------------
   x_logical <- x[,!is_numeric,drop=F]
   x_logical <- as.matrix(x_logical)

   pvalues_logical <- numeric()
   cramer_v <- numeric()
   conting_logical <- conting_template

   if(ncol(x_logical)!=0){
      if(verbose){ message('Performing fisher tests for logical features...') }
      x_logical <- as.matrix(as.data.frame(x_logical, check.names=F))
      #table(x_logical[,'fusion.TMPRSS2_ERG'])

      conting_logical <- contingencyMatrix(x_logical, y, use.totals=F)
      pvalues_logical <- fisherTest.data.frame(
         conting_logical,
         alternative = alternative[colnames(x_logical)]
      )

      if(verbose){ message("Calculating Cramer's V for logical features...") }
      cramer_v <- cramerV.data.frame(conting_logical)
   }

   ## Aggregate stats from numeric/logical data --------------------------------
   ## Initialized dataframe
   tests <- data.frame(
      feature=c(colnames(x_numeric), colnames(x_logical)),
      feature_type=c(rep('numeric',ncol(x_numeric)), rep('logical',ncol(x_logical)))
   )
   rownames(tests) <- tests$feature

   ## Add test data
   tests$alternative <- alternative[rownames(tests)]
   tests$pvalue <- c(pvalues_numeric, pvalues_logical)
   tests$cliff_delta <- c(cliff_delta, rep(0,ncol(x_logical)))
   tests$cramer_v <- c(rep(0,ncol(x_numeric)), cramer_v)

   ## Merge cliff delta and cramer's v
   which_not_numeric <- tests$feature_type!='numeric'
   tests$eff_size <- tests$cliff_delta
   tests$eff_size[which_not_numeric] <- tests$cramer_v[which_not_numeric]
   tests$eff_size_metric <- 'cliff_delta'
   tests$eff_size_metric[which_not_numeric] <- 'cramer_v'

   ## Add case/ctrl cohort stats  --------------------------------
   if(verbose){ message("Calculating summary stats...") }
   if(show.avg){

      calcAvg <- function(x.numeric, x.logical){
         avg_numeric <- numeric()
         if(ncol(x.numeric)!=0){
            if(avg.numeric.func=='iqm'){
               avg_numeric <- colMeansTrimmed(x.numeric, trim=0.25, na.rm=T)
            } else if(avg.numeric.func=='trimmed_mean'){
               avg_numeric <- colMeansTrimmed(x.numeric, trim=0.1, na.rm=T)
            } else if(avg.numeric.func=='mean'){
               avg_numeric <- colMeans(x.numeric, na.rm=T)
            } else {
               avg_numeric <- matrixStats::colMedians(x.numeric, na.rm=T)
            }
         }

         avg_logical <- numeric()
         if(ncol(x.logical)!=0){
            avg_logical <- colMeans(x.logical, na.rm=T) ## Using colMeans() is equal taking the proportion of TRUE
         }

         c(avg_numeric, avg_logical)
      }

      tests$avg_case <- calcAvg(x_numeric[y,,drop=F], x_logical[y,,drop=F])
      tests$avg_ctrl <- calcAvg(x_numeric[!y,,drop=F], x_logical[!y,,drop=F])
      tests$avg_metric <- avg.numeric.func
      tests$avg_metric[which_not_numeric] <- 'prop'
   }

   if(show.conting){
      tests <- cbind(
         tests,
         rbind(
            contingencyMatrix(x=x_numeric, y=y),
            conting_logical
         )
      )
   }

   ## Show sample size --------------------------------
   if(show.sample.size){
      tests$n_case <- sum(y, na.rm=T)
      tests$n_ctrl <- sum(!y, na.rm=T)
   }

   ## Post-processing --------------------------------
   tests <- tests[order(tests$pvalue),]

   ## Select features
   is_pass_feature <- structure(
      rep(FALSE, nrow(tests)),
      names=tests$feature
   )

   is_pass_feature[
      with(tests,{
         pvalue < max.pvalue &
         (
            (alternative %in% c('two.sided','greater') & cliff_delta >=  min.cliff.delta) |
            (alternative %in% c('two.sided','less')    & cliff_delta <= -min.cliff.delta)
         ) |

         (
            (alternative %in% c('two.sided','greater') & cramer_v >=  min.cramer.v) |
            (alternative %in% c('two.sided','less')    & cramer_v <= -min.cramer.v)
         )
      })
   ] <- TRUE

   ## Return feature names or filtered feature matrix
   keep_features <- names(is_pass_feature)[is_pass_feature]
   if(length(keep_features) < min.features){
      keep_features <- names(is_pass_feature)[1:min.features]
   } else if(!is.null(sel.top.n.features)){
      top_n_features <- min(length(keep_features), sel.top.n.features, ncol(x))
      keep_features <- keep_features[ 1:top_n_features ]
   }

   tests$is_pass_feature <- is_pass_feature
   tests$is_keep_feature <- is_pass_feature

   tests$is_keep_feature[ is_pass_feature & !(tests$feature %in% keep_features) ] <- FALSE

   if(!is.null(whitelist)){
      tests$is_keep_feature[ tests$feature %in% whitelist ] <- TRUE
   }

   ## Cleanup output
   NULL -> rownames(tests) -> tests$cliff_delta -> tests$cramer_v

   if(!order.by.pvalue){
      tests <- tests[match(colnames(x), tests$feature),]
   }

   return(tests)
}



#' Calculate Cliff's delta values (for continuous data)
#'
#' @rdname cliffDelta
#'
#' @description Calculate effect sizes using Cliff's delta values for pairwise continuous variable
#' comparisons (e.g. comparisons that would be done by a wilcox test or t-test)
#'
#' @param x A numeric vector for the 1st group of observations. Alternatively, a matrix, where
#' comparisons will be performed by row.
#' @param y Same as x but for the 2nd group of observations.
#'
#' @return A numeric vector of Cliff's delta values
#' @export
#'
cliffDelta <- function (x, ...) {
   UseMethod("cliffDelta", x)
}

USE_CLIFF_DELTA_R <- FALSE
tryCatch(
   expr = {
      Rcpp::cppFunction('double cliffd(NumericVector x, NumericVector y){
         int len_x = x.size();
         int len_y = y.size();

         int sign_sum = 0;
         for(int i = 0; i < len_x; i++){
            for(int j = 0; j < len_y; j++){
               if(x[i] < y[j]){
                  sign_sum--;
               } else if(x[i] > y[j]){
                  sign_sum++;
               }
            }
         }

         return sign_sum / ((double)len_x*len_y);
      }
      ')
   },
   error = function(e){
      warning('C++ version of cliff delta could not be compiled. Using R version')
      USE_CLIFF_DELTA_R <<- TRUE
   }
)

#' @rdname cliffDelta
#' @method cliffDelta default
#' @export
cliffDelta.default <- function(x, y, use.r.implementation=USE_CLIFF_DELTA_R){

   if(!is.numeric(x) | !is.numeric(y)){ stop('x and y must be numeric matrices') }

   if(!use.r.implementation){
      cliffd(x,y)
   } else {
      signs <- sign(outer(x, y, FUN="-"))
      sum(signs, na.rm=T) / length(signs)
   }
}

#' @rdname cliffDelta
#' @method cliffDelta matrix
#' @export
cliffDelta.matrix  <- function(x, y, use.r.implementation=USE_CLIFF_DELTA_R){

   if(ncol(x)!=ncol(y)){ stop('x and y must have the sample number of columns') }
   x <- as.matrix(x); dimnames(x) <- NULL
   y <- as.matrix(y); dimnames(y) <- NULL

   if(!is.numeric(x) & !is.logical(x)){ stop('x must be a numeric or logical matrix') }
   if(!is.numeric(y) & !is.logical(y)){ stop('y must be a numeric or logical matrix') }

   if(!use.r.implementation){
      ## Cpp implementation
      ## Calculate cliff delta for every col between x and y

      unlist(lapply(1L:ncol(x), function(i){
         #i=1
         cliffd(x[,i],y[,i])
      }), use.names=F)

   } else {
      ## R implementation (~6-7x slower than C++ implementation)
      n_comparisons <- nrow(x)*nrow(y)

      ## Rows of x are compared with cols of y
      ## Therefore need to transpose
      y <- t(y)

      ## For every row of x, calculate how many x values are larger than values in the y matrix
      gt_sums <- rowSums(
         apply(x,1,function(i){
            #i=x[1,]
            rowSums(i > y)
         })
      )

      ## ... same for less than
      lt_sums <- rowSums(
         apply(x,1,function(i){
            rowSums(i < y)
         })
      )

      sign_sums <- gt_sums - lt_sums
      sign_sums / n_comparisons
   }
}

#' @rdname cliffDelta
#' @method cliffDelta data.frame
#' @export
cliffDelta.data.frame <- cliffDelta.matrix

# if(F){
#    set.seed(1)
#    x <- matrix(
#       sample(1:10, 120000,replace=T),
#       ncol=600
#    )
#
#    y <- matrix(
#       sample(2:15, 240000,replace=T),
#       ncol=600
#    )
#
#
#    ##
#    system.time({
#       cliffDelta(x,y,use.r.implementation=F)
#    })
#
#    system.time({
#       cliffDelta(x,y,use.r.implementation=T)
#    })
#
#    ## --------------------------------
#    set.seed(1)
#    x <- matrix(
#       sample(1:10, 120,replace=T),
#       ncol=15
#    )
#
#    y <- matrix(
#       sample(2:15, 120,replace=T),
#       ncol=15
#    )
#
#    ## Orig
#    cliffDelta(x,y,use.r.implementation=T)
#
#    # rowSums(
#    #    apply(x,1,function(i){
#    #       rowSums(i > t(y))
#    #    })
#    # )
#
#    # ## 1
#    # cliffDeltaR <- function(x,y){
#    #    ## Make the vector:
#    #    ## c(row1, row2, row3)
#    #    x2 <- c(t(x))
#    #
#    #    ## Make the vector:
#    #    ## c(row1, row1, row1,
#    #    ##   row2, row2, row2,
#    #    ##   row3, row3, row3)
#    #    y2 <- c(y[rep(1:nrow(y),each=nrow(y)),])
#    #
#    #    ## Use recycling of x2 to compare with y2
#    #    gt_sums <- colSums(matrix(x2>y2, ncol=ncol(x)))
#    #    lt_sums <- colSums(matrix(x2<y2, ncol=ncol(x)))
#    #    (gt_sums - lt_sums) / (nrow(x)*nrow(y))
#    # }
#
#    # ## 2
#    # cliffDeltaR <- function(x,y){
#    #
#    #    nrow_y <- nrow(y)
#    #    gt_lt <- lapply(1L:ncol(x), function(i){
#    #       #i=1
#    #       y_ <- rep(y[,i], each=nrow_y)
#    #       c(
#    #          sum(x[,i] > y_),
#    #          sum(x[,i] < y_)
#    #       )
#    #    })
#    #    gt_lt <- do.call(rbind, gt_lt)
#    #
#    #    (gt_lt[,1] - gt_lt[,2]) / (nrow(x)*nrow(y))
#    # }
#    #
#    # system.time({
#    #    cliffDelta(x,y,use.r.implementation=F)
#    # })
#    #
#    # system.time({
#    #    cliffDeltaR(x,y)
#    # })
#    #
#    # system.time({
#    #    cliffDelta(x,y,use.r.implementation=T)
#    # })
# }

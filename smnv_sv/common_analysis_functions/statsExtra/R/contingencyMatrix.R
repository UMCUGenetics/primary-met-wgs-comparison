#' Create multiple contingency matrices
#'
#' @param x A logical or numeric vector/matrix
#' @param y A logical vector indicating the case/control group (i.e. response)
#' @param use.totals Should case/control totals be used instead of case.false and ctrl.false? Only
#' applies when `x` is logical
#' @param avg.func Can be mean, iqm, or median. When `x` is numeric, the averages in the case and
#' ctrl groups are calculated. The resulting contigency matrix is then: (i) the number of case
#' samples above/below the ctrl average, and (ii) the number of ctrl samples above/below the case
#' average
#' @param show.warnings Show warnings?
#'
#' @return A matrix with the columns: case.true, case.false, ctrl.true, ctrl.false. Each row is
#' thus one contingency matrix
#' @export
#'
contingencyMatrix <- function(x, y, use.totals=FALSE, avg.func=c('mean','iqm','median'), show.warnings=TRUE){

   ## Init --------------------------------
   if(!is.logical(x) & !is.numeric(x)){ stop('x must be a logical or numeric matrix/vector') }
   if(!is.logical(y)){ stop('y must be a logical vector') }

   if(!is.matrix(x)){ x <- matrix(x, ncol=1) }

   case <- x[y,,drop=F]
   ctrl <- x[!y,,drop=F]

   avg.func <- match.arg(avg.func, c('mean','iqm','median'))

   ## Convert input, numeric to logical --------------------------------
   avg_func <- switch(
      avg.func,
      mean =   function(...){ colMeans(..., na.rm=T) },
      iqm =    function(...){ colMeansTrimmed(..., na.rm=T, trim=0.25) },
      median = function(...){ colMedians(..., na.rm=T) }
   )

   if(is.numeric(x)){
      case.avg <- avg_func(case)
      ctrl.avg <- avg_func(ctrl)

      case <- case >= ctrl.avg
      ctrl <- ctrl <  case.avg
   }

   ## Main --------------------------------
   if(anyNA(x)){
      if(show.warnings){
         warning('`x` contains NAs. These were excluded from the contingency counts')
      }
      case.NAs <- apply(case, 2, function(i){ sum(is.na(i)) })
      ctrl.NAs <- apply(ctrl, 2, function(i){ sum(is.na(i)) })

      case[is.na(case)] <- FALSE
      ctrl[is.na(ctrl)] <- FALSE
   } else {
      case.NAs <- 0
      ctrl.NAs <- 0
   }

   case.true <- colSums(case)
   case.false <- nrow(case) - case.true - case.NAs

   ctrl.true <- colSums(ctrl)
   ctrl.false <- nrow(ctrl) - ctrl.true - ctrl.NAs

   if(!use.totals | is.numeric(x)){
      cbind(
         case.true, case.false,
         ctrl.true, ctrl.false
      )
   } else {
      cbind(
         case.true,
         case.total=nrow(case) - case.NAs, ## Use case totals instead
         ctrl.true,
         case.total=nrow(ctrl) - ctrl.NAs ## Use ctrl totals instead
      )
   }
}


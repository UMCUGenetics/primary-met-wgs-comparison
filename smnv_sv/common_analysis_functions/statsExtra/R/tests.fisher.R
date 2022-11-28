#' Fast fisher's exact tests
#'
#' @rdname fisherTest
#'
#' @param case.true Case group (treated) responders
#' @param case.false Case group (treated) non responders
#' @param ctrl.true Control group (untreated) responders
#' @param ctrl.false Control group (untreated) non responders
#' @param alternative Indicates the alternative hypothesis and must be one of "two.sided", "greater"
#' or "less"
#' @param m A matrix containing the 4 columns corresponding to a contingency matrix:
#'
#' @param verbose Show progress?
#' @param ... Arguments that can be passed to fisherTest.default()
#'
#' @return A numeric vector of pvalues
#' @export
#'
fisherTest <- function (x, ...) {
   UseMethod("fisherTest", x)
}

#' @rdname fisherTest
#' @method fisherTest default
#' @export
fisherTest.default <- function(
   case.true, case.false, ctrl.true, ctrl.false,
   alternative='two.sided', verbose=F
){

   ## Based on:
   ## https://stats.stackexchange.com/questions/454248/fisher-test-alternative-problem-in-r
   x <- case.true
   m <- case.true + case.false
   n <- ctrl.true + ctrl.false
   k <- case.true + ctrl.true

   if(verbose){
      pb <- txtProgressBar(max=length(case.true), style=3L)
      counter <- 0L
   }
   unlist(
      Map(function(x,m,n,k,alternative){

         if(verbose){
            counter <<- counter + 1L
            setTxtProgressBar(pb, counter)
         }

         a <- 0L:min(m, k)
         prob <- dhyper(a, m, n, k)

         sel_probs <- switch(
            alternative,
            #two.sided=prob[prob<=prob[a==x]],
            two.sided=prob[prob<=prob[x+1L]],
            greater=prob[a>=x],
            less=prob[a<=x]
         )

         sum(sel_probs)
      }, x,m,n,k,alternative),

      use.names=F
   )
}

#' @rdname fisherTest
#' @method fisherTest matrix
#' @export
fisherTest.matrix <- function(m, ...){
   if(ncol(m)!=4){ stop('Input matrix must have 4 columns corresponding to a flattened contingency matrix') }
   fisherTest.default(
      m[,1L], m[,2L], m[,3L], m[,4L],
      ...
   )
}

#' @rdname fisherTest
#' @method fisherTest data.frame
#' @export
fisherTest.data.frame <- fisherTest.matrix


####################################################################################################
#' Fast chi-squared test for 2x2 matrices (similar to a fisher test)
#'
#' @rdname fisherChi2
#'
#' @param case.true Group 'x' responders
#' @param case.false Group 'x' non responders
#' @param ctrl.true Group 'y' responders
#' @param ctrl.false Group 'y' non responders
#' @param correct a logical indicating whether to apply continuity correction when computing the
#' test statistic for 2 by 2 tables: one half is subtracted from all |O - E| differences; however,
#' the correction will not be bigger than the differences themselves.
#' @param return.statistic Return the raw chi-square statistic
#' @param m A matrix containing the 4 columns corresponding to a contingency matrix
#'
#' @param verbose Show progress?
#' @param ... Arguments that can be passed to fisherChi2.default()
#'
#' @return A numeric vector of pvalues or chi-squared statistics
#' @export
#'
fisherChi2 <- function (x, ...) {
   UseMethod("fisherChi2", x)
}

#' @rdname fisherChi2
#' @method fisherChi2 default
#' @export
fisherChi2.default <- function(case.true, case.false, ctrl.true, ctrl.false, correct=TRUE, return.statistic=FALSE){
   ## Adapted from stats::chisq.test()
   # if(F){
   #    #x=matrix(c(50,60,30,200),nrow=2)
   #    case.true=c(50,50,10)
   #    case.false=c(60,60,20)
   #    ctrl.true=c(30,30,200)
   #    ctrl.false=c(200,200,1000)
   # }

   m <- cbind(case.true, case.false, ctrl.true, ctrl.false)
   n <- rowSums(m)

   sums <- data.frame(
      row1=case.true + ctrl.true,
      row2=case.false + ctrl.false,
      col1=case.true + case.false,
      col2=ctrl.true + ctrl.false
   )

   E <- do.call(rbind, Map(function(row1,row2,col1,col2, n){
      as.vector(outer(c(row1,row2), c(col1,col2), "*")/n)
   }, sums$row1,sums$row2,sums$col1,sums$col2, n))

   m_minus_E <- abs(m - E)

   YATES <- 0
   if(correct){
      YATES <- apply( cbind(0.5, m_minus_E),1,min )
   }

   STATISTIC <- rowSums((m_minus_E - YATES)^2/E)
   if(return.statistic){ return(STATISTIC) }

   pchisq(STATISTIC, 1, lower.tail = FALSE)
}

#' @rdname fisherChi2
#' @method fisherChi2 matrix
#' @export
fisherChi2.matrix <- function(m, ...){
   if(ncol(m)!=4){ stop('Input matrix must have 4 columns corresponding to a flattened contingency matrix') }
   fisherChi2.default(m[,1],m[,2],m[,3],m[,4])
}

#' @rdname fisherChi2
#' @method fisherChi2 data.frame
#' @export
fisherChi2.data.frame <- fisherChi2.matrix


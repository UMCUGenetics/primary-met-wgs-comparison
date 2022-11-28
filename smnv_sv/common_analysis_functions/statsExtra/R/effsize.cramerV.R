#' Calculate Cramer's V
#'
#' @rdname cramerV
#'
#' @description Calculate effect sizes using Cramer's V values for pairwise categorical variable
#' comparisons (e.g. comparisons that would be done by a fisher's exact test)
#'
#' @param case.true Case group responders
#' @param case.false Case group non responders
#' @param ctrl.true Control group responders
#' @param ctrl.false Control group non responders
#' @param show.sign Add `-` to values where effect decreases?
#'
#' @return A numeric vector
#' @export
#'
cramerV <- function (x, ...) {
   UseMethod("cramerV", x)
}

#' @rdname cramerV
#' @method cramerV default
#' @export
cramerV.default <- function(case.true, case.false, ctrl.true, ctrl.false, show.sign=T){
   # if(F){
   #    case.true= c(290 ,60)
   #    case.false=c(289 ,247)
   #    ctrl.true= c(3   ,0)
   #    ctrl.false=c(5034,5309)
   # }

   ## The 2x2 contigency matrix looks like this
   ## case.true   ctrl.true
   ## case.false  ctrl.false

   ## Calculate expect values --------------------------------
   ## Based on this tutorial: http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r
   ## From chisq.test()
   # x <- matrix(c(case.true[1], case.false[1], ctrl.true[1], ctrl.false[1]), nrow=2)
   # n <- sum(x)
   # sr <- rowSums(x)
   # sc <- colSums(x)
   # E <- outer(sr, sc, "*")/n

   ## Vectorized version
   observed <- cbind(case.true, case.false, ctrl.true, ctrl.false)
   n <- rowSums(observed)

   rowsums1 <- case.true  + ctrl.true
   rowsums2 <- case.false + ctrl.false
   colsums1 <- case.true  + case.false
   colsums2 <- ctrl.true  + ctrl.false

   expected <- cbind(
      case.true  = rowsums1*colsums1,
      case.false = rowsums2*colsums1,
      ctrl.true  = rowsums1*colsums2,
      ctrl.false = rowsums2*colsums2
   ) / n

   ## Chi-squared statistic
   chi2 <- rowSums(
      (observed - expected)^2 / expected
   )

   ## Calculate Cramers V --------------------------------
   ## Based on this tutorial: https://www.real-statistics.com/chi-square-and-f-distributions/effect-size-chi-square/
   ## V = sqrt(chi2/(n*df))
   ## where df* = min(r – 1, c – 1) and r = the number of rows and c = the number of columns in the contingency table.
   V <- sqrt(chi2/n)
   V[is.na(V)] <- 0

   if(show.sign){
      ## Add sign based on log(odds ratio)
      observed_1 <- observed + 1
      ratio_case <- observed_1[,'case.true'] / observed_1[,'case.false']
      ratio_ctrl <- observed_1[,'ctrl.true'] / observed_1[,'ctrl.false']
      log_odds <- log2(ratio_case / ratio_ctrl)
      V[log_odds<0] <- -V[log_odds<0]
   }

   return(V)
}

#' @rdname cramerV
#' @method cramerV matrix
#' @export
cramerV.matrix <- function(m){
   if(ncol(m)!=4){ stop('Input matrix must have 4 columns corresponding to a flattened contingency matrix') }
   cramerV.default(m[,1],m[,2],m[,3],m[,4])
}

#' @rdname cramerV
#' @method cramerV data.frame
#' @export
cramerV.data.frame <- cramerV.matrix

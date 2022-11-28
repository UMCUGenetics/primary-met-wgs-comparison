#' Calculate trimmed means by column
#'
#' @param m A numeric matrix
#' @param trim The fraction (0 to 0.5) of observations to be trimmed from each end of `m` before the
#' mean is computed. Values of trim outside that range are taken as the nearest endpoint
#' @param na.rm Should missing values (including NaN) be omitted from the calculations?
#'
#' @return A numeric vector
#' @export
#'
colMeansTrimmed <- function(m, trim=0.25, na.rm=F){
   if(!is.matrix(m) | !is.numeric(m)){ stop('`m` must be a numeric matrix') }
   if(nrow(m)==1){
      #warning('Matrix only contains one sample. Returning the original values')
      return(m[1,])
   }

   m <- apply(m,2,sort)

   n <- nrow(m)
   lo <- floor(n * trim) + 1
   hi <- n + 1 - lo

   m <- m[lo:hi,,drop=F]
   colMeans(m, na.rm=na.rm)
}

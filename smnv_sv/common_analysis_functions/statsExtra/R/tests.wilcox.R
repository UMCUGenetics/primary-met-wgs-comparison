#' Perform multiple wilcox tests
#'
#' @rdname wilcoxTest
#'
#' @description Perform column-wise wilcox tests. Adapted from \
#' https://github.com/karoliskoncevicius/matrixTests/
#'
#' @param x numeric matrix.
#' @param y numeric matrix for the second group of observations.
#' @param alternative alternative hypothesis to use for each row/column of x.
#' A single string or a vector with values for each observation.
#' Values must be one of "two.sided" (default), "greater" or "less".
#' @param mu true values of the location shift for the null hypothesis.
#' A single number or numeric vector with values for each observation.
#' @param exact logical or NA (default) indicator whether an exact p-value
#' should be computed (see Details).
#' A single value or a logical vector with values for each observation.
#' @param correct logical indicator whether continuity correction should be
#' applied in the cases where p-values are obtained using normal approximation.
#' A single value or logical vector with values for each observation.
#' @param pvalue.only Only return pvalues?
#'
#' @return If pvalue.only==TRUE, a numeric vector of pvalues
#' Else, a data.frame where each row contains the results of a wilcoxon test
#' performed on the corresponding row/column of x.
#' The columns will vary depending on the type of test performed.\cr\cr
#' They will contain a subset of the following information:\cr
#' 1. obs.x - number of x observations\cr
#' 2. obs.y - number of y observations\cr
#' 3. obs.tot - total number of observations\cr
#' 4. obs.paired - number of paired observations (present in x and y)\cr
#' 5. statistic - Wilcoxon test statistic\cr
#' 6. pvalue - p-value\cr
#' 7. alternative - chosen alternative hypothesis\cr
#' 8. location.null - location shift of the null hypothesis\cr
#' 9. exact - indicates if exact p-value was computed\cr
#' 10. correct - indicates if continuity correction was performed
#'
wilcoxTest <- function (x, ...) {
   UseMethod("wilcoxTest", x)
}

#' @rdname wilcoxTest
#' @method wilcoxTest matrix
#' @export
wilcoxTest.matrix <- function(
   x, y, alternative="two.sided", mu=0, exact=NA, correct=TRUE,
   pvalue.only=TRUE
){
   ## Function from row_wilcoxon_twosample()
   ## Translate to col_wilcoxon_twosample()
   ## Therefore need to transpose matrices
   x <- t(x)
   y <- t(y)

   force(x)
   force(y)

   # if(is.vector(x))
   #    x <- matrix(x, nrow=1)
   # if(is.vector(y))
   #    y <- matrix(y, nrow=1)

   # if (is.data.frame(x) && all(sapply(x, is.numeric)))
   #    x <- data.matrix(x)
   # if(is.data.frame(y) && all(sapply(y, is.numeric)))
   #    y <- data.matrix(y)

   assert_numeric_mat_or_vec(x)
   assert_numeric_mat_or_vec(y)

   if(nrow(y)==1L & nrow(x)>1L) {
      y <- matrix(y, nrow=nrow(x), ncol=ncol(y), byrow=TRUE)
   }

   assert_equal_nrow(x, y)

   if(length(alternative)==1)
      alternative <- rep(alternative, length.out=nrow(x))
   assert_character_vec_length(alternative, 1, nrow(x))

   choices <- c("two.sided", "less", "greater")
   alternative <- choices[pmatch(alternative, choices, duplicates.ok=TRUE)]
   assert_all_in_set(alternative, choices)

   if(length(mu)==1)
      mu <- rep(mu, length.out=nrow(x))
   assert_numeric_vec_length(mu, 1, nrow(x))
   assert_all_in_open_interval(mu, -Inf, Inf)

   if(length(exact)==1)
      exact <- rep(exact, length.out=nrow(x))
   assert_logical_vec_length(exact, 1, nrow(x))
   assert_all_in_set(exact, c(TRUE, FALSE, NA))

   if(length(correct)==1)
      correct <- rep(correct, length.out=nrow(x))
   assert_logical_vec_length(correct, 1, nrow(x))
   assert_all_in_set(correct, c(TRUE, FALSE))


   hasinfx <- is.infinite(x)
   x[hasinfx] <- NA
   hasinfx <- rowSums(hasinfx) > 0

   hasinfy <- is.infinite(y)
   y[hasinfy] <- NA
   hasinfy <- rowSums(hasinfy) > 0

   nxs  <- rep.int(ncol(x), nrow(x)) - matrixStats::rowCounts(is.na(x))
   nys  <- rep.int(ncol(y), nrow(y)) - matrixStats::rowCounts(is.na(y))

   naexact <- is.na(exact)
   exact[naexact] <- (nxs[naexact] < 50) & (nys[naexact] < 50)

   r <- matrixStats::rowRanks(cbind(x - mu, y), ties.method="average")

   statistic <- rowSums(r[,seq_len(ncol(x)),drop=FALSE], na.rm=TRUE) - nxs * (nxs + 1)*0.5

   nties   <- rowTies(r)
   hasties <- rowSums(nties>0) > 0

   wres <- rep(NA_integer_, nrow(x))
   inds <- exact & !hasties
   wres[inds]  <- do_wilcox_2_exact(statistic[inds], nxs[inds], nys[inds], alternative[inds])
   wres[!inds] <- do_wilcox_2_approx(
      statistic[!inds], nxs[!inds], nys[!inds], alternative[!inds],
      nties[!inds,,drop=FALSE], correct[!inds]
   )

   w1 <- hasinfx
   showWarning(w1, 'had infinite "x" observations that were removed')

   w2 <- hasinfy
   showWarning(w2, 'had infinite "y" observations that were removed')

   w3 <- nxs < 1
   showWarning(w3, 'had less than 1 remaining finite "x" observation')

   w4 <- nys < 1
   showWarning(w4, 'had less than 1 remaining finite "y" observation')

   w5 <- exact & hasties
   showWarning(w5, 'had ties: cannot compute exact p-values with ties')

   statistic[w3 | w4] <- NA

   exact <- exact & !hasties
   correct <- correct & !exact


   rnames <- rownames(x)
   if(!is.null(rnames)) rnames <- make.unique(rnames)

   if(pvalue.only){ return(wres) }

   data.frame(
      obs.x=nxs, obs.y=nys, obs.tot=nxs+nys, statistic=statistic,
      pvalue=wres, alternative=alternative, location.null=mu,
      exact=exact, corrected=correct,
      stringsAsFactors=FALSE, row.names=rnames
   )
}

#' @rdname wilcoxTest
#' @method wilcoxTest data.frame
#' @export
wilcoxTest.data.frame <- wilcoxTest.matrix

#' @rdname wilcoxTest
#' @method wilcoxTest default
#' @export
wilcoxTest.default <- function(x, y, ...){
   x <- matrix(x, ncol=1L)
   y <- matrix(y, ncol=1L)
   wilcoxTest.matrix(x,y, ...)
}

####################################################################################################
## Asserts
assert_numeric_mat_or_vec <- function(x) {
   name <- as.character(substitute(x))
   if(is.null(x) || !is.numeric(x) | !(is.matrix(x) | is.vector(x)))
      stop(paste0('"', name, '"', ' must be a numeric matrix or vector'))
}

assert_logical_vec_length <- function(x, ...) {
   name   <- as.character(substitute(x))
   lens   <- unlist(list(...))
   lnames <- as.character(substitute(list(...)))[-1]
   lnames <- paste(lnames, collapse=' or ')
   if(!(length(x) %in% lens) | !is.logical(x) | (NCOL(x) > 1 & NROW(x) > 1))
      stop(paste0('"', name, '"', ' must be a logical vector with length ', lnames))
}

assert_character_vec_length <- function(x, ...) {
   name   <- as.character(substitute(x))
   lens   <- unlist(list(...))
   lnames <- as.character(substitute(list(...)))[-1]
   lnames <- paste(lnames, collapse=' or ')
   if(!(length(x) %in% lens) | !is.character(x) | (NCOL(x) > 1 & NROW(x) > 1))
      stop(paste0('"', name, '"', ' must be a character vector with length ', lnames))
}

assert_numeric_vec_length <- function(x, ...) {
   name   <- as.character(substitute(x))
   lens   <- unlist(list(...))
   lnames <- as.character(substitute(list(...)))[-1]
   lnames <- paste(lnames, collapse=' or ')
   if(!(length(x) %in% lens) | !is.numeric(x) | (NCOL(x) > 1 & NROW(x) > 1))
      stop(paste0('"', name, '"', ' must be a numeric vector with length ', lnames))
}

assert_all_in_set <- function(x, vals) {
   name <- as.character(substitute(x))
   vnames <- paste(vals, collapse=", ")
   if(is.null(x) | !all(x %in% vals))
      stop(paste0('all "', name, '" values must be in: ', vnames))
}

assert_all_in_open_interval <- function(x, min, max) {
   name <- as.character(substitute(x))
   if(is.null(x) | any(anyNA(x) | x<=min | x>=max))
      stop(paste0('all "', name, '" values must be greater than ', min, ' and lower than ', max))
}

assert_equal_nrow <- function(x, y) {
   namex <- as.character(substitute(x))
   namey <- as.character(substitute(y))
   if(nrow(x) != nrow(y))
      stop(paste0('"', namex, '" and "', namey, '" must have the same number of rows'))
}

## General
showWarning <- function(isWarning, err) {
   if(any(isWarning, na.rm=TRUE)) {
      parentFun <- deparse(as.list(sys.call(-1))[[1]])
      grandFun  <- as.list(sys.call(-2))
      if(length(grandFun) > 0) {
         grandFun <- deparse(grandFun[[1]])
         if(grandFun %in% getNamespaceExports("matrixTests")) {
            parentFun <- grandFun
         }
      }
      pref <- "row"
      if(grepl("^col_", parentFun)) pref <- "column"
      n <- sum(isWarning, na.rm=TRUE)
      i <- match(TRUE, isWarning)
      err <- paste0(parentFun, ": ", n, ' of the ', pref, 's ', err, ".",
                    '\nFirst occurrence at ', pref, ' ', i
      )
      warning(err, call.=FALSE)
   }
}

rowTies <- function(x) {
   dupRows <- apply(x, 1, anyDuplicated, incomparables=NA) != 0
   if(any(dupRows)) {
      dups <- matrix(FALSE, nrow=nrow(x), ncol=ncol(x))
      dups[dupRows,] <- t(apply(x[dupRows,,drop=FALSE], 1, duplicated, incomparables=NA))
      dups <- cbind(which(dups, arr.ind=TRUE), val=x[dups])
      dups <- dups[order(dups[,1], dups[,3]),,drop=FALSE]

      sp <- split(dups[,3], dups[,1])
      sp <- lapply(sp, function(x) rle(x)$length+1)
      cl <- lapply(sp, seq_along)

      dups <- cbind(rep(unique(dups[,1]), lengths(sp)),
                    unlist(cl),
                    unlist(sp)
      )

      res <- matrix(0L, nrow=nrow(x), ncol=max(dups[,2]))
      res[dups[,1:2,drop=FALSE]] <- dups[,3]
   } else {
      res <- matrix(0L, nrow=nrow(x), ncol=1)
   }
   res
}

## Tests
do_wilcox_2_exact <- function(stat, nx, ny, alt) {
   res <- rep(NA_integer_, length(stat))

   case <- stat > (nx*ny*0.5)


   inds <- alt=="two.sided" & case
   if(any(inds)) {
      res[inds] <- stats::pwilcox(stat[inds]-1, nx[inds], ny[inds], lower.tail=FALSE)
      res[inds] <- pmin(2*res[inds], 1)
   }

   inds <- alt=="two.sided" & !case
   if(any(inds)) {
      res[inds] <- stats::pwilcox(stat[inds], nx[inds], ny[inds])
      res[inds] <- pmin(2*res[inds], 1)
   }

   inds <- alt=="greater"
   if(any(inds)) {
      res[inds] <- stats::pwilcox(stat[inds]-1, nx[inds], ny[inds], lower.tail=FALSE)
   }

   inds <- alt=="less"
   if(any(inds)) {
      res[inds] <- stats::pwilcox(stat[inds], nx[inds], ny[inds])
   }

   res
}

do_wilcox_2_approx <- function(stat, nx, ny, alt, nties, correct) {
   res <- rep(NA_integer_, length(stat))

   z <- stat - nx*ny*0.5
   correction <- rep(0, length(stat))
   correction[correct & alt=="two.sided"] <- sign(z[correct & alt=="two.sided"]) * 0.5
   correction[correct & alt=="greater"]   <- 0.5
   correction[correct & alt=="less"   ]   <- -0.5
   z <- z - correction

   sigma <- sqrt((nx*ny/12) * ((nx+ny+1) - rowSums(nties^3 - nties, na.rm=TRUE) / ((nx+ny) * (nx+ny-1))))
   z <- z/sigma


   inds <- alt=="two.sided"
   if(any(inds)) {
      res[inds] <- 2 * pmin(stats::pnorm(z[inds]), stats::pnorm(z[inds], lower.tail=FALSE))
   }

   inds <- alt=="greater"
   if(any(inds)) {
      res[inds] <- stats::pnorm(z[inds], lower.tail=FALSE)
   }

   inds <- alt=="less"
   if(any(inds)) {
      res[inds] <- stats::pnorm(z[inds])
   }

   res
}

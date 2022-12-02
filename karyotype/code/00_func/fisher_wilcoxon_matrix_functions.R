#  ----------------------- Custom functions for matrix calculations of Fishers and Wilcoxon test

# (stolen from Luan)
# calculates cramers V includng the directionality from a contingency table
# for our purpose: 
# case.true = mut_primary, case.false = wt_primary
# ctrl.true = mut_metastatic, ctrl.false = wt_metastatic
cramerV.default <- function(case.true, case.false, ctrl.true, ctrl.false, show.sign=T){
  
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

# function that takes a Nx4 Matrix as input and calculates the Cramers V per row
cramerV.matrix <- function(m){
  if(ncol(m)!=4){ stop('Input matrix must have 4 columns corresponding to a flattened contingency matrix') }
  cramerV.default(m[,1],m[,2],m[,3],m[,4])
}

# (stolen from Luan)
# custom Fishers exact test function which takes the contingency table as input
# for our purpose: 
# case.true = mut_primary, case.false = wt_primary
# ctrl.true = mut_metastatic, ctrl.false = wt_metastatic
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
        two.sided=prob[prob<=prob[x+1L]],
        greater=prob[a>=x],
        less=prob[a<=x]
      )
      
      sum(sel_probs)
    }, x,m,n,k,alternative),
    
    use.names=F
  )
}

# function that takes the Nx4 contingency matrix as input and calculates the Fishers exact test p-value per row
fisherTest.matrix <- function(m, ...){
  if(ncol(m)!=4){ stop('Input matrix must have 4 columns corresponding to a flattened contingency matrix') }
  fisherTest.default(
    m[,1L], m[,2L], m[,3L], m[,4L],
    ...
  )
}

# (stolen from Luan)
# custom Fishers exact test function odds ratio calculation
# for our purpose: 
# case.true = mut_primary, case.false = wt_primary
# ctrl.true = mut_metastatic, ctrl.false = wt_metastatic
oddsRatio.default <- function(case.true, case.false, ctrl.true, ctrl.false, trans.method='none', logistic.growth=0.2){
  
  ## Calc preliminary odds ratios
  x <- case.true/case.false
  y <- ctrl.true/ctrl.false
  
  or <- x/y
  or_inv <- y/x
  
  ## Scale odds ratios < 1 to be negative
  out <- or
  out[out<1] <- NA
  out[is.na(out)] <- -(or_inv[is.na(out)])
  
  ## Deal with 0/0 cases
  out[is.na(out)] <- 0
  
  if(trans.method=='linear'){
    return(out)
  }
  
  if(trans.method=='logistic'){
    return(
      2 / (1 + exp(1)^(-logistic.growth*out) ) - 1
    )
  }
  
  return(or)
}

# function that takes a Nx4 Matrix as input and calculates the Odds ratio of a 2x2 contingency table per row
oddsRatio.matrix <- function(m, ...){
  if(ncol(m)!=4){ stop('Input matrix must have 4 columns corresponding to a flattened contingency matrix') }
  oddsRatio.default(m[,1],m[,2],m[,3],m[,4])
}

# stolen from luan
# columnwise wilcoxon test
# x is a numeric matrix
# y is a numeric matrix
# wilcoxon tests will be performed columnwise: col1 of x compared to col 1 of y and so on.
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

# stolen from luan
# cliff delta effect size for wilcoxon test
USE_CLIFF_DELTA_R <- TRUE

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

#  ----------------------------------------
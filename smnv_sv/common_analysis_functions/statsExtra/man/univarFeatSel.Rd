% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/univarFeatSel.R
\name{univarFeatSel}
\alias{univarFeatSel}
\alias{univarFeatSel.default}
\title{Performs pairwise testing and selects significant features}
\usage{
univarFeatSel(x, y, order.by.pvalue = TRUE, verbose = F, ...)

univarFeatSel.default(
  x,
  y,
  alternative = NULL,
  max.pvalue = 0.01,
  min.cliff.delta = 0.1,
  min.cramer.v = 0.1,
  sel.top.n.features = NULL,
  min.features = 2,
  show.conting = FALSE,
  show.avg = TRUE,
  avg.numeric.func = "iqm",
  show.sample.size = FALSE,
  whitelist = NULL,
  order.by.pvalue = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{x}{A dataframe of features where rows are samples}

\item{y}{A vector of class labels for each sample in `x`. For binary classification a logical
vector. For multiclass classification a character or factor vector.}

\item{order.by.pvalue}{If FALSE, rows of the output dataframe will be ordered by the same
order of the features (i.e. colnames) in `x`. If TRUE, it will be sorted by pvalue}

\item{verbose}{Show progress messages?}

\item{...}{Arguments that can be passed to `univarFeatSel.default()`}

\item{alternative}{A vector containing 'two.sided','greater' or 'less', corresponding to each
feature}

\item{max.pvalue}{pvalue threshold for keeping features. Default: 0.01}

\item{min.cliff.delta}{0 to 1. Cliff delta threshold for keeping continuous features. Applies to
both +ve and -ve values}

\item{min.cramer.v}{0 to 1. Cliff delta threshold for keeping categorical features. Applies
to both +ve and -ve values}

\item{sel.top.n.features}{Limit the total number of features that are selected}

\item{min.features}{Minimum number of features to keep. Prevents outputting no features}

\item{show.conting}{Show contingency matrix for each feature?}

\item{show.avg}{Show case and control group averages per feature?}

\item{avg.numeric.func}{Can be 'iqm' (interquartile mean), 'mean', or 'median'}

\item{show.sample.size}{Show the sample size of case and control groups?}

\item{whitelist}{A character vector of feature names in x to keep, regardless of statistical
enrichment}
}
\value{
A dataframe containing the pvalue, effect size, and other summary statistics for each
feature
}
\description{
For numerical variables, wilcoxon tests are performed. For categorical variables,
fisher exact tests are performed. The first factor level is assumed to be the negative outcome,
while the other levels are grouped together as the positive outcome. For example,
with the factor `as.factor(c('none','loh+pathogenic','deep_deletion'))`, 'none' is considered the
negative outcome.

When y is a factor (multiclass classification), multiple one-vs-rest pairwise tests (i.e. one for
each class label) are performed for each feature. A feature is kept if any of the pairwise tests
give a significant pvalue/qvalue.
}

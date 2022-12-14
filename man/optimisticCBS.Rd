% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ChangePoints.R
\name{optimisticCBS}
\alias{optimisticCBS}
\title{Optimistic Circular Binary Segmentation}
\usage{
optimisticCBS(
  X,
  optimization = "L2",
  method = "advanced",
  select = "mixed",
  alpha = 0.01,
  nr_perms = 10000,
  bd = NULL,
  min_seg = 2
)
}
\arguments{
\item{X}{data vector or matrix}

\item{optimization}{minimization of \code{"L1"} or \code{"L2"} loss}

\item{method}{chose between the fast \code{advanced} or the complete \code{full} search}

\item{select}{Model selection approach "perm" (permutation test only), "mixed" (first Wilcoxon test, then permutation test), "wilc" (Wilcoxon test only)}

\item{alpha}{Significance level for testing change points}

\item{nr_perms}{Number of permutation for permutation tests}

\item{bd}{stopping boundary for early stopping}

\item{min_seg}{minimal segment length. Note that change points still can be closer together than \code{min_seg}}
}
\value{
Integer vector containing the change points found in the data
}
\description{
This is the main function of the package. It uses \code{CheckValid}
to find the best change point candidates and then calls itself recursively to find
further change points on subsegments until no more significant change points are found.
}
\examples{
set.seed(1)
X <- rep(c(0,1,0),c(100,20,80)) + rnorm(200,0,0.6)
optimisticCBS(X)
}

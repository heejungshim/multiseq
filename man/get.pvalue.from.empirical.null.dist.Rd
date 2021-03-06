% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get.pvalue.from.empirical.null.dist}
\alias{get.pvalue.from.empirical.null.dist}
\title{Compute a p-value using the empirical distribution of test statistic under the null.}
\usage{
get.pvalue.from.empirical.null.dist(
  statistic.null,
  observed.statistic,
  big.sig = TRUE,
  seed = NULL
)
}
\arguments{
\item{statistic.null}{a vector of test statistics under the null.}

\item{observed.statistic}{a vector of observed test statistics.}

\item{big.sig}{bool indicating whether bigger statistic is more significant.}

\item{seed}{a number; seed number.}
}
\value{
a vector of p-values corresponding to the observed test statistics.
}
\description{
This function takes empirical distribution of statistic under the null (\code{statistic.null}) and 
a series of observed test statistics (\code{observed.statistic}), and returns a series of p-values corresponding 
to the observed test statistics. The resulting p-values will be often used in the  `qvalue` R package (Storey et al., 2020) 
to compute the Storey’s FDR (Storey, 2011). The p-values computed from the empirical null distribution are not always 
continuous, potentially leading to the inaccurate estimation of FDR in the `qvalue` package which uses the insignificant p values 
to estimate the proportion of null tests. So this function uses randomization to produce continuous p-values that are uniformly 
distributed under the null. Let t be observed test statistic. p-value is defined by P(T >= t | H0), but this is not uniformly 
distributed unless test statistic is continuous. We randomize them using P(T > t | H0) + U*P(T = t | H0), where U ~ uniform(0,1). 
The proposed p-value is between P(T > t | H0) and P(T > t | H0) + P(T = t | H0). This p-value can be empirically computed 
by #[T > t] + U*(#[T = t] + 1) / (total number of statistic under the null + 1). Here `+1`` in the denominator is due to observed t.
}

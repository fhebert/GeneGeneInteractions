\name{OmniGGInterTest}
\alias{OmniGGInterTest}
\title{
  Omnibus gene-gene interaction test
}
\description{
\code{OmniGGInterTest} can be used to perform an omnibus gene-gene interaction test. Given two genes and a binary phenotype, two interaction models are constructed (one considering a continuous modelization of the SNPs and another considering a dummy coding) and corresponding test statistics vectors are computed. On each vector, the minP, higher criticism, L2-norm and Hotelling methods are applied. Finally, an overall p-value is computed by combining the 8 tests using Simes' method.
}
\usage{
OmniGGInterTest(X1, X2, Y, N = 1000, perm.method = "parametric", 
                test.method = c("minP", "HC", "L2norm", "Hotelling"),U = NULL)
}
\arguments{
  \item{X1}{a n x p1 matrix of genotypes (coded as 0/1/2)}
  \item{X2}{a n x p2 matrix of genotypes (coded as 0/1/2)}
  \item{Y}{a vector of length n containing the values of the binary phenotype (0/1)}
  \item{N}{the desired number of permutations of the phenotype (default = 1000)}
  \item{perm.method}{a character giving the chosen resampling method for the phenotype: paramatric bootstrap (\code{perm.method="parametric"}) or simple permutations (any other character)}
  \item{test.method}{a vector giving the methods to be combined in the omnibus test, among \code{"minP"}, \code{"HC"}, \code{"L2norm"} and \code{"Hotelling"}. By default, the four methods are included}
  \item{U}{an optional n x q matrix of covariates}
}
\value{
  A list containing the following elements:
  \item{res_omni}{the omnibus p-value, obtained by combining the minP, higher criticism, L2-norm 
  and Hotelling's tests using the continuous and dummy codings}
  \item{res_omni_MGFR}{the p-value of an omnibus test combining the MGF-R test using the continuous and dummy codings}
  \item{res_continuous}{the p-values of the minP, higher criticism, L2-norm, Hotelling and MGF-R methods 
  using the continuous coding}
  \item{res_dummy}{the p-values of the minP, higher criticism, L2-norm, Hotelling and MGF-R methods
  using the dummy coding}
  \item{p0_continuous}{the resampling-based p-values of the minP, higher criticism, L2-norm, Hotelling and MGF-R
  methods using the continuous coding}
  \item{p0_dummy}{the resampling-based p-values of the minP, higher criticism, L2-norm, Hotelling and MGF-R 
  methods using the dummy coding}
  \item{fit_continuous}{a list containing the vector of test statistics, corresponding correlation matrix and vectors of test statistics obtained on permutations using the continuous coding}
  \item{fit_dummy}{a list containing the vector of test statistics, corresponding correlation matrix and vectors of test statistics obtained on permutations using the dummy coding}
}

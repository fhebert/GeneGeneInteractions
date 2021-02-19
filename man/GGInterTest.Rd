\name{GGInterTest}
\alias{GGInterTest}
\title{
  Gene-gene interaction test
}
\description{
\code{GGInterTest} can be used to perform global gene-gene interaction tests. Given two genes and a binary phenotype, two interaction models are constructed (one considering a continuous modelization of the SNPs and another considering a dummy coding) and corresponding test statistics vectors are computed. On each vector, the minP, higher criticism, L2-norm, Hotelling and MGF-R tests are applied. Finally, the p-value corresponding to each test is given.
}
\usage{
GGInterTest(X1, X2, Y, N = 1000, perm.method = "parametric", 
inter.mod = c("continuous", "dummy"), U = NULL)
}
\arguments{
  \item{X1}{a n x p1 matrix of genotypes (coded as 0/1/2)}
  \item{X2}{a n x p2 matrix of genotypes (coded as 0/1/2)}
  \item{Y}{a vector of length n containing the values of the binary phenotype (0/1)}
  \item{N}{the desired number of permutations of the phenotype (default = 1000)}
  \item{perm.method}{a character giving the chosen resampling method for the phenotype: paramatric bootstrap (\code{perm.method="parametric"}) or simple permutations (any other character)}
  \item{inter.mod}{a character giving the type of coding to use to modelize the interaction effect. 
  Possible values are: \code{"continuous"}, \code{"dummy"}, \code{c("continuous","dummy")} 
  (default: computes both models)}
  \item{U}{an optional n x q matrix of covariates}
}
\value{
  A list containing the following elements:
  \item{res_continuous}{the p-values of the minP, higher criticism L2-norm, Hotelling and MGF-R tests using the continuous coding}
  \item{res_dummy}{the p-values of the minP, higher criticism L2-norm, Hotelling and MGF-R tests using the dummy coding}
}

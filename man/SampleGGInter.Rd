\name{SampleGGInter}
\alias{SampleGGInter}
\title{Simulation of samples under a gene-gene interaction model}
\usage{
SampleGGInter(X1, X2, beta0, beta1, I.beta1, mod.beta1 = rep("A",length(I.beta1)),
              beta2, I.beta2, mod.beta2 = rep("A",length(I.beta2)), gamma, I.gamma, 
              n0, n1, U = NULL, betaU = NULL)
}
\arguments{
\item{X1}{a N x p1 matrix of genotypes (coded as 0/1/2) for a given gene for a population}

\item{X2}{a N x p2 matrix of genotypes (coded as 0/1/2) for a given gene for a population}

\item{beta0}{the intercept of the logistic model}

\item{beta1}{the vector of non-zero parameters for the main effects for the first gene}

\item{I.beta1}{the vector of coordinates corresponding to the coefficients in \code{beta1}}

\item{mod.beta1}{a vector of characters specifying the genetic association model for each SNP position given in \code{I.beta1}. Possibles values: \code{"A"} (additive model), \code{"D"} (dominant) and \code{"R"} (recessive)}

\item{beta2}{the vector of non-zero parameters for the main effects for the second gene}

\item{I.beta2}{the vector of coordinates corresponding to the coefficients in \code{beta2}}

\item{mod.beta2}{a vector of characters specifying the genetic association model for each SNP position given in \code{I.beta2}. Possibles values: \code{"A"} (additive model), \code{"D"} (dominant) and \code{"R"} (recessive)}

\item{gamma}{a K x 4 matrix or a vector of length K of gene-gene interaction effects parameters. 
If it is a K x 4 matrix, each row is a vector of parameters for a pair of SNPs having an interaction 
effect and a dummy coding of the SNPs is used. The four parameters in a row correspond to the 
heterozygous/heterozygous, heterozygous/homozygous, homozygous/heterozygous and homozygous/homozygous 
configurations, respectively. If it is a vector, each coordinate is the interaction parameter for 
a pair of SNPs having an interaction effect. In this case, a continuous coding is used and the 
interaction effect is multiplicative}

\item{I.gamma}{a K x 2 matrix of coordinates for the interacting SNPs. Each row is a pair of SNPs. For a given row \code{(i,j)}, the i-th SNP of the first gene and the j-th SNP of the second gene have an interaction effect on the phenotype. The interaction parameters for the pair of SNPs in the k-th row of \code{I.gamma} are given by the k-th row (or coordinate) of \code{gamma}}

\item{n0}{the desired number of controls}

\item{n1}{the desired number of cases}

\item{U}{an optional n x q matrix of covariates}

\item{betaU}{an optional vector of length q containing the parameters of the covariates effects}
}
\value{
A list with 3 elements: 
\item{SNP1}{the sample matrix of genotypes for the first gene}
\item{SNP2}{the sample matrices of genotypes for the second gene}
\item{Y}{the vector of phenotype values}
\item{Covariates}{the sample matrix of covariates}
}
\description{
\code{SampleGGInter} is used to create a sample of genotype values for two genes and of binary phenotype values. 
The phenotype is generated conditionally to the genes using a logistic model.
}

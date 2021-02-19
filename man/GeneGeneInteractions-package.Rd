\name{GeneGeneInteractions-package}
\alias{GeneGeneInteractions-package}
\alias{GeneGeneInteractions}
\docType{package}
\encoding{UTF-8}
\title{
  Global tests for the detection of gene-gene interaction effects on binary phenotypes
}
\description{
  The GeneGeneInteractions package provides functions to compute global p-values for the detection of 
  gene-gene interaction effects on a binary phenotype. It also provides a function to generate samples
  of genotype and phenotype values under a gene-gene interaction logistic model.
}
\author{
Florian Hebert, Mathieu Emily, David Causeur

Maintainer: Florian Hebert <florian.hebert@agrocampus-ouest.fr>
}
\keyword{Gene-gene interactions}
\seealso{
  \code{\link{GGInterScoreTest}}
  \code{\link{SampleGGInter}}
  \code{\link{OmniGGInterTest}}
  \code{\link{GGInterTest}}
}
\examples{
#Matrices of genotypes for 2 genes, for a population of 100,000 individuals
X1 = matrix(rbinom(10*100000,2,0.4),ncol=10)
X2 = matrix(rbinom(10*100000,2,0.4),ncol=10)

#Sample of genotype and phenotype values. No main effects, interaction effect between:
#SNP 3 of gene 1 and SNP 7 of gene 2 
#SNP 4 of gene 1 and SNP 6 of gene 2


#Dominant-dominant interaction effect
coef = rep(0.2,4)
I = matrix(c(3,7,4,6),ncol=2,byrow=TRUE)
tmp = SampleGGInter(X1,X2,-3,0,1,"A",0,1,"A",
      matrix(rep(coef,2),ncol=4,byrow=TRUE),I,1000,1000)

#Score test statistics (continuous and dummy-coding based models)
test.cont = GGInterScoreTest(tmp$SNP1,tmp$SNP2,tmp$Phenotype,inter.mod="continuous")
test.dummy = GGInterScoreTest(tmp$SNP1,tmp$SNP2,tmp$Phenotype,inter.mod="dummy")

#Global tests
res1 = GGInterTest(tmp$SNP1,tmp$SNP2,tmp$Phenotype,inter.mod="continuous")
res2 = GGInterTest(tmp$SNP1,tmp$SNP2,tmp$Phenotype,inter.mod="dummy")
res3 = GGInterTest(tmp$SNP1,tmp$SNP2,tmp$Phenotype)

#Omnibus tests
res = OmniGGInterTest(tmp$SNP1,tmp$SNP2,tmp$Phenotype)


#Recessive-recessive interaction effect
coef = c(0,0,0,0.5)
I = matrix(c(1,3,2,4,3,7,4,6),ncol=2,byrow=TRUE)
tmp = SampleGGInter(X1,X2,-3,0,1,"A",0,1,"A",
      matrix(rep(coef,nrow(I)),ncol=4,byrow=TRUE),I,1000,1000)

#Score test statistics (continuous and dummy-coding based models)
test.cont = GGInterScoreTest(tmp$SNP1,tmp$SNP2,tmp$Phenotype,inter.mod="continuous")
test.dummy = GGInterScoreTest(tmp$SNP1,tmp$SNP2,tmp$Phenotype,inter.mod="dummy")

#Global tests
res1 = GGInterTest(tmp$SNP1,tmp$SNP2,tmp$Phenotype,inter.mod="continuous")
res2 = GGInterTest(tmp$SNP1,tmp$SNP2,tmp$Phenotype,inter.mod="dummy")
res3 = GGInterTest(tmp$SNP1,tmp$SNP2,tmp$Phenotype)

#Omnibus tests
res = OmniGGInterTest(tmp$SNP1,tmp$SNP2,tmp$Phenotype)
}

SampleGGInter = function(X1,X2,beta0,beta1,I.beta1,mod.beta1=rep("A",length(I.beta1)),
                         beta2,I.beta2,mod.beta2=rep("A",length(I.beta2)),
                         gamma,I.gamma,n0,n1,U=NULL,betaU=NULL){
  Y = PopulationPhenotypeGGInter(X1,X2,beta0,beta1,I.beta1,mod.beta1,
                                 beta2,I.beta2,mod.beta2,gamma,I.gamma,U=U,coef.U=betaU)
  i = which(Y==0)
  ind = sample(c(sample(i,n0,FALSE),sample((1:length(Y))[-i],n1,FALSE)))
  XX1 = X1[ind,,drop=FALSE]
  XX2 = X2[ind,,drop=FALSE]
  Y = matrix(Y[ind],ncol=1)
  U = U[ind,,drop=FALSE]
  return(list(SNP1=XX1,SNP2=XX2,Phenotype=Y,Covariates=U))
}

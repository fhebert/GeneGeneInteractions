ParametricResampling = function(probs,N=1000){
  Y0 = sapply(1:N,function(i){rbinom(length(probs),1,probs)})
  return(Y0)
}

PCsel = function(X,thresh){
  X = scale(X)
  S = cor(X)
  ev = eigen(S)
  ind = which(ev$values>=thresh)
  XX = X%*%ev$vectors[,ind]
  A = matrix(1,ncol=1,nrow=nrow(X))%*%matrix(sqrt(ev$values[ind]),nrow=1)
  XX = XX/A
  return(XX)
}

PopulationPhenotypeGGInter = function(X1,X2,mu,alpha,I.alpha,mod.alpha=rep("A",length(I.alpha)),
                                      beta,I.beta,mod.beta=rep("A",length(I.beta)),gamma,I.gamma,U,coef.U){
  if(!is.matrix(I.gamma)){
    stop("I.gamma must be a matrix with 2 columns")
  }
  if(ncol(I.gamma)!=2){
    stop("I.gamma must be a matrix with 2 columns")
  }
  if(any(I.gamma[,1]>ncol(X1))){
    stop("Indices in the first column of I.gamma must be between 1 and ncol(X1)")
  }
  if(any(I.gamma[,2]>ncol(X2))){
    stop("Indices in the second column of I.gamma must be between 1 and ncol(X2)")
  }
  if(!(is.vector(gamma)|is.matrix(gamma))){
    stop("gamma must either be a vector or a matrix with 4 columns")
  }
  if(is.matrix(gamma)){
    if(ncol(gamma)!=4){
      stop("gamma must either be a vector or a matrix with 4 columns")
    }
  }
  if(!is.null(U)){
    U.eff = U%*%coef.U
  }else{
    U.eff = 0
  }
  if(is.matrix(gamma)){
    if(nrow(I.gamma)!=nrow(gamma)){
      stop("gamma and I.gamma must be matrices with equal number of rows")
    }
    XX = matrix(0,nrow=nrow(X1),ncol=4*nrow(I.gamma))
    ind.XX = 1:4
    for(i in 1:nrow(gamma)){
      ind.tmp = I.gamma[i,]
      x11 = X1[,ind.tmp[1]]==1
      x12 = X1[,ind.tmp[1]]==2
      x21 = X2[,ind.tmp[2]]==1
      x22 = X2[,ind.tmp[2]]==2
      XX[,ind.XX] = cbind(x11*x21,x11*x22,x12*x21,x12*x22)
      ind.XX = ind.XX+4
    }
    XX = as(XX,"dgCMatrix")
    if(any(mod.alpha=="R")){
      indR = I.alpha[which(mod.alpha=="R")]
      X1R = X1[,indR]
      ind0 = which(X1R<2)
      X1R[ind0] = 0
      X1R[-ind0] = 1
      X1[,indR] = X1R
    }
    if(any(mod.alpha=="D")){
      indD = I.alpha[which(mod.alpha=="D")]
      X1D = X1[,indD]
      ind0 = which(X1D<1)
      X1D[ind0] = 0
      X1D[-ind0] = 1
      X1[,indD] = X1D
    }
    if(any(mod.beta=="R")){
      indR = I.beta[which(mod.beta=="R")]
      X2R = X2[,indR]
      ind0 = which(X2R<2)
      X2R[ind0] = 0
      X2R[-ind0] = 1
      X2[,indR] = X2R
    }
    if(any(mod.beta=="D")){
      indD = I.beta[which(mod.beta=="D")]
      X2D = X2[,indD]
      ind0 = which(X2D<1)
      X2D[ind0] = 0
      X2D[-ind0] = 1
      X2[,indD] = X2D
    }
    Y = mu+X1[,I.alpha,drop=FALSE]%*%alpha+X2[,I.beta,drop=FALSE]%*%beta+XX%*%as.vector(t(gamma))+U.eff
  }else{
    if(nrow(I.gamma)!=length(gamma)){
      stop("length of gamma must equal the number of rows of I.gamma")
    }
    XX = matrix(0,nrow=nrow(X1),ncol=nrow(I.gamma))
    for(i in 1:nrow(I.gamma)){
      ind.tmp = I.gamma[i,]
      x1 = X1[,ind.tmp[1]]
      x2 = X2[,ind.tmp[2]]
      XX[,i] = x1*x2
    }
    XX = as(XX,"dgCMatrix")
    if(any(mod.alpha=="R")){
      indR = I.alpha[which(mod.alpha=="R")]
      X1R = X1[,indR]
      ind0 = which(X1R<2)
      X1R[ind0] = 0
      X1R[-ind0] = 1
      X1[,indR] = X1R
    }
    if(any(mod.alpha=="D")){
      indD = I.alpha[which(mod.alpha=="D")]
      X1D = X1[,indD]
      ind0 = which(X1D<1)
      X1D[ind0] = 0
      X1D[-ind0] = 1
      X1[,indD] = X1D
    }
    if(any(mod.beta=="R")){
      indR = I.beta[which(mod.beta=="R")]
      X2R = X2[,indR]
      ind0 = which(X2R<2)
      X2R[ind0] = 0
      X2R[-ind0] = 1
      X2[,indR] = X2R
    }
    if(any(mod.beta=="D")){
      indD = I.beta[which(mod.beta=="D")]
      X2D = X2[,indD]
      ind0 = which(X2D<1)
      X2D[ind0] = 0
      X2D[-ind0] = 1
      X2[,indD] = X2D
    }
    Y = mu+X1[,I.alpha,drop=FALSE]%*%alpha+X2[,I.beta,drop=FALSE]%*%beta+XX%*%gamma+U.eff
  }
  Y = as.vector(1/(1+exp(-Y)))
  Y = rbinom(length(Y),1,Y)
  return(Y)
}


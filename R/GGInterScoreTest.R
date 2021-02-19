eigenK = function(Sigma1,Sigma2,thresh=1){
  ev1 = eigen(Sigma1)
  ev2 = eigen(Sigma2)
  ind = which(ev1$values>thresh)
  ev1$values = ev1$values[ind]
  ev1$vectors = ev1$vectors[,ind,drop=FALSE]
  ind = which(ev2$values>thresh)
  ev2$values = ev2$values[ind]
  ev2$vectors = ev2$vectors[,ind,drop=FALSE]
  lambda = rep(ev1$values,each=length(ev2$values))
  mu = rep(ev2$values,length(ev1$values))
  theta = lambda*mu
  U = ev1$vectors
  V = ev2$vectors
  p1 = ncol(U)
  p2 = ncol(V)
  k = nrow(U)
  U = U[,rep(1:p1,each=p2),drop=FALSE]
  U = U[rep(1:k,each=nrow(V)),,drop=FALSE]
  V = V[,rep(1:p2,p1),drop=FALSE]
  V = V[rep(1:nrow(V),k),,drop=FALSE]
  P = U*V
  ind = order(theta,decreasing=TRUE)
  theta = theta[ind]
  P = P[,ind,drop=FALSE]
  return(list(vectors=P,values=theta))
}

GGInterScoreTest = function(X1,X2,Y,inter.mod="continuous",perm.method="parametric",N=1000,U=NULL){
  n = length(Y)
  if(inter.mod=="continuous"){
    XX1 = matrix(apply(X1,2,function(x){rep(x,ncol(X2))}),nrow=n)
    XX2 = matrix(rep(X2,ncol(X1)),nrow=n)
    S = as(XX1*XX2,"dgCMatrix")
    A1 = X1
    A2 = X2
  }
  if(inter.mod=="dummy"){
    X1 = as(X1,"dgCMatrix")
    X2 = as(X2,"dgCMatrix")
    X1.1 = X1==1
    X1.2 = X1==2
    X2.1 = X2==1
    X2.2 = X2==2
    ind1 = unique(c(which(apply(X1.1,2,sd)==0),which(apply(X1.2,2,sd)==0)))
    ind2 = unique(c(which(apply(X2.1,2,sd)==0),which(apply(X2.2,2,sd)==0)))
    if(length(ind1)>0){
      X1.1 = X1.1[,-ind1]
      X1.2 = X1.2[,-ind1]
    }
    if(length(ind2)>0){
      X2.1 = X2.1[,-ind2]
      X2.2 = X2.2[,-ind2]
    }
    p1 = ncol(X1)-length(ind1)
    p2 = ncol(X2)-length(ind2)
    S1 = as(matrix(apply(cbind(X1.1,X1.2),2,function(x){rep(x,2*p2)}),nrow=n),"dgCMatrix")
    S2 = as(matrix(rep(cbind(X2.1,X2.2),2*p1),nrow=n),"dgCMatrix")
    S = S1*S2
    A1 = as.matrix(cbind(X1.1,X1.2))
    A2 = as.matrix(cbind(X2.1,X2.2))
  }
  W = cbind(1,U,PCsel(A1,0.001*ncol(A1)),PCsel(A2,0.001*ncol(A2)))
  hatY0 = LogitProbs(W,Y)
  Z = as.vector(crossprod(Y-hatY0,S))
  s2Y = sum((Y-hatY0)^2)/length(Y)
  WtS = t(W)%*%(S)
  A = as.matrix(crossprod(S)-t(WtS)%*%solve(crossprod(W))%*%WtS)
  Gamma = s2Y*A
  d = sqrt(diag(Gamma))
  Z = Z/d
  D = matrix(rep(d,ncol(S)),ncol=ncol(S))
  Sigma = Gamma/(D*t(D))
  if(perm.method=="parametric"){
    Y0 = ParametricResampling(hatY0,N)
    hatY0 = PermLogitProbs(W,Y0)
  }else{
    Y0 = sapply(1:N,function(i){sample(Y)})
    hatY0 = PermLogitProbs(W,Y0)
  }
  d0 = Y0-hatY0
  Z0 = as.matrix(crossprod(d0,S))
  sY0 = sqrt(crossprod(rep(1,length(Y)),d0^2)[1,]/length(Y))
  dA = sqrt(diag(A))
  Z0 = Z0/(matrix(sY0,ncol=1)%*%matrix(dA,nrow=1))
  Z0 = as.matrix(Z0)
  diag(Sigma) = 1
  ind = unique(c(which(is.na(Z)),which(colSums(is.na(Z0))>0)))
  if(length(ind)>0){
    Z[ind] = 0
    Z0[,ind] = 0
  }
  ev = eigenK(cor(A1),cor(A2),1)
  return(list(Z=Z,Z0=Z0,Sigma=Sigma,Y0=Y0,inter.mod=inter.mod,X1=A1,X2=A2,
              S=as.matrix(S),W=W,undefined.inter=ind,ev=ev))
}

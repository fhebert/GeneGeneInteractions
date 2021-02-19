L2norm = function(Z,Z0){
  stat = sum(Z^2)
  stat0 = (Z0^2)%*%matrix(1,ncol=1,nrow=ncol(Z0))
  p = mean(stat0>=stat)
  return(list(p=p,TZ=stat,TZ0=stat0))
}

minP = function(Z,Z0){
  stat = max(Z^2)
  stat0 = apply(Z0^2,1,max)
  p = mean(stat0>=stat)
  return(list(p=p,TZ=stat,TZ0=stat0))
}

HC = function(Z,Z0,alpha0=0.5){
  m = length(Z)
  i = 1:floor(alpha0*m)
  pZ = sort(2*(1-pnorm(abs(Z))))[i]
  stat = max(sqrt(m)*(i/m-pZ)/sqrt(pZ*(1-pZ)))
  pZ0 = apply(2*(1-pnorm(abs(Z0))),1,sort)[i,,drop=FALSE]
  stat0 = apply(sqrt(m)*(matrix(rep(i/m,ncol(pZ0)),ncol=ncol(pZ0))-pZ0)/sqrt(pZ0*(1-pZ0)),2,max)
  p = mean(stat0>=stat)
  return(list(p=p,TZ=stat,TZ0=stat0))
}

Hotelling = function(Z,Z0,ev){
  if(length(ev$values)==1){
    A = (1/sqrt(ev$values))*(ev$vectors%*%t(ev$vectors))
  }else{
    A = ev$vectors%*%diag(1/sqrt(ev$values))%*%t(ev$vectors)
  }
  stat = sum((Z%*%A)^2)
  stat0 = rowSums((Z0%*%A)^2)
  p = mean(stat<=stat0)
  return(list(p=p,TZ=stat,TZ0=stat0))
}

simes = function(p,p0){
  m = length(p)
  p = min(m*sort(p)/(1:m))
  p0 = t(apply(p0,1,sort))*m/matrix(rep(1:m,nrow(p0)),ncol=m,byrow=TRUE)
  p0 = apply(p0,1,min)
  p = mean(p0<=p)
  return(p)
}

tsimes = function(p,p0){
  m = length(p)/2
  p = sort(p)[1:m]
  p = min(m*p/(1:m))
  p0 = t(apply(p0,1,sort))[,1:m]*m/matrix(rep(1:m,nrow(p0)),ncol=m,byrow=TRUE)
  p0 = apply(p0,1,min)
  p = mean(p0<=p)
  return(p)
}

fisher = function(p,p0){
  Z = -2*sum(log(p))
  Z0 = -2*rowSums(log(p0))
  p = mean(Z<=Z0)
  return(p)
}

tfisher = function(p,p0){
  m = length(p)/2
  p = sort(p)[1:m]
  p0 = t(apply(p0,1,sort))[,1:m]
  Z = -2*sum(log(p))
  Z0 = -2*rowSums(log(p0))
  p = mean(Z<=Z0)
  return(p)
}

OptT2 = function(vtt,lambda,mgamma) {
  m = length(lambda)
  mlambda = rep(1,nrow(mgamma))%*%t(lambda)
  tmp = sqrt(mgamma/mlambda)
  mtmp = rowSums(tmp)
  tmp = tmp/(mtmp%*%t(rep(1,m)))
  tmp = rowSums(tmp*mlambda*mgamma)
  rowSums(mgamma)%*%t(1/(2*vtt))+tmp%*%t(sign(vtt)*abs(m-(1/(2*vtt))*sum(1/lambda)))
}

MGFR = function(Z,Z0,Sigma=NULL,eigSigma=NULL,vtt=NULL){
  if(is.null(eigSigma)){
    eigSigma = eigen(Sigma)
    ind = which(eigSigma$values>10^(-8))
    lambda = eigSigma$values[ind]
    Omega = t(eigSigma$vectors[,ind,drop=FALSE])
  }else{
    lambda = eigSigma$values
    Omega = t(eigSigma$vectors)
  }
  gamma = ((Omega%*%Z)[,1])^2/lambda
  m = length(Z)
  Z0star = tcrossprod(Z0,Omega)
  mgamma = (Z0star^2)/(rep(1,nrow(Z0))%*%t(lambda))
  mlambda = rep(1,nrow(Z0))%*%t(lambda)
  mchi20 = mgamma
  mlambda = rep(1,nrow(mgamma))%*%t(lambda)
  if(is.null(vtt)){
    vtt = seq(-mean(1/lambda),mean(1/lambda),length=1000)
  }
  vpval = rep(0,length(vtt))
  mT20 = OptT2(vtt,lambda,mchi20)
  indNA = which(rowSums(is.na(mT20))>0)
  if(length(indNA)>0){
    mT20[indNA,] = 0
  }
  vT2 = OptT2(vtt,lambda,matrix(gamma,nrow=1))
  vpval = rep(0,length(vtt))
  mpval0 = matrix(0,nrow=nrow(mgamma),ncol=length(vtt))
  for (j in 1:length(vtt)) {
    vpval[j] = mean(abs(mT20[,j])>abs(vT2[j]))
    mpval0[,j] = 1-(order(order(abs(mT20[,j])))-1)/nrow(mT20)
  }
  minpval0 = apply(mpval0,1,min)
  minpval = min(vpval)
  p = mean(minpval0<=minpval)
  return(list(p=p,p0=minpval0))
}

OmniGGInterTest = function(X1,X2,Y,N=1000,perm.method="parametric",
                           test.method=c("minP","HC","L2norm","Hotelling"),
                           U=NULL){
  res1 = GGInterScoreTest(X1,X2,Y,N=N,perm.method=perm.method,U=U)
  res2 = GGInterScoreTest(X1,X2,Y,inter.mod="dummy",N=N,perm.method=perm.method,U=U)
  tmp = res_continuous = res_dummy = NULL
  if(is.element("minP",test.method)){
    res1_minP = minP(res1$Z,res1$Z0)
    res2_minP = minP(res2$Z,res2$Z0)
    tmp = cbind(tmp,order(order(res1_minP$TZ0)),order(order(res2_minP$TZ0)))
    pp = res1_minP$p
    names(pp) = "minP"
    res_continuous = c(res_continuous,pp)
    pp = res2_minP$p
    names(pp) = "minP"
    res_dummy = c(res_dummy,pp)
  }
  if(is.element("HC",test.method)){
    res1_HC = HC(res1$Z,res1$Z0)
    res2_HC = HC(res2$Z,res2$Z0)
    tmp = cbind(tmp,order(order(res1_HC$TZ0)),order(order(res2_HC$TZ0)))
    pp = res1_HC$p
    names(pp) = "HC"
    res_continuous = c(res_continuous,pp)
    pp = res2_HC$p
    names(pp) = "HC"
    res_dummy = c(res_dummy,pp)
  }
  if(is.element("L2norm",test.method)){
    res1_L2 = L2norm(res1$Z,res1$Z0)
    res2_L2 = L2norm(res2$Z,res2$Z0)
    tmp = cbind(tmp,order(order(res1_L2$TZ0)),order(order(res2_L2$TZ0)))
    pp = res1_L2$p
    names(pp) = "L2norm"
    res_continuous = c(res_continuous,pp)
    pp = res2_L2$p
    names(pp) = "L2norm"
    res_dummy = c(res_dummy,pp)
  }
  ev1 = res1$ev
  ev2 = res2$ev
  if(is.element("Hotelling",test.method)){
    res1_Hotelling = Hotelling(res1$Z,res1$Z0,ev1)
    res2_Hotelling = Hotelling(res2$Z,res2$Z0,ev2)
    tmp = cbind(tmp,order(order(res1_Hotelling$TZ0)),order(order(res2_Hotelling$TZ0)))
    pp = res1_Hotelling$p
    names(pp) = "Hotelling"
    res_continuous = c(res_continuous,pp)
    pp = res2_Hotelling$p
    names(pp) = "Hotelling"
    res_dummy = c(res_dummy,pp)
  }
  p0 = 1-tmp/N
  res = simes(c(res_dummy,res_continuous),p0)
  
  res1_MGFR = MGFR(res1$Z,res1$Z0,eigSigma=ev1)
  res2_MGFR = MGFR(res2$Z,res2$Z0,eigSigma=ev2)
  p0_MGFR = cbind(order(order(res1_MGFR$p0))/N,order(order(res2_MGFR$p0))/N)
  p_MGFR = c(res1_MGFR$p,res2_MGFR$p)
  res_MGFR = simes(p_MGFR,p0_MGFR)
  names(p_MGFR) = rep("MGFR",2)
  res_continuous = c(res_continuous,p_MGFR[1])
  res_dummy = c(res_dummy,p_MGFR[2])
  p0 = cbind(p0,p0_MGFR)
  p0_continuous = p0[,seq(1,ncol(p0),by=2)]
  p0_dummy = p0[,-seq(1,ncol(p0),by=2)]
  colnames(p0_continuous) = colnames(p0_dummy) = names(res_continuous)
  
  res = list(res_omni=res,res_omni_MGFR=res_MGFR,res_continuous=res_continuous,res_dummy=res_dummy,
             p0_continuous=p0_continuous,p0_dummy=p0_dummy,
             fit_continuous=res1,fit_dummy=res2)
  return(res)
}
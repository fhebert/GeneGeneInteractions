GGInterTest = function(X1,X2,Y,N=1000,perm.method="parametric",inter.mod=c("continuous","dummy"),U=NULL){
  if(!any(inter.mod%in%c("continuous","dummy"))){
    stop("Possible values for inter.mod are 'continuous', 'dummy' and c('continuous','dummy')")
  }
  res_continuous = res_dummy = NULL
  if(any(inter.mod=="continuous")){
    res1 = GGInterScoreTest(X1,X2,Y,N=N,perm.method=perm.method,U=U)
    res1_minP = minP(res1$Z,res1$Z0)$p
    res1_HC = HC(res1$Z,res1$Z0)$p
    res1_L2 = L2norm(res1$Z,res1$Z0)$p
    ev1 = res1$ev
    res1_Hotelling = Hotelling(res1$Z,res1$Z0,ev1)$p
    res1_MGFR = MGFR(res1$Z,res1$Z0,eigSigma=ev1)$p
    res_continuous = c(res1_minP,res1_HC,res1_L2,res1_Hotelling,res1_MGFR)
    names(res_continuous) = c("minP","HC","L2","Hotelling","MGFR")
  }
  if(any(inter.mod=="dummy")){
    res2 = GGInterScoreTest(X1,X2,Y,inter.mod="dummy",N=N,perm.method=perm.method,U=U)
    res2_minP = minP(res2$Z,res2$Z0)$p
    res2_HC = HC(res2$Z,res2$Z0)$p
    res2_L2 = L2norm(res2$Z,res2$Z0)$p
    ev2 = res2$ev
    res2_Hotelling = Hotelling(res2$Z,res2$Z0,ev2)$p
    res2_MGFR = MGFR(res2$Z,res2$Z0,eigSigma=ev2)$p
    res_dummy = c(res2_minP,res2_HC,res2_L2,res2_Hotelling,res2_MGFR)
    names(res_dummy) = c("minP","HC","L2","Hotelling","MGFR")
  }
  res = list(res_continuous=res_continuous,res_dummy=res_dummy)
  return(res)
}

//[[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;


// log-likelihood of a logistic model
double logLnbeta(Eigen::MatrixXd X, Eigen::VectorXd Y, Eigen::VectorXd beta){
  Eigen::MatrixXd Xbeta = X*beta;
  double res = (Y.array()*Xbeta.array()-(1+Xbeta.array().exp()).log()).sum();
  return(res);
}

// Efficient cross-product X'X
Eigen::MatrixXd CrossProd(Eigen::MatrixXd X){
  int n = X.cols();
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(n,n);
  res.selfadjointView<Eigen::Lower>().rankUpdate(X.transpose());
  res.triangularView<Eigen::StrictlyUpper>() = res.transpose(); 
  return (res);
}

// Efficient weighted cross-product X'diag(w)X
Eigen::MatrixXd WeightedCrossProd(Eigen::MatrixXd X, Eigen::VectorXd w){
  int n = X.cols();
  w = w.cwiseSqrt();
  Eigen::MatrixXd WX = X.array().colwise()*w.array();
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(n,n);
  res.selfadjointView<Eigen::Lower>().rankUpdate(WX.transpose());
  res.triangularView<Eigen::StrictlyUpper>() = res.transpose(); 
  return (res);
}

// Estimated probabilities of a logistic model
//[[Rcpp::export]]
Eigen::VectorXd LogitProbs(Eigen::MatrixXd X, Eigen::VectorXd Y){
  int p = X.cols();
  Eigen::VectorXd beta_init = Eigen::VectorXd(p);
  Eigen::MatrixXd Xt = X.transpose();
  Eigen::MatrixXd A = Eigen::MatrixXd::Zero(p,p);
  beta_init.fill(0);
  Eigen::VectorXd pbeta = (1+(-X*beta_init).array().exp()).inverse();
  A = (WeightedCrossProd(X,pbeta.array()*(1-pbeta.array()))).inverse();
  Eigen::VectorXd new_beta = beta_init+A*(Xt*(Y-pbeta));
  double d0 = -logLnbeta(X,Y,beta_init);
  double d1 = -logLnbeta(X,Y,new_beta);
  double d = std::abs(d1-d0)/(std::abs(d1)+0.1);
  int niter = 1;
  while((d>0.00000001)&(niter<=50)){
    beta_init = new_beta;
    pbeta = (1+(-X*beta_init).array().exp()).inverse();
    A = (WeightedCrossProd(X,pbeta.array()*(1-pbeta.array()))).inverse();
    new_beta.noalias() += A*(Xt*(Y-pbeta));
    d0 = -logLnbeta(X,Y,beta_init);
    d1 = -logLnbeta(X,Y,new_beta);
    d = std::abs(d1-d0)/(std::abs(d1)+0.1);
    niter++;
  }
  return(pbeta);
}

// Probabilities of a logistic model for each column of Y0
//[[Rcpp::export]]
Eigen::MatrixXd PermLogitProbs(Eigen::MatrixXd X, Eigen::MatrixXd Y0){
  int n = X.rows();
  int N = Y0.cols();
  Eigen::MatrixXd res = Eigen::MatrixXd::Zero(n,N);
  Eigen::VectorXd Y = Eigen::VectorXd(n);
  for(int i=0;i<N;i++){
    Y = Y0.col(i);
    res.col(i) = LogitProbs(X,Y);
  }
  return(res);
}

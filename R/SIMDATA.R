

library(mvtnorm)
library(Matrix)


#### set 1:  
gen_Beta <- function(p, pact, set,  beta_val=NULL){
  beta <- rep(0, p)
  if (set==1) beta[1:pact] <- c(3, 1.5, 1, 0.5)
  if (set==2) beta[1:pact] <- c(3, 1.5, 1, 0.5, -3, -1.5, -1, -0.5) 
  if (set==3) beta[1:pact] = beta_val
  return(beta)
}


gen_X <- function(n, p, pact, case, rho1,  rho12, rho2, v){
  ar1 <- function(n, rho) rho^(abs(matrix(1:n-1, nrow=n, ncol=n, byrow=TRUE)-(1:n-1)))
  if (case == 1) sigm <- diag(p)
  if (case == 2){
    sigm <- matrix(rho12, p, p)
    sigm[1:(pact), 1:(pact)] <- ar1(pact, rho1)
    sigm[(pact+1):(p), (pact+1):(p)] <- ar1((p-pact), rho2)
    sigm=  as.matrix(nearPD(sigm)$mat)
  }
  X <- mvnfast::rmvn(n, mu = rep(0, p), sigma = sigm) 
  return(X)
}



simdata <- function(n, p, pact, v, case, set,  rho1,  rho12, rho2, beta_val, link){
  X0 = as.matrix(apply(gen_X(2*n, p, pact, case,  rho1,  rho12, rho2, v=v),2 , scale))
  Bc = gen_Beta(p, pact, set,  beta_val)
  eta = X0%*%Bc
  Z = numeric()
  E0 = numeric()
  for(i in 1:length(eta)){
    if (link =="probit"){
      Z[i] = pnorm(eta[i])
    }
    if (link =="tlink"){
      Z[i] = pt(eta[i], v)
    }
    if (link =="logit"){
      Z[i] = plogis(eta[i])
    }
    
    E0[i] = rbinom(1, 1, Z[i]) 
  }
  
  train_id <- sample(seq_len(length(E0)), size = ceiling(0.5 * length(E0)) )
  Xtrain = X0[train_id, ]; Xtest = X0[-train_id, ] 
  Etrain = E0[train_id]; Etest = E0[-train_id]; Ztrain = Z[train_id]
  
  return(list(Xtrain=as.matrix(Xtrain), Xtest=as.matrix(Xtest),
              Etrain=Etrain,Etest =Etest, Ztrain= Ztrain, B=Bc, p=p, pact= pact,
              case = case, set = set, link = link))
}


library(mvtnorm)
library(truncnorm)
library(Matrix)
library(invgamma)
library(mvnfast)

Maj_sp=function(betahat,Zgamma,cutoff){
  Z_spr = ifelse(Zgamma>=cutoff, 1, 0)
  beta_spr = ifelse(Z_spr == 1,betahat,0)
  return(list("beta_spr"= beta_spr, "Z_spr" = Z_spr))
}
################################################################################
### Define Exact Gibbs Function
################################################################################
# Xia edit: remove Z and set initialization in the codes; update tau1_ with r, s hyperparameter
exactGibbs <- function(X, Y, gamma, tau0_2, r, s, q, nburn=2000, niter=2000){
  #####
  # X is a nxp matrix
  # Y is a 1xn vector with 0's or 1's
  # Z is a 1xn vector
  # gamma is a 1xp vector with 0's or 1's
  # tau0 is a real number 
  # tau1 is a real number
  # q is a real number
  #####
  n <- nrow(X)         
  p <- ncol(X)
  k = sum(gamma)
# Xia edit: initialize Z
  Z = rep(0, n)
  for(i in 1:n){
        if (Y[i]==1) Z[i] <- truncnorm::rtruncnorm(1, a=0, mean=0, sd= 1)
        if (Y[i]==0) Z[i] <- truncnorm::rtruncnorm(1, b=0, mean=0, sd= 1)}

  # Xia edit: preset the matrix to hold the MCMC results; remove out Z; add Sizeout
  outbeta<- matrix(NA, nrow = niter, ncol = p)  
  gammamat<- matrix(NA, nrow = niter, ncol = p)
  outtau1 <- c()
  outgamma <- rep(0, p)
  Sizeout<- numeric()
  
  inv_tau0_2 <- 1/tau0_2
# Xia edit: add the prior distribuiton on the slab variance as in T
  tau1_2 = invgamma::rinvgamma(1, shape = r, rate = s) 
  inv_tau1_2 <- 1/tau1_2
  # xia edit then the following is not a constant
####  const <- (q*sqrt(tau0_2))/((1-q)*sqrt(tau1_2))
  
  # start updating
  for (itr in 1:(nburn+niter)){
    #print(paste("iteration", itr))
    # Check Processing
    if (itr %% 2000 == 0) print(paste("finish iteration", itr))
    
    ##### update beta
    D_gamma <- diag(gamma*inv_tau1_2+(1-gamma)*inv_tau0_2) 
    V =   (t(X)%*%X)+D_gamma
    sigma1 <- solve(V + diag(rep(.Machine$double.eps, p)))
    #sigma1 =   as.matrix(nearPD(sigma1)$mat) 
    mean1 <- sigma1%*%t(X)%*%Z
    beta <- t(unlist(mvnfast::rmvn(1, mu = mean1, sigma = sigma1))) 
    
    ##### update gamma
    for (j in 1:p){
      temp1 <- beta[j]^2*0.5*(inv_tau0_2-inv_tau1_2)
#xia edit:  update it within MCMC
      const <- (q*sqrt(tau0_2))/((1-q)*sqrt(tau1_2))
      log_d_j <- log(const)+temp1
      if (exp(log_d_j)==Inf) gamma[j] <- 1
      else gamma[j] <- rbinom(1, 1, exp(log_d_j)/(1+exp(log_d_j)))
    }

       # Xia edit: add the sampling on tau1_2 as in T
    ##################### Update tau_1^2 #######################
    tau1_2 = rinvgamma(1, shape = k/2 + r, rate = s + sum(gamma*beta*beta)/2) 
    inv_tau1_2 = 1/tau1_2
    # xia add: the size of the model
    k = sum(gamma)
    #### update Z
    
    for (i in 1:n){
      if (Y[i]==1) Z[i] <- rtruncnorm(1, a=0, mean=X[i,]%*%beta, sd=1)
      if (Y[i]==0) Z[i] <- rtruncnorm(1, b=0, mean=X[i,]%*%beta, sd=1)
    }
    
    if (itr > nburn) {
      outtau1[(itr-nburn)]<- tau1_2
      outbeta[(itr-nburn),] <- beta
      gammamat[(itr-nburn), ] = unlist(gamma)
      outgamma <- outgamma+gamma
      Sizeout[itr-nburn] <- k
    }
  }
  
  est = Maj_sp(betahat= colSums(outbeta)/niter,Zgamma= outgamma/niter,cutoff = 0.5)
  return(list("betahat"=est$beta_spr, "gammahat" = est$Z_spr,  "outgamma"=outgamma/niter, 
               "Size" = Sizeout, "betamcmc"=outbeta, "tau1" = outtau1,"gammamat"=gammamat)) 

}



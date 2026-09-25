library(mvtnorm)
library(truncnorm)
library(invgamma)
library(mvnfast)

Maj_sp=function(betahat,Zgamma,cutoff){
  Z_spr = ifelse(Zgamma>=cutoff, 1, 0)
  beta_spr = ifelse(Z_spr == 1,betahat,0)
  return(list("beta_spr"= beta_spr, "Z_spr" = Z_spr))
}
################################################################################
### Define skinny Gibbs Function
################################################################################

# Xia edit: remove Z and set initialization in the codes; update tau1_ with r, s hyperparameter
###skinnyGibbs <- function(X, Y, gamma, tau0_2, tau1_2, q, nburn=2000, niter=2000){
skinnyGibbs <- function(X, Y, gamma, tau0_2, r, s, q, nburn=2000, niter=2000){
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
###  outZ <- rep(0, n)
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
    idx_active <- which(gamma==1)
    idx_inactive <- which(gamma==0)
    n_active <- length(idx_active)
    
    beta <- rep(0, p)
    if (n_active==0) beta <- rnorm(n=p, mean=0, sd=1/sqrt(n-1+inv_tau0_2)) 
    else {
      # Generate active beta
      X_active <- matrix(X[, idx_active], nrow=n)
      sigma1 <- solve(t(X_active)%*%X_active+inv_tau1_2*diag(n_active))
      mean1 <- sigma1%*%t(X_active)%*%Z
      beta[idx_active] <- rmvnorm(n=1, mean=mean1, sigma=sigma1)
      # Generate inactive beta 
      beta[idx_inactive] <- rnorm(n=p-n_active, mean=0, sd=1/sqrt(n-1+inv_tau0_2))
    }      
    
    ##### update gamma
    if (n_active==0) temp2_2_2 <- Z
    else{
      X_active <- matrix(X[, idx_active], nrow=n) # guarantee it is matrix form. 
      beta_active <- beta[idx_active]
      temp2_2_2 <- Z-X_active%*%beta_active # dim(temp2_2_2) is nx1
    }
    for (j in 1:p){
      temp1 <- beta[j]^2*0.5*(inv_tau0_2-inv_tau1_2)
      temp2_1 <- beta[j]*X[,j]
      if (j %in% idx_active){
        temp2_2_1 <- Z-matrix(X_active[, -which(idx_active==j)], 
                              nrow=n)%*%beta_active[-which(idx_active==j)]
        temp2 <- temp2_1%*%temp2_2_1
      }
      else temp2 <- temp2_1%*%temp2_2_2
# Xia edit: replace log const      
###      log_d_j_star <- log(const)+temp1+temp2
      log_d_j_star <- log(q*sqrt(tau0_2)) - log((1 - q)*sqrt(tau1_2)) + temp1 + temp2
      if (exp(log_d_j_star)==Inf) gamma[j] <- 1
      else gamma[j] <- rbinom(1, 1, exp(log_d_j_star)/(1+exp(log_d_j_star)))
    }

    # Xia edit: add the sampling on tau1_2 as in T
    ##################### Update tau_1^2 #######################
    tau1_2 = rinvgamma(1, shape = n_active/2 + r, rate = s + sum(gamma*beta*beta)/2) 
    inv_tau1_2 = 1/tau1_2
    #### update Z
      # Xia: add the test if there is an empty model
    k = sum(gamma)
    if(k == 0) {
      for(i in 1:n){
        if (Y[i]==1) Z[i] <- truncnorm::rtruncnorm(1, a=0, mean=0, sd= 1)
        if (Y[i]==0) Z[i] <- truncnorm::rtruncnorm(1, b=0, mean=0, sd= 1)
                   }
       }else{
     idx_active <- which(gamma==1)
     X_active <- matrix(X[, idx_active], nrow=n) # guarantee it is matrix form. 
     beta_active <- beta[idx_active]
      for (i in 1:n){
      if (Y[i]==1) Z[i] <- rtruncnorm(1, a=0, mean=X_active[i,]%*%beta_active, sd=1)
      if (Y[i]==0) Z[i] <- rtruncnorm(1, b=0, mean=X_active[i,]%*%beta_active, sd=1)
                    }
            }
    
    if (itr > nburn) {
#xia edit: change outbeta as a matrix holding MCMC output; remove the output of Z we do not need it.
      # also add sizeout to be consistent with t.
      outtau1[(itr-nburn)]<- tau1_2
      outbeta[(itr-nburn),] <- beta
      gammamat[(itr-nburn), ] = unlist(gamma)
      outgamma <- outgamma+gamma
#####      outZ <- outZ+Z
      Sizeout[itr-nburn] <- k
    }
  }

  #xia edit: make outbeta to colSums  
###  est = Maj_sp(betahat= outbeta/niter,Zgamma= outgamma/niter,cutoff = 0.5)
### return(list("betahat"=est$beta_spr, "gammahat" = est$Z_spr,"outbeta"=outbeta , "outgamma"=outgamma, "outZ"=outZ/niter))
  est = Maj_sp(betahat= colSums(outbeta)/niter,Zgamma= outgamma/niter,cutoff = 0.5)
  return(list("betahat"=est$beta_spr, "gammahat" = est$Z_spr,  "outgamma"=outgamma/niter, 
               "Size" = Sizeout, "betamcmc"=outbeta, "tau1" = outtau1,"gammamat"=gammamat)) 

}




# rm(list=ls())

library(MASS)
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
SGlogit <- function(X, E, gammaZ , b0, r, s, qr, tau1_2, nburn=2000, niter=4000){
  #####
  # X (Design Matrix):      n*p matrix
  # E (Binary Response):    1*n vector with 0's or 1's
  # gammaZ (Variable Indicator): 1*p vector with 0's or 1's
  # tau0_2 (Inactive Var):  a real value 
  # tau1_2 (Active Var):    a real value
  # qr:                      a real value
  #####
  
  n <- nrow(X)
  p <- ncol(X)
  
  tau0_2= b0/n
  inv_tau0_2 = 1/tau0_2
  tau1_2 = invgamma::rinvgamma(1, shape = r, rate = s) 
  inv_tau1_2 = 1/tau1_2
  v <- 7.3
  w_2 <- (pi^2)*(v-2)/(3*v)
  Omega_sq = invgamma::rinvgamma(n, shape = v/2, rate = w_2*v/2 )
  k = sum(gammaZ)
  beta = rep(0,p)
  Z = rep(0, n)
  for (i in 1:n){
    if (E[i]==1) Z[i] <- rtruncnorm(n=1, a=0, mean=0, sd=sqrt(Omega_sq[i]))
    if (E[i]==0) Z[i] <- rtruncnorm(n=1, b=0, mean=0, sd=sqrt(Omega_sq[i]))}
  

  outbeta<- matrix(NA, nrow = niter, ncol = p)  
  gammamat<- matrix(NA, nrow = niter, ncol = p)
  outtau1 <- c()
  outgammaZ <- rep(0, p)
  Sizeout<- numeric()
  
  # start updating
  for (itr in 1:(nburn+niter)){
    
    if (itr %% 2000 == 0) print(paste("iteration", itr))
    
    idx_active <- which(gammaZ==1)
    idx_inactive <- which(gammaZ==0)
  
    ##### update beta
    beta <- rep(0, p)
    if (k == 0){
      beta <- rnorm(p, mean = 0, sd = 1 / sqrt(n-1 + inv_tau0_2))
    }else{
      X_active <- as.matrix(X[, idx_active])
      sigma1 <-  solve(crossprod(X_active, diag(1/Omega_sq))%*%X_active +  inv_tau1_2*diag(k))
      mean1 <- sigma1 %*% crossprod(X_active, diag(1/Omega_sq))%*%Z
      beta[idx_active] <- t(unlist(rmvn(1, mu = mean1, sigma = sigma1)))
      beta[idx_inactive] <- rnorm(n = p-k, mean = 0, sd = 1/sqrt(n - 1 + inv_tau0_2) )
      
    }
    
    
    ##### update Z latent variable amd omegasq
    if(k == 0) {
      for(i in 1:n){
        if (E[i]==1) Z[i] <- truncnorm::rtruncnorm(1, a=0, mean=0, sd= sqrt(Omega_sq[i]))
        if (E[i]==0) Z[i] <- truncnorm::rtruncnorm(1, b=0, mean=0, sd= sqrt(Omega_sq[i]))
        #######################. UPDATE Omega_sq ######################
        Omega_sq[i] =  invgamma::rinvgamma(1, shape = (v+1)/2, rate = (w_2*v + (Z[i])^2)/2 )
      }
    }else{
      for(i in 1:n){
        if (E[i]==1) Z[i] <- truncnorm::rtruncnorm(1, a=0, mean=sum(X_active[i,]*beta[idx_active]), sd= sqrt(Omega_sq[i]))
        if (E[i]==0) Z[i] <- truncnorm::rtruncnorm(1, b=0, mean=sum(X_active[i,]*beta[idx_active]), sd= sqrt(Omega_sq[i]))
        #######################. UPDATE Omega_sq ######################
        Omega_sq[i] =  invgamma::rinvgamma(1, shape = (v+1)/2, rate = (w_2*v + (Z[i] - sum(X_active[i,]*beta[idx_active]))^2)/2 )
      }
    }
    
    ##################### Update tau_1^2 #######################
    tau1_2 = rinvgamma(1, shape = k/2 + r, rate = s + sum(gammaZ*beta*beta)/2) 
    inv_tau1_2 = 1/tau1_2
    ##### update gammaZ (Variable Indicator)
    
    if(k > 0) { temp_2_12 <- Z - as.matrix(X[, idx_active]) %*% beta[idx_active]
    }else{ temp_2_12 <- Z - numeric(n) }
    
    for(j in seq_len(p)) {
      beta_active = beta[idx_active]
      temp_0 <- 0.5 * beta[j]^2 * (inv_tau0_2 - inv_tau1_2)
      temp_2_1 <- beta[j]*X[,j]*(1/Omega_sq)
      
      if(j %in% idx_active){
        temp_2_11 <- (Z- as.matrix(X_active[, -which(idx_active==j)])%*%beta_active[-which(idx_active==j)])
        temp_2 <- sum(temp_2_1*temp_2_11)
      }else{
        temp_2 <- sum(temp_2_1*temp_2_12)
      }
      temp_3 <- 0.5*(beta[j])^2*sum(X[,j]*(1 - (1/Omega_sq))*X[,j])
      
      log_d_j <- log(qr*sqrt(tau0_2)) - log((1 - qr)*sqrt(tau1_2)) + temp_0 + temp_2 + temp_3
      
      if (is.infinite(exp(log_d_j)) ==T) {
        gammaZ[j] <- 1
      } else { 
        gammaZ[j] <- rbinom(1, 1, exp(log_d_j)/(1 + exp(log_d_j)))
      }
    }
    
    
    k = sum(gammaZ)
    
    
    ###Results 
    if(itr > nburn) {
      outtau1[(itr-nburn)]<- tau1_2
      gammamat[(itr-nburn), ] = unlist(gammaZ)
      outbeta[(itr-nburn), ] = unlist(beta)
      outgammaZ <- outgammaZ+gammaZ
      Sizeout[itr-nburn] <- k
    }
  }
  
  est = Maj_sp(betahat= colSums(outbeta)/niter,Zgamma= outgammaZ/niter,cutoff = 0.5)
  
  return(list("betahat"=est$beta_spr, "gammaZhat"  = est$Z_spr,  "gammamat"=gammamat, 
              "Size" = Sizeout,  "betamcmc" = outbeta,"tau1" = outtau1,
              "outgammaZ"=outgammaZ/niter))
}

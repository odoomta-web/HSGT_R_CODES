

library(MASS)
library(invgamma)
library(mvnfast)


Maj_sp=function(betahat,Zgamma,cutoff){
  Z_spr = ifelse(Zgamma>=cutoff, 1, 0)
  beta_spr = ifelse(Z_spr == 1,betahat,0)
  return(list("beta_spr"= beta_spr, "Z_spr" = Z_spr))
}

# Xia edit: remove Z
SkinyGibbsT = function(X, E, gammaZ, b0, v, r, s, qr, nburn=2000, niter=4000 ){
  
  n <- nrow(X)
  p <- ncol(X)
  tau0_sq= b0/n
  inv_tau0sq = 1/tau0_sq 
  tau1_sq = invgamma::rinvgamma(1, shape = r, rate = s) 
  inv_tau1sq = 1/tau1_sq
  Omega_sq = invgamma::rinvgamma(n, shape = v/2, rate = v/2 )
  k = sum(gammaZ)
  beta = rep(0,p)
# Xia edit: initialize Z
  Z = rep(0, n)
  for(i in 1:n){
        if (E[i]==1) Z[i] <- truncnorm::rtruncnorm(1, a=0, mean=0, sd= sqrt(Omega_sq[i]))
        if (E[i]==0) Z[i] <- truncnorm::rtruncnorm(1, b=0, mean=0, sd= sqrt(Omega_sq[i]))}

  # Xia edit: preset the matrix to hold the MCMC results  
 # outbeta <- numeric()
  outbeta<- matrix(NA, nrow = niter, ncol = p)  
  gammamat<- matrix(NA, nrow = niter, ncol = p)
  outtau1 <- c()
  outgammaZ <- rep(0, p)
  Sizeout<- numeric()
  
  for(itr in 1:(nburn+niter)){
    
    if(itr %% 2000 == 0)print(paste("finished iteration", itr))
    
    idx_active <- which(gammaZ == 1)
    idx_inactive <- which(gammaZ == 0)
    
    ################### UPDATE Beta
    if (k == 0){
      beta <- rnorm(p, mean = 0, sd = 1 / sqrt(n-1 + inv_tau0sq))
    }else{
      X_active <- as.matrix(X[, idx_active])
      #t(X_active)%*%diag(1/Omega_sq)%*%X_active 
      sigma1 <-  solve(crossprod(X_active, diag(1/Omega_sq))%*%X_active +  inv_tau1sq*diag(k))
      mean1 <- sigma1 %*% crossprod(X_active, diag(1/Omega_sq))%*%Z
      beta[idx_active] <- t(unlist(rmvn(1, mu = mean1, sigma = sigma1)))
      beta[idx_inactive] <- rnorm(n = p-k, mean = 0, sd = 1/sqrt(n - 1 + inv_tau0sq) )
      
    }
    
    ##################### Update tau_1^2 #######################
    tau1_sq = rinvgamma(1, shape = k/2 + r, rate = s + sum(gammaZ*beta*beta)/2) 
    inv_tau1sq = 1/tau1_sq
    #################### UPDATE Z and UPDATE Omega_sq ################
    #Xia edits
  
    if(k == 0) {
###      Z = Z
      for(i in 1:n){
        if (E[i]==1) Z[i] <- truncnorm::rtruncnorm(1, a=0, mean=0, sd= sqrt(Omega_sq[i]))
        if (E[i]==0) Z[i] <- truncnorm::rtruncnorm(1, b=0, mean=0, sd= sqrt(Omega_sq[i]))
        #######################. UPDATE Omega_sq ######################
        Omega_sq[i] = rinvgamma(1, shape = (v+1)/2, rate = (v + (Z[i])^2)/2 )
      }
    }else{
      for(i in 1:n){
        if (E[i]==1) Z[i] <- truncnorm::rtruncnorm(1, a=0, mean=sum(X_active[i,]*beta[idx_active]), sd= sqrt(Omega_sq[i]))
        if (E[i]==0) Z[i] <- truncnorm::rtruncnorm(1, b=0, mean=sum(X_active[i,]*beta[idx_active]), sd= sqrt(Omega_sq[i]))
        #######################. UPDATE Omega_sq ######################
        Omega_sq[i] = rinvgamma(1, shape = (v+1)/2, rate = (v + (Z[i] - sum(X_active[i,]*beta[idx_active]))^2)/2 )
      }
    }
    
    ###Update GammaZ
    if(k > 0) { temp_2_12 <- Z - as.matrix(X[, idx_active]) %*% beta[idx_active]
          }else{ temp_2_12 <- Z - numeric(n) }
    
    for(j in seq_len(p)) {
      beta_active = beta[idx_active]
      temp_0 <- 0.5 * beta[j]^2 * (inv_tau0sq - inv_tau1sq)
      temp_2_1 <- beta[j]*X[,j]*(1/Omega_sq)
      
      if(j %in% idx_active){
        temp_2_11 <- (Z- as.matrix(X_active[, -which(idx_active==j)])%*%beta_active[-which(idx_active==j)])
        temp_2 <- sum(temp_2_1*temp_2_11)
      }else{
        temp_2 <- sum(temp_2_1*temp_2_12)
      }
      temp_3 <- 0.5*(beta[j])^2*sum(X[,j]*(1 - (1/Omega_sq))*X[,j])
      
      log_d_j <- log(qr*sqrt(tau0_sq)) - log((1 - qr)*sqrt(tau1_sq)) + temp_0 + temp_2 + temp_3
      
      if (is.infinite(exp(log_d_j)) ==T) {
        gammaZ[j] <- 1
      } else { 
        gammaZ[j] <- rbinom(1, 1, exp(log_d_j)/(1 + exp(log_d_j)))
      }
    }
  
 
    k = sum(gammaZ)
  
    if (itr > nburn) {
      outtau1[(itr-nburn)]<- tau1_sq
      gammamat[(itr-nburn), ] = unlist(gammaZ)
      outbeta[(itr-nburn), ] = unlist(beta)
      outgammaZ <- outgammaZ+gammaZ
      Sizeout[itr-nburn] <- k
    }
    
  }
  
  est = Maj_sp(betahat= colSums(outbeta)/niter,Zgamma= outgammaZ/niter,cutoff = 0.5)
  
  return(list("betahat"=est$beta_spr, "gammaZhat"  = est$Z_spr,  "outgammaZ"=outgammaZ/niter, 
              "Size" = Sizeout,  "betamcmc" = outbeta,  "tau1" = outtau1,"gammamat"=gammamat))
  
}


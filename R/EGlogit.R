library(MASS)
library(invgamma)
library(mvnfast)


Maj_sp=function(betahat,Zgamma,cutoff){
  Z_spr = ifelse(Zgamma>=cutoff, 1, 0)
  beta_spr = ifelse(Z_spr == 1,betahat,0)
  return(list("beta_spr"= beta_spr, "Z_spr" = Z_spr))
}

EGlogit = function(X, E, gammaZ, b0, r, s, qr, nburn=2000, niter=4000){
  
  n <- nrow(X);p <- ncol(X)
  tau0_sq= b0/n
  inv_tau0sq = 1/tau0_sq 
  tau1_sq = invgamma::rinvgamma(1, shape = r, rate = s) 
  inv_tau1sq = 1/tau1_sq
  v <- 7.3
  w_2 <- (pi^2)*(v-2)/(3*v)
  Omega_sq = invgamma::rinvgamma(n, shape = v/2, rate = w_2*v/2 )
  k = sum(gammaZ)
  eps <- 100*.Machine$double.eps * diag(p)
  Z = rep(0, n)
  for(i in 1:n){
    if (E[i]==1) Z[i] <- truncnorm::rtruncnorm(1, a=0, mean=0, sd= sqrt(Omega_sq[i]))
    if (E[i]==0) Z[i] <- truncnorm::rtruncnorm(1, b=0, mean=0, sd= sqrt(Omega_sq[i]))}
  
  outbeta<- matrix(NA, nrow = niter, ncol = p)  
  gammamat<- matrix(NA, nrow = niter, ncol = p)
  outtau1 <- c()
  outgammaZ <- rep(0, p)
  Sizeout<- numeric()
  
  for(itr in 1:(nburn+niter)){
    
    if(itr %% 2000 == 0)print(paste("iteration ", itr))
    
    ### UPDATE Beta
    D = diag(inv_tau0sq*(1-gammaZ) + inv_tau1sq*gammaZ)
    sigma1 = solve( t(X)%*%diag(1/Omega_sq)%*%X+ D + eps*diag(p))
    mean1 <- sigma1%*%t(X)%*%diag(1/Omega_sq)%*%Z
    beta <- t(unlist(mvnfast::rmvn(1, mu = mean1, sigma = sigma1)))
    
    ################ UPDATE Z and Omega_sq ##########################
    for(i in 1:n){
      ### UPDATE Z
      if (E[i]==1) Z[i] <- truncnorm::rtruncnorm(1, a=0, mean=sum(X[i,]*beta), sd= sqrt(Omega_sq[i]))
      if (E[i]==0) Z[i] <- truncnorm::rtruncnorm(1, b=0, mean=sum(X[i,]*beta), sd=sqrt(Omega_sq[i]))
      ### UPDATE Omega_sq
      Omega_sq[i] = rinvgamma(1, shape = (v+1)/2, rate = (w_2*v + (Z[i] - sum(X[i,]*beta))^2 )/2)
    }
    
    #################. Update tau_1^2 #######################
    tau1_sq = rinvgamma(1, shape =  k/2 + r, rate =  s + sum(gammaZ*beta*beta)/2)
    inv_tau1sq = 1/tau1_sq
    ################## Update GammaZ #######################
    
    const = (qr*sqrt(tau0_sq))/((1 - qr)*sqrt(tau1_sq))
    for (j in 1:p){
      temp1 <- beta[j]^2*0.5*(inv_tau0sq-inv_tau1sq)
      log_dj <- log(const) + temp1
      if(is.infinite(exp(log_dj))== T){
        gammaZ[j] <- 1
      }else{
        gammaZ[j] <- rbinom(1, 1, exp(log_dj)/(1+exp(log_dj)))
      }
      
    }
    
    # update the model size(Xia added)   
    k = sum(gammaZ)
    
    if (itr > nburn) {
      outtau1[(itr-nburn)]<- tau1_sq
      gammamat[(itr-nburn), ] = unlist(gammaZ)
      outbeta[(itr-nburn), ] = unlist(beta)
      outgammaZ <- outgammaZ+gammaZ
      Sizeout[itr-nburn] <- k
    }
  }
  
  est = Maj_sp(betahat=colSums(outbeta)/niter ,Zgamma= outgammaZ/niter,cutoff = 0.5)
  
  return(list("betahat"=est$beta_spr, "gammaZhat" = est$Z_spr,  "gammamat"=gammamat, 
              "Size" = Sizeout,  "betamcmc" = outbeta,"tau1" = outtau1,
              "outgammaZ"=outgammaZ/niter))
  
}



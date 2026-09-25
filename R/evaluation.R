


################################################################################
### Define evaluation Function
################################################################################


#if(!require(caret)) install.packages("caret")
library(caret)

glm_t  = function(X, y, v){ 
  log_lik_t <- function(beta, X, y, v){
    if(is.vector(X) ==T){eta= beta*X} 
    if(is.vector(X) ==F){eta = as.matrix(X)%*%beta}
    p <-  pt(eta, v)
    #xia edit: change to avoid log 0 issue
####    ll <- sum(y * log(p) + (1 - y) * log(1 - p))
    ll <- sum(log(p^y) + log((1 - p)^(1 - y)))
    return(-ll)  # Negative log-likelihood for minimization
  }
  if(is.vector(X) ==T){initial_beta= 0}  # Starting values
  if(is.vector(X) ==F){initial_beta= rep(0, ncol(as.matrix(X)) )}
  betapar =optim(par = initial_beta, fn = log_lik_t, X = X, v=v, y = y, method = "BFGS")
  
  return(coefficients = betapar$par)
}


#######################
evaluation <- function(actual, predicted,betaest, beta,  X, E, v,  method = c("logit","probit", "T")){
  
  true.idx <- which(actual==1)
  false.idx <- which(actual==0)
  positive.idx <- which(predicted==1)
  negative.idx <- which(predicted==0)
  
  TP <- length(intersect(true.idx, positive.idx))
  FP <- length(intersect(false.idx, positive.idx))
  FN <- length(intersect(true.idx, negative.idx))
  TN <- length(intersect(false.idx, negative.idx))
  
  Sensitivity <- TP/(TP+FN)
  if ((TP+FN)==0) Sensitivity <- 1
  
  Specific <- TN/(TN+FP)
  if ((TN+FP)==0) Specific <- 1
  
  MCC.denom <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  if (MCC.denom==0) MCC.denom <- 1
  MCC <- (TP*TN-FP*FN)/MCC.denom
  if ((TN+FP)==0) MCC <- 1
  
  
  if(method == "T"){
   p_pred = pt(X%*%betaest, df=v)
  }
  
  if(method == "probit"){
    p_pred <- pnorm(X%*%betaest)
  } 
  if(method == "logit"){
   p_pred <- plogis(X%*%betaest)
  } 

  E_pred = ifelse(p_pred >=0.5, 1, 0)
  
  
  MSPE <- mean((p_pred-E)^2)
  MAE <-  mean(abs(p_pred-E))
  confusion_matrix <- confusionMatrix(as.factor(c(E_pred)), as.factor(E))
  RMSE_Beta = sum((betaest - beta)^2)/sum(beta*beta)
  
  
  beta_result = list("Sensitivity"=Sensitivity, "Specific"=Specific, "MCC"=MCC, 
                     "TP"=TP, "FP"=FP, "TN"=TN, "FN"=FN, "RMSE_Beta"=RMSE_Beta,
                     "MSPE" = MSPE, "Accuracy" = confusion_matrix$overall[1] )
  
  return("resultbeta" = beta_result)
}





rm(list=ls())
# setwd("~/TSkinny/Case1/codes/")
# All source()/load() paths below are relative to the project root.
source("R/evaluation.R")
source("R/SkinyGibbsT.R")
source("R/ExactGibbsT.R")
source("R/skinny_gibbs_probit.R")
source("R/exact_gibbs_probit.R")
source("R/SGlogit.R")
source("R/EGlogit.R")
##############################################


nsim=50
data_index<-seq(1, (2*nsim))

# the index of the datasets                                                       
#ii=data_index[NAI]     

SkinnyT_results = matrix(rep(0, nsim*10), ncol = 10) ## results summary
ExactT_results =  matrix(rep(0,  nsim*10), ncol = 10) ## results summary
SGLOGIT  = matrix(rep(0, nsim*10), ncol =10) ## results summary
EGLOGIT  = matrix(rep(0, nsim*10), ncol =10) ## results summary
SGprobit = matrix(rep(0, nsim*10), ncol =10) ## results summary
EGprobit = matrix(rep(0, nsim*10), ncol =10) ## results summary
LassoC = matrix(rep(0, nsim*10), ncol =10) ## results summary
LassoL = matrix(rep(0, nsim*10), ncol =10) ## results summary
rLassoL = matrix(rep(0, nsim*10), ncol =10) ## results summary
LassoP = matrix(rep(0, nsim*10), ncol =10) ## results summary


SkinnyT_results_ext = matrix(rep(0, nsim*10), ncol = 10) ## results summary
ExactT_results_ext  =  matrix(rep(0,  nsim*10), ncol = 10) ## results summary
SGLOGIT_ext  = matrix(rep(0, nsim*10), ncol =10) ## results summary
EGLOGIT_ext  = matrix(rep(0, nsim*10), ncol =10) ## results summary
SGprobit_ext = matrix(rep(0, nsim*10), ncol =10) ## results summary
EGprobit_ext = matrix(rep(0, nsim*10), ncol =10) ## results summary
LassoC_ext = matrix(rep(0, nsim*10), ncol =10) ## results summary
LassoL_ext = matrix(rep(0, nsim*10), ncol =10) ## results summary
rLassoL_ext = matrix(rep(0, nsim*10), ncol =10) ## results summary
LassoP_ext = matrix(rep(0, nsim*10), ncol =10) ## results summary

SkinnyT_3 = matrix(rep(0, nsim*10), ncol = 10) ## results summary
ExactT_3 =  matrix(rep(0,  nsim*10), ncol = 10) ## results summary
SkinnyT_7 = matrix(rep(0, nsim*10), ncol = 10) ## results summary
ExactT_7 =  matrix(rep(0,  nsim*10), ncol = 10) ## results summary
SkinnyT_3ext = matrix(rep(0, nsim*10), ncol = 10) ## results summary
ExactT_3ext  =  matrix(rep(0,  nsim*10), ncol = 10) ## results summary
SkinnyT_7ext = matrix(rep(0, nsim*10), ncol = 10) ## results summary
ExactT_7ext  =  matrix(rep(0,  nsim*10), ncol = 10) ## results summary



#data_index
for(ii in data_index){
  
  # Path is relative to the project root and matches the folder structure
  # written by Simulation.R (data/p<p>/n<n>p<p>case<case>/dflist_<i>.RData).
  load(paste0("data/p100/n100p100case1/dflist_",ii,".RData"))
  b0 = 1 ##
  
  if (ii<=nsim){dflist = dflist }
  if (ii>nsim){ dflist=dflist_ext}
  ### setup dataset and parameters
  n0=100 ###number of observations
  pall=dflist$p
  pact = dflist$pact
  set=dflist$set
  case=dflist$case
  qn=dflist$qn
  
  ### Initial number of active covariates (Upper bound)
  if(pact==4) K = 10

  
  v = df = 1
  ind = sample(1:pall, K)
  Z0=rep(0,pall); Z0[ind] = 1 ##gammaZ inital
  ### Model fitting
  
  modelskiny = SkinyGibbsT(X = dflist$Xtrain, E = dflist$Etrain, gammaZ=Z0, b0=1, v=df, r=2, s=1, qr=qn,
                           nburn=2000, niter=4000)

  modelexact =  ExactGibbsT(X = dflist$Xtrain, E = dflist$Etrain, gammaZ=Z0, b0=1, v=df, r=2, s=1, qr=qn,
                            nburn=2000, niter=4000)

  Slogit <- SGlogit(X= as.matrix(dflist$Xtrain), E= unlist(dflist$Etrain),gammaZ=Z0,  b0=1,  r=2, s=1, qr =qn,
                    nburn=2000, niter=4000)

  Elogit <- EGlogit(X= as.matrix(dflist$Xtrain), E= unlist(dflist$Etrain), gamma=Z0, b0=1,r=2, s=1, q =qn,
                    nburn=2000, niter=4000)

  SProbit <- skinnyGibbs(X= as.matrix(dflist$Xtrain), Y= unlist(dflist$Etrain),gamma=Z0,  tau0_2=b0/n0,  r=2, s=1, q =qn,
                         nburn=2000, niter=4000)

  EProbit <- exactGibbs(X= as.matrix(dflist$Xtrain), Y= unlist(dflist$Etrain), gamma=Z0, tau0_2=b0/n0,r=2, s=1, q =qn,
                        nburn=2000, niter=4000)

  lassologit <- glmnet::cv.glmnet(x= as.matrix(dflist$Xtrain), y = unlist(dflist$Etrain), family = "binomial",
                                   alpha = 1, standardize = F)
  lassologit_gam = ifelse(coef(lassologit, s = "lambda.1se")[-1]!= 0, 1, 0)

  lassoprobit <- glmnet::cv.glmnet(x= as.matrix(dflist$Xtrain), y = unlist(dflist$Etrain), alpha = 1,
                                   family = binomial(link = "probit"), standardize = F)
  lassoprobit_gam = ifelse(coef(lassoprobit, s = "lambda.1se")[-1]!= 0, 1, 0)

  lasso_T <- glmnet::cv.glmnet(x= as.matrix(dflist$Xtrain), y = unlist(dflist$Etrain), alpha = 1,
                                    family = binomial(link = "cauchit"), standardize = F)
  lasso_T_gam = ifelse(coef(lasso_T, s = "lambda.1se")[-1]!= 0, 1, 0)

  rlasso <- hdm::rlassologit(x= as.matrix(dflist$Xtrain), y = unlist(dflist$Etrain), post = TRUE,
                                 intercept = FALSE, model = TRUE)
  rlasso_gam = ifelse(coef(rlasso, s = "lambda.1se")!= 0, 1, 0)
  
  
  modelskiny3 = SkinyGibbsT(X = dflist$Xtrain, E = dflist$Etrain, gammaZ=Z0, b0=1, v=3, r=2, s=1, qr=qn, 
                           nburn=2000, niter=4000)  
  
  modelexact3 =  ExactGibbsT(X = dflist$Xtrain, E = dflist$Etrain, gammaZ=Z0, b0=1, v=3, r=2, s=1, qr=qn,
                            nburn=2000, niter=4000)
  
  modelskiny7 = SkinyGibbsT(X = dflist$Xtrain, E = dflist$Etrain, gammaZ=Z0, b0=1, v=7, r=2, s=1, qr=qn, 
                           nburn=2000, niter=4000)  
  
  modelexact7 =  ExactGibbsT(X = dflist$Xtrain, E = dflist$Etrain, gammaZ=Z0, b0=1, v=7, r=2, s=1, qr=qn,
                            nburn=2000, niter=4000)
  
  
  ### Evaluation#########
  gammaz_act = c(rep(1,dflist$pact),rep(0,dflist$p-dflist$pact))
  
  eval_skiny = evaluation(actual=gammaz_act, predicted = modelskiny$gammaZhat, betaest = modelskiny$betahat, beta = dflist$B,
                          X=as.matrix(dflist$Xtest), E = dflist$Etest,v=v, method = "T")

  eval_exact = evaluation(actual=gammaz_act, predicted = modelexact$gammaZhat, betaest = modelexact$betahat, beta = dflist$B,
                          X=as.matrix(dflist$Xtest), E = dflist$Etest,v =v, method = "T")

  eval_Slogit = evaluation(actual=gammaz_act, predicted = Slogit$gammaZhat, betaest= Slogit$betahat, beta = dflist$B,
                           X=as.matrix(dflist$Xtest), E = dflist$Etest,v=v, method = "logit")

  eval_Elogit = evaluation(actual=gammaz_act, predicted = Elogit$gammaZhat, betaest= Elogit$betahat, beta = dflist$B,
                           X=as.matrix(dflist$Xtest), E = dflist$Etest,v=v, method = "logit")

  eval_Sprobit = evaluation(actual=gammaz_act, predicted = SProbit$gammahat, betaest= SProbit$betahat, beta = dflist$B,
                            X=as.matrix(dflist$Xtest), E = dflist$Etest,v=v, method = "probit")

  eval_EGprobit = evaluation(actual=gammaz_act, predicted = EProbit$gammahat, betaest= EProbit$betahat, beta = dflist$B,
                             X=as.matrix(dflist$Xtest), E = dflist$Etest,v=v, method = "probit")

  eval.lassologit = evaluation(actual = gammaz_act, predicted=lassologit_gam, betaest=coef(lassologit, s = "lambda.1se")[-1],
                            beta = dflist$B, X=as.matrix(dflist$Xtest), E = dflist$Etest,v=v, method = "logit")

  eval.rlasso = evaluation(actual = gammaz_act, predicted=rlasso_gam, betaest=coef(rlasso, s = "lambda.1se"),
                                 beta = dflist$B, X=as.matrix(dflist$Xtest), E = dflist$Etest,v=v, method = "logit")

  eval.lasso_T = evaluation(actual = gammaz_act, predicted=lasso_T_gam, betaest=coef(lasso_T, s = "lambda.1se")[-1],
                            beta = dflist$B, X=as.matrix(dflist$Xtest), E = dflist$Etest,v=v, method = "T")

  eval.lassoprobit = evaluation(actual = gammaz_act, predicted=lassoprobit_gam, betaest=coef(lassoprobit, s = "lambda.1se")[-1],
                            beta = dflist$B, X=as.matrix(dflist$Xtest), E = dflist$Etest,v=v, method = "probit")
  
  eval_skiny3 = evaluation(actual=gammaz_act, predicted = modelskiny3$gammaZhat, betaest = modelskiny3$betahat, beta = dflist$B,
                          X=as.matrix(dflist$Xtest), E = dflist$Etest,v=3, method = "T")
  
  eval_exact3 = evaluation(actual=gammaz_act, predicted = modelexact3$gammaZhat, betaest = modelexact3$betahat, beta = dflist$B,
                          X=as.matrix(dflist$Xtest), E = dflist$Etest,v =3, method = "T")
  
  eval_skiny7 = evaluation(actual=gammaz_act, predicted = modelskiny7$gammaZhat, betaest = modelskiny7$betahat, beta = dflist$B,
                          X=as.matrix(dflist$Xtest), E = dflist$Etest,v=7, method = "T")
  
  eval_exact7 = evaluation(actual=gammaz_act, predicted = modelexact7$gammaZhat, betaest = modelexact7$betahat, beta = dflist$B,
                          X=as.matrix(dflist$Xtest), E = dflist$Etest,v =7, method = "T")
  
  if(ii<=nsim){
    SkinnyT_results[ii, ] <- c(unlist(eval_skiny) )
    ExactT_results[ii, ] <- c(unlist(eval_exact) )
    SGprobit[ii, ]<- c(unlist(eval_Sprobit))
    EGprobit[ii, ]<- c(unlist(eval_EGprobit))
    SGLOGIT[ii, ]<- c(unlist(eval_Slogit))
    EGLOGIT[ii, ]<- c(unlist(eval_Elogit))
    LassoC[ii, ] = c(unlist(eval.lasso_T))
    LassoL[ii, ] = c(unlist(eval.lassologit))
    rLassoL[ii, ] = c(unlist(eval.rlasso))
    LassoP[ii, ] = c(unlist(eval.lassoprobit))
    
    SkinnyT_3[ii, ] <- c(unlist(eval_skiny3) )
    ExactT_3[ii, ] <- c(unlist(eval_exact3) )
    SkinnyT_7[ii, ] <- c(unlist(eval_skiny7) )
    ExactT_7[ii, ] <- c(unlist(eval_exact7) )
  }
  
  if (ii > nsim){
    SkinnyT_results_ext[ii-nsim, ] <- c(unlist(eval_skiny) )
    ExactT_results_ext[ii-nsim, ] <- c(unlist(eval_exact) )
    SGprobit_ext[ii-nsim, ]<- c(unlist(eval_Sprobit))
    EGprobit_ext[ii-nsim, ]<- c(unlist(eval_EGprobit))
    SGLOGIT_ext[ii-nsim, ]<- c(unlist(eval_Slogit))
    EGLOGIT_ext[ii-nsim, ]<- c(unlist(eval_Elogit))
    LassoC_ext[ii-nsim, ] = c(unlist(eval.lasso_T))
    LassoL_ext[ii-nsim, ] = c(unlist(eval.lassologit))
    rLassoL_ext[ii-nsim, ] = c(unlist(eval.rlasso))
    LassoP_ext[ii-nsim, ] = c(unlist(eval.lassoprobit))
    SkinnyT_3ext[ii-nsim, ] <- c(unlist(eval_skiny3) )
    ExactT_3ext[ii-nsim, ] <- c(unlist(eval_exact3) )
    SkinnyT_7ext[ii-nsim, ] <- c(unlist(eval_skiny7) )
    ExactT_7ext[ii-nsim, ] <- c(unlist(eval_exact7) )
  }
 
  
}




results_Nooutlier <- data.frame(
  Skiny_T = apply(SkinnyT_results, 2, mean),
  Exact_T = apply(ExactT_results, 2, mean),
  SGLOGIT = apply(SGLOGIT, 2, mean),
  EGLOGIT = apply(EGLOGIT, 2, mean),
  SGProbit = apply(SGprobit, 2, mean),
  EGProbit = apply(EGprobit, 2, mean),
  LASSOC = apply(LassoC, 2, mean),
  LASSOL = apply(LassoL, 2, mean),
  rLASSOL= apply(rLassoL, 2, mean),
  LassoP= apply(LassoP, 2, mean),
   Skiny_T3 = apply(SkinnyT_3, 2, mean),  
   Exact_T3 = apply(ExactT_3, 2, mean),
   Skiny_T7 = apply(SkinnyT_7, 2, mean),  
   Exact_T7 = apply(ExactT_7, 2, mean),
  row.names = c("SEN", "SPE" ,"MCC", "TP", "FP", "TN", "FN", "rMSE", "MSPE", "ACC")) 


results_Outlier <- data.frame(
  Skiny_T = apply(SkinnyT_results_ext, 2, mean),
  Exact_T = apply(ExactT_results_ext, 2, mean),
  SGLOGIT = apply(SGLOGIT_ext, 2, mean),
  EGLOGIT = apply(EGLOGIT_ext, 2, mean),
  SGProb = apply(SGprobit_ext, 2, mean),
  EGProbit = apply(EGprobit_ext[1:18,], 2, mean),
  LASSOC = apply(LassoC_ext, 2, mean),
  LASSOL = apply(LassoL_ext, 2, mean),
  rLASSOL= apply(rLassoL_ext, 2, mean),
  LassoP= apply(LassoP_ext, 2, mean),
  
  Skiny_T3 = apply(SkinnyT_3ext, 2, mean),  
  Exact_T3 = apply(ExactT_3ext, 2, mean),
  Skiny_T7 = apply(SkinnyT_7ext, 2, mean),  
  Exact_T7 = apply(ExactT_7ext, 2, mean),
  row.names = c("SEN", "SPE" ,"MCC", "TP", "FP", "TN", "FN", "rMSE", "MSPE", "ACC")) 


results_Nooutlier
results_Outlier



#write.csv(round(results_Nooutlier,4), paste0("results/tlink/result_100_case1Nooutlier.csv"))
#write.csv(round(results_Outlier,4),   paste0("results/tlink/result_100_case1Outlier.csv"))





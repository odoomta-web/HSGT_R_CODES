

rm(list=ls())
##############################################
# Run this script with the project root as the working directory.
source("R/SIMDATA.R")


nsim= 50 # number of simulation replicates
n0 = 100 # number of observations 
df = 1 # the t-degree of freedom
pall = c(100, 200, 500) #number of covariates
pact0 = c(4,8) #number of active covariates 
case0=c(1,2)# 1 = independent 2 = correlated
set0=c(1, 2, 3) # 1 = 4 active; 2 = 8 active, 3 = constant beta val
# constant beta
beta_val_all=rep(1.5, max(pact0))

Kval = c(10, 20) ## Initial number of active covariates (Upper bound)
#correlations under case0 = 2
rho10 = 0.50; rho20 =0.50 ; rho120 = 0.5 ;
link = c("probit","tlink","logit")


nlink = length(link)
ncase = length(case0)
npact = length(pact0)
npall = length(pall)
nset = length(set0)

# c5 # link (1 = probit, 2= t-link, 3= logit)
# c1 # indepndent or not  (1 = independent, 2 = correlated)
# c2 active covariates.  (1 = 4 = actives, 2 = 8 actives)
# c3 number of covaraites ( 1 = 100ps, 2=200ps, 3=500ps)
# c4 beta values ( 1= postives , 2= +ves and -ves, 3 =constant value)
# i # of replicates. (20 replicates)

for (c5 in c(3)){
  for (c1 in c(1,2)){
    for (c2 in c(1)){
      for (c3 in c(1)){
        for (c4 in c(1)){
          if(set0[c4] == 3){beta_val = beta_val_all[pact0[c2]]}else{
            beta_val = NULL
          }
          
          ## newfolder results
          pathres = paste0("results/p",pall[c3],"/new")
          newfolder11= paste0("n",n0,"p",pall[c3],"case",case0[c1])
          newpath11 = file.path(dirname(pathres), newfolder11)
          dir.create(newpath11)
          ## newfolder data
          path0dat = paste0("data/p",pall[c3],"/new")
          newfolder0= paste0("n",n0,"p",pall[c3],"case",case0[c1])
          newpath0 = file.path(dirname(path0dat), newfolder0)
          dir.create(newpath0)
          
          
          for (i in 1:(nsim)){
            
           
            dflist = simdata(n=n0, p=pall[c3], pact = pact0[c2], 
                             v =df, 
                             set =set0[c4], case = case0[c1], rho12=rho120,
                             rho1 =rho10, rho2 =rho20, 
                             beta_val = beta_val,
                             link = link[c5])
            K = Kval[c2]
            choicep <-function(x){return(x - K + qnorm(0.90)*sqrt((x/pall[c3])*(1- x/pall[c3])))}
            cp = uniroot(choicep, c(1,K))$root
            qn = cp/pall[c3]
        
            dflist$index = i
            dflist$qn = qn
            
            # create an outlier case
            dflist_ext = dflist
            dflist_ext$Xtrain[which(dflist$Etrain==1)[1]] = -10
             
            
           save(dflist, file=paste0("data/p",pall[c3],"/n",n0,"p",pall[c3], "case",case0[c1],"/dflist_",i,".RData"))
           save(dflist_ext, file=paste0("data/p",pall[c3],"/n",n0,"p",pall[c3], "case",case0[c1],"/dflist_",(i+nsim),".RData"))

          } # i # of replicates
        } #c4 beta values
      }# c3 number of covaraites
    } # c2 active covariates
  } # c1 # indepndent or not
} # c5 # link



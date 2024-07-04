#Fitting simulated data for model/parameter recovery
rm(list=ls())
library(rstan)
load("simFits/envE2.RData")
nSub <- 50
simRewards <- env.E2[[1]][1:nSub,,,]
simStim <- env.E2[[2]][1:nSub,,,]
nRd <- dim(simRewards)[2]
nT <- dim(simRewards)[3]
nC <- dim(simRewards)[4]

#Parameters to be used to simulate behaviour
set.seed(290696)
theta <- runif(nSub,-50,50)
omega <- runif(nSub,-50,50)
tauS <- runif(nSub,0.1,10)
dA <- runif(nSub,0,1)
dB <- runif(nSub,0,1)
sticky <- runif(nSub, 0,20)

#Compile simulation model
simExp2.stan <- rstan::stan_model(file= "STANmodels/simulationExp2.stan")

#Prep STAN inputs for simulation 
param.AsySticky <- list(nSub=nSub,nRd=nRd,nT=nT,nC=nC,
                        stim=simStim,R=simRewards,
                        initS0=1,
                        dA=dA,dB=dB,
                        theta=theta,
                        omega=omega,
                        sticky=sticky,
                        tauS=tauS)

param.AsyFull <- list(nSub=nSub,nRd=nRd,nT=nT,nC=nC,
                      stim=simStim,R=simRewards,
                      initS0=1,
                      dA=dA,dB=dB,
                      theta=theta,
                      omega=omega,
                      sticky=rep(0,nSub),
                      tauS=tauS)

param.AsyOmega <- list(nSub=nSub,nRd=nRd,nT=nT,nC=nC,
                       stim=simStim,R=simRewards,
                       initS0=1,
                       dA=dA,dB=dB,
                       theta=rep(0,nSub),
                       omega=omega,
                       sticky=rep(0,nSub),
                       tauS=tauS)

param.AsyTheta <- list(nSub=nSub,nRd=nRd,nT=nT,nC=nC,
                       stim=simStim,R=simRewards,
                       initS0=1,
                       dA=dA,dB=dB,
                       theta=theta,
                       omega=rep(0,nSub),
                       sticky=rep(0,nSub),
                       tauS=tauS)

param.Sticky <- list(nSub=nSub,nRd=nRd,nT=nT,nC=nC,
                     stim=simStim,R=simRewards,
                     initS0=1,
                     dA=rep(1,nSub),dB=rep(1,nSub),
                     theta=theta,
                     omega=omega,
                     sticky=sticky,
                     tauS=tauS)

param.Full <- list(nSub=nSub,nRd=nRd,nT=nT,nC=nC,
                   stim=simStim,R=simRewards,
                   initS0=1,
                   dA=rep(1,nSub),dB=rep(1,nSub),
                   theta=theta,
                   omega=omega,
                   sticky=rep(0,nSub),
                   tauS=tauS)

param.Omega <- list(nSub=nSub,nRd=nRd,nT=nT,nC=nC,
                    stim=simStim,R=simRewards,
                    initS0=1,
                    dA=rep(1,nSub),dB=rep(1,nSub),
                    theta=rep(0,nSub),
                    omega=omega,
                    sticky=rep(0,nSub),
                    tauS=tauS)

param.Theta <- list(nSub=nSub,nRd=nRd,nT=nT,nC=nC,
                    stim=simStim,R=simRewards,
                    initS0=1,
                    dA=rep(1,nSub),dB=rep(1,nSub),
                    theta=theta,
                    omega=rep(0,nSub),
                    sticky=rep(0,nSub),
                    tauS=tauS)

#Save Parmaeters used for simulation  
save(param.AsySticky,param.AsyFull,param.AsyOmega,param.AsyTheta,param.Sticky,param.Full,param.Omega,param.Theta,
     file = "simFits/genParams.RData")

#Generate Behaviour Under Given Model
load("simFits/genParams.RData")
gen.AsySticky <- rstan::sampling(simExp2.stan, data = param.AsySticky, chains = 1, iter = 1, 
                                 seed=29061996,algorithm='Fixed_param')
gen.AsyFull <- rstan::sampling(simExp2.stan, data = param.AsyFull, chains = 1, iter = 1, 
                                 seed=29061996,algorithm='Fixed_param')
gen.AsyOmega <- rstan::sampling(simExp2.stan, data = param.AsyOmega, chains = 1, iter = 1, 
                                 seed=29061996,algorithm='Fixed_param')
gen.AsyTheta <- rstan::sampling(simExp2.stan, data = param.AsyTheta, chains = 1, iter = 1, 
                                 seed=29061996,algorithm='Fixed_param')
gen.Sticky <- rstan::sampling(simExp2.stan, data = param.Sticky, chains = 1, iter = 1, 
                                 seed=29061996,algorithm='Fixed_param')
gen.Full <- rstan::sampling(simExp2.stan, data = param.Full, chains = 1, iter = 1, 
                              seed=29061996,algorithm='Fixed_param')
gen.Omega <- rstan::sampling(simExp2.stan, data = param.Omega, chains = 1, iter = 1, 
                              seed=29061996,algorithm='Fixed_param')
gen.Theta <- rstan::sampling(simExp2.stan, data = param.Theta, chains = 1, iter = 1, 
                              seed=29061996,algorithm='Fixed_param')
#Save generated behaviour
save(gen.AsySticky,gen.AsyFull,gen.AsyOmega,gen.AsyTheta,gen.Sticky,gen.Full,gen.Omega,gen.Theta,
     file = "simFits/genBehaviour.RData")

load("simFits/genBehaviour.RData")

#Function to fit models. Takes generated object, extracts and fits using model PATH 
func.Fit <- function(nSub,nRd,nT,nC,genObj,stanModel){
  exGenObj <- rstan::extract(genObj)
  dataGenObj <- list(nSub=nSub,nRd=nRd,nT=nT,nC=nC,
                     choice=exGenObj$y_pred[1,,,],
                     stim=exGenObj$stimulation[1,,,],
                     Rmu=exGenObj$Rmu_ts[1,,,,])
  fitObj <- rstan::stan(file=stanModel,data=dataGenObj,
                        chains = 4, cores=4, 
                        warmup = 2000, iter=4000,
                        seed=29061996,
                        verbose=F, save_warmup=F,
                        control = list(adapt_delta=.99))
  return(fitObj)
}

#Loop to fit all models to all generated datasets

Models <- c("AsySticky","AsyFull","AsyOmega","AsyTheta","Sticky","Full","Omega","Theta")

for (gen in 1:length(Models)){
  for (fit in 1:length(Models)){
    genObj.Name <- paste0("gen.",Models[gen])
    genObj <- get(genObj.Name)
    fitPath <- paste0("STANmodels/hierarchicalE2_",Models[fit],".stan")
    
    print(paste("Fitting Model",fitPath,"to",genObj.Name))
    
    fit.file <- func.Fit(nSub,nRd,nT,nC,genObj,fitPath)
    saveRDS(fit.file,file=paste0("simFits/gen",Models[gen],"_fit",Models[fit],".fit"))
    
  }
}



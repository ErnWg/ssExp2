#Fitting simulated data for model/parameter recovery
rm(list=ls())
library(rstan)
load("simFits/experimentSim.RData")
nSub <- 10
simRewards <- Env[[3]][1:nSub,,,]
simStim <- Env[[4]][1:nSub,,,]
nRd <- dim(simRewards)[2]
nT <- dim(simRewards)[3]
nB <- dim(simRewards)[4]
opt1 <- Env[[1]][1:nSub,,]
opt2 <- Env[[2]][1:nSub,,]

#Parameters to be used to simulate behaviour
theta <- param.simFull$theta[1:nSub]
omega <- param.simFull$omega[1:nSub]
tau <- param.simFull$tau[1:nSub]

#Function to fit models. Takes generated object, extracts and fits using model PATH 
func.Fit <- function(nSub,nRd,nP1,nP2,nT,nB,opt1,opt2,genObj,stanModel){
  exGenObj <- rstan::extract(genObj)
  dataGenObj <- list(nSub=nSub,nRd=nRd,nP1=nP1,nP2=nP2,nT=nT,nB=nB,
                     opt1=opt1,opt2=opt2,
                     choice=exGenObj$y_pred[1,1:nSub,,],
                     Rmu=exGenObj$Rmu_ts[1,1:nSub,,,],
                     Smu=exGenObj$Smu_ts[1,1:nSub,,,],
                     Ssig=exGenObj$Ssig_ts[1,1:nSub,,,])
  fitObj <- rstan::stan(file=stanModel,data=dataGenObj,
                        chains = 4, cores=4, 
                        warmup = 2000, iter=4000,
                        seed=29061996,
                        verbose=F, save_warmup=F,
                        control = list(adapt_delta=.85))
  return(fitObj)
}

#Loop to fit all models to all generated datasets

Models <- c("Full","Omega","Theta")

for (gen in 1:length(Models)){
  for (fit in 1:length(Models)){
    genObj.Name <- paste0("sim.",Models[gen])
    genObj <- get(genObj.Name)
    fitPath <- paste0("STANmodels/hierarchical",Models[fit],".stan")
    
    print(paste("Fitting Model",fitPath,"to",genObj.Name))
    
    fit.file <- func.Fit(nSub,nRd,nP1=2,nP2=2,nT,nB,opt1,opt2,genObj,fitPath)
    saveRDS(fit.file,file=paste0("simFits/gen",Models[gen],"_fit",Models[fit],".fit"))
    
  }
}



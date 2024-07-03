#This version creates an experimental environment more aligned with the norbury task
rm(list=ls())

packages <- c("plyr", 'dplyr','tidyr',  "ggplot2", 'pander', 'emmeans','sjPlot', 'lmerTest', 'data.table', "ggbeeswarm",
              'lmerTest',  'brms', 'BayesFactor', "gridExtra", "lsr",'ggridges','cowplot','entropy','zoo','Hmisc', 
              'DescTools', 'magrittr', 'purrr', 'forcats', 'modelr', 'tidybayes', 'rstan' , 'loo', 'R.matlab','ggpubr')
invisible(lapply(packages, require, character.only = TRUE))

sigma <- 2.5
R.true <- c(0,10)
S.true <- c(0.15,0.5,0.85)
nP1 <- 30
nP2 <- 30
nT <- nP1 + nP2
nB = 6

dfParams <- expand.grid(reward=R.true,stim=S.true)
dfParams$rewardNoise <- sigma

dfComb <- data.frame(t(combn(1:nB,2))) #X1 is right bandit, X2 is left bandit

createExp <- function(nSub,nRd,repP1,repP2,dfParams,dfComb){
  nB <- nrow(dfParams)
  nT <- nrow(dfComb) * (repP1 + repP2)
  rewards <- array(NA, dim=c(nSub,nRd,nT,nB))
  stim <- array(NA, dim=c(nSub,nRd,nT,nB))
  opt1 <- array(NA, dim=c(nSub,nRd,nT))
  opt2 <- array(NA, dim=c(nSub,nRd,nT))
  for (s in 1:nSub){
    for (rd in 1:nRd){
      #Create trial choices
      trials.P1 <- do.call("rbind",replicate(repP1,dfComb,simplify=F))
      trials.P1 <- trials.P1[sample(nrow(trials.P1)), ]
      trials.P2 <- do.call("rbind",replicate(repP2,dfComb,simplify=F))
      trials.P2 <- trials.P2[sample(nrow(trials.P2)), ]
      trials.full <- rbind(trials.P1,trials.P2)
      opt1[s,rd,] <- unlist(trials.full$X1)
      opt2[s,rd,] <- unlist(trials.full$X2)
      
      #Create Reward and Stim Environment
      for (b in 1:nB){
        rewards[s,rd,,b] <- round(rnorm(nT,dfParams$reward[b],dfParams$rewardNoise[b]))
        stim[s,rd,,b] <- rbinom(nT,size=1,dfParams$stim[b])
      }
    }
  }
  return(list(opt1,opt2,rewards,stim))
}


nSub <- 100
set.seed(29061996)
Env <- createExp(nSub,4,2,2,dfParams,dfComb)

Sim.opt1 <- Env[[1]]
Sim.opt2 <- Env[[2]]
Sim.rewards <- Env[[3]]
Sim.stim <- Env[[4]]


nRd <- 4
nP1 <- 30
nP2 <- 30
nT <- 60
nB <- 6
noise <- 2.5
opt1 <- Sim.opt1
opt2 <- Sim.opt2
stim <- Sim.stim
rewards <- Sim.rewards
initR0 <- 5
set.seed(29061996)
theta <- runif(nSub,-50,50)
omega <- runif(nSub,-50,50)
tau <- runif(nSub,0.1,10)
dA <- runif(nSub,0,1)
dB <- runif(nSub,0,1)

simBehaviour.stan <- rstan::stan_model(file= "STANmodels/v2simulation.stan")

#Parameters for Full model simulation 
param.simFull <- list(nSub=nSub,nRd=nRd,nP1=nP1,nP2=nP2,nT=nT,nB=nB,
                      noise=sigma,
                      opt1=Sim.opt1,opt2=Sim.opt2,
                      stim=Sim.stim,rewards=Sim.rewards,
                      initR0=5,
                      dA=rep(1,nSub),dB=rep(1,nSub),
                      theta=theta,omega=omega,tau=tau)

sim.Full <- rstan::sampling(simBehaviour.stan, data = param.simFull, chains = 1, iter = 1, algorithm='Fixed_param')

#Parameters for Omega model simulation
param.simOmega <- list(nSub=nSub,nRd=nRd,nP1=nP1,nP2=nP2,nT=nT,nB=nB,
                      noise=sigma,
                      opt1=Sim.opt1,opt2=Sim.opt2,
                      stim=Sim.stim,rewards=Sim.rewards,
                      initR0=5,
                      dA=rep(1,nSub),dB=rep(1,nSub),
                      theta=rep(0,nSub),omega=omega,tau=tau)

sim.Omega <- rstan::sampling(simBehaviour.stan, data = param.simOmega, chains = 1, iter = 1, algorithm='Fixed_param')

#Parameters for Theta model simulation
param.simTheta <- list(nSub=nSub,nRd=nRd,nP1=nP1,nP2=nP2,nT=nT,nB=nB,
                       noise=sigma,
                       opt1=Sim.opt1,opt2=Sim.opt2,
                       stim=Sim.stim,rewards=Sim.rewards,
                       initR0=5,
                       dA=rep(1,nSub),dB=rep(1,nSub),
                       theta=theta,omega=rep(0,nSub),tau=tau)

sim.Theta <- rstan::sampling(simBehaviour.stan, data = param.simTheta, chains = 1, iter = 1, algorithm='Fixed_param')

save(Env,param.simFull,sim.Full,param.simOmega,sim.Omega,param.simTheta,sim.Theta, file = "simFits/experimentSim.RData")

#load("simFits/experimentSim.RData")

simEx.Full <- rstan::extract(sim.Full)
simEx.Omega <- rstan::extract(sim.Omega)
simEx.Theta <- rstan::extract(sim.Theta)

#iN 

#Frequency of options
choice.Full <- simEx.Full$chosenObj[1,,,]
prop.table(table(choice.Full))

choice.Omega <- simEx.Omega$chosenObj[1,,,]
prop.table(table(choice.Omega))

choice.Theta <- simEx.Theta$chosenObj[1,,,]
prop.table(table(choice.Theta))


Rmu <- simEx.Full$Rmu_ts[1,,,,]
Rsig <- simEx.Full$Rsig_ts[1,,,,]
Smu <- simEx.Full$Smu_ts[1,,,,]
Ssig <- simEx.Full$Ssig_ts[1,,,,]
simEx.Full$stimDeliver[1,1,1,]




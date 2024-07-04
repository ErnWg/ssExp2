#Author: Ern Wong
#This code runs simulations + parameter/model recovery analysis on an extension of the Norbury task.
#Briefly, its a 3 armed extension.

#Load packages
rm(list=ls())
source("functions.R")
packages <- c("plyr", 'dplyr','tidyr',  "ggplot2", 'pander', 'emmeans','sjPlot', 'lmerTest', 'data.table', "ggbeeswarm",
              'lmerTest',  'brms', 'BayesFactor', "gridExtra", "lsr",'ggridges','cowplot','entropy','zoo','Hmisc', 
              'DescTools', 'magrittr', 'purrr', 'forcats', 'modelr', 'tidybayes', 'rstan' , 'loo', 'R.matlab','ggpubr','abind')
invisible(lapply(packages, require, character.only = TRUE))
#Set up STAN environment
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())


#Define environments as a dataframe
paramEnv1 <- data.frame(Environment = c(1,1,1),
                        Choice = c(1,2,3),
                        Rmean = c(0,10,10),
                        Rvar = c(0,0,0),
                        prob = c(0.5,0.9,0.1),
                        stimLearn = c(1,1,1))

paramEnv2 <- data.frame(Environment = c(2,2,2),
                        Choice = c(1,2,3),
                        Rmean = c(0,10,10),
                        Rvar = c(0,0,0),
                        prob = c(0.9,0.1,0.5),
                        stimLearn = c(1,1,1))

paramEnv3 <- data.frame(Environment = c(3,3,3),
                        Choice = c(1,2,3),
                        Rmean = c(0,10,10),
                        Rvar = c(0,0,0),
                        prob = c(0.1,0.5,0.9),
                        stimLearn = c(1,1,1))

paramEnv4 <- data.frame(Environment = c(4,4,4),
                        Choice = c(1,2,3),
                        Rmean = c(10,0,0),
                        Rvar = c(0,0,0),
                        prob = c(0.5,0.9,0.1),
                        stimLearn = c(1,1,1))

paramEnv5 <- data.frame(Environment = c(5,5,5),
                        Choice = c(1,2,3),
                        Rmean = c(10,0,0),
                        Rvar = c(0,0,0),
                        prob = c(0.9,0.1,0.5),
                        stimLearn = c(1,1,1))

paramEnv6 <- data.frame(Environment = c(6,6,6),
                        Choice = c(1,2,3),
                        Rmean = c(10,0,0),
                        Rvar = c(0,0,0),
                        prob = c(0.1,0.5,0.9),
                        stimLearn = c(1,1,1))

paramEnv <- rbind(paramEnv1,paramEnv2,paramEnv3,paramEnv4,paramEnv5,
                  paramEnv6)

set.seed(290696)
env.E2 <- createEnv2(100,4,20,paramEnv)
save(env.E2, file = "simFits/envE2.RData")

#Prep for simulation 
nSub = 10
load("simulations/simExp2.RData")
simRewards <-simExp2[[1]][1:nSub,,,]
simStim <- simExp2[[2]][1:nSub,,,]
#Rmean <- simExp2[[3]]
#Rnoise <- sqrt(simExp2[[4]])
#stimLearn <-t(replicate(nSub,simExp2[[5]]))
set.seed(290696)
theta <- runif(nSub,-50,50)
omega <- runif(nSub,-50,50)
tauS <- runif(nSub,0.1,10)
dA <- runif(nSub,0,1)
dB <- runif(nSub,0,1)
sticky <- runif(nSub, 0,20)

param.Exp2 <- list(nSub=nSub,nRd=24,nT=20,nC=3,
                      stim=simStim,R=simRewards,
                      initS0=1,
                      dA=dA,dB=dB,
                      theta=theta,
                      omega=omega,
                      sticky=rep(0,nSub),
                      tauS=tauS)

simExp2.stan <- rstan::stan_model(file= "STANmodels/simulationExp2.stan")
sim.Exp2 <- rstan::sampling(simExp2.stan, data = param.Exp2, chains = 1, iter = 1, algorithm='Fixed_param')
simEx.Exp2 <- rstan::extract(sim.Exp2)

data.Exp2 <- list(nSub=10,nRd=24,nT=20,nC=3,
                     choice=simEx.Exp2$y_pred[1,,,],
                     stim=simEx.Exp2$stimulation[1,,,],
                     Rmu=simEx.Exp2$Rmu_ts[1,,,,])

fit.Exp2 <- rstan::stan(file='STANmodels/hierarchicalE2_Asy.stan', data = data.Exp2, chains = 4, cores=4, verbose=F)
saveRDS(fit.Exp2,"Exp2_Asy.fit")
fit.Exp2 <- readRDS("Exp2_full.fit")

ex.exp2 <- rstan::extract(fit.Exp2)
full <- data.frame(simTheta = theta,
                   simOmega = omega,
                   rTheta = colMeans(ex.exp2$theta),
                   rOmega = colMeans(ex.exp2$omega))

library("ggpubr")
ggscatter(full, x = "simOmega", y = "rOmega", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Simulated Omega", ylab = "Recovered Omega")

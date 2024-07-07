rm(list=ls())

packages <- c("plyr", 'dplyr','tidyverse',  "ggplot2", 'pander', 'emmeans','sjPlot', 'lmerTest', 'data.table', "ggbeeswarm",
              'lmerTest',  'brms', 'BayesFactor', "gridExtra", "lsr",'ggridges','cowplot','entropy','zoo','Hmisc', 
              'DescTools', 'magrittr', 'purrr', 'forcats', 'modelr', 'tidybayes', 'rstan' , 'loo', 'R.matlab','ggpubr','gridExtra')
invisible(lapply(packages, require, character.only = TRUE))

load("simFits/experimentSim.RData")
nSub <- 50
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


genFull_fitFull <- readRDS("simFits/genFull_fitFull.fit")
genFull_fitOmega <- readRDS("simFits/genFull_fitOmega.fit")
genFull_fitTheta <- readRDS("simFits/genFull_fitTheta.fit")

genOmega_fitFull <- readRDS("simFits/genOmega_fitFull.fit")
genOmega_fitOmega <- readRDS("simFits/genOmega_fitOmega.fit")
genOmega_fitTheta <- readRDS("simFits/genOmega_fitTheta.fit")

genTheta_fitFull <- readRDS("simFits/genTheta_fitFull.fit")
genTheta_fitOmega <- readRDS("simFits/genTheta_fitOmega.fit")
genTheta_fitTheta <- readRDS("simFits/genTheta_fitTheta.fit")

#Extract

ex.genFull_fitFull <- rstan::extract(genFull_fitFull)
ex.genOmega_fitOmega <- rstan::extract(genOmega_fitOmega)
ex.genTheta_fitTheta <- rstan::extract(genTheta_fitTheta)

rFull <- data.frame(simTheta = theta,
                    simOmega = omega,
                    simTau = tau,
                    rTheta = colMeans(ex.genFull_fitFull$theta),
                    rOmega = colMeans(ex.genFull_fitFull$omega),
                    rTau = colMeans(ex.genFull_fitFull$tau))

rOmega <- data.frame(simOmega = omega,
                     simTau = tau,
                     rOmega = colMeans(ex.genOmega_fitOmega$omega),
                     rTau = colMeans(ex.genOmega_fitOmega$tau))

rTheta <- data.frame(simTheta = theta,
                     simTau = tau,
                     rTheta = colMeans(ex.genTheta_fitTheta$theta),
                     rTau = colMeans(ex.genTheta_fitTheta$tau))




ggscatter(rFull, x = "simTheta", y = "rTheta", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Simulated Theta", ylab = "Recovered Theta")

ggscatter(rFull, x = "simOmega", y = "rOmega", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Simulated Omega", ylab = "Recovered Omega")

ggscatter(rFull, x = "simTau", y = "rTau", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Simulated Tau", ylab = "Recovered Tau")


#Model Comparison

loo.Obj <- loo(genFull_fitFull)

loo.Obj$pointwise[,4]


loo(genFull_fitOmega)
loo(genFull_fitTheta)







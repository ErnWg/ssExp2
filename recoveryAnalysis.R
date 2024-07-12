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
          xlab = "Simulated \u03b8", ylab = "Recovered \u03b8",
          color="#0070c0",cor.coef.size = 7) +
  theme_classic() + 
  theme(text=element_text(size=20))

ggscatter(rFull, x = "simOmega", y = "rOmega", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Simulated \u03c9", ylab = "Recovered \u03c9",
          color="#cc79a7",cor.coef.size = 7) +
          theme_classic() + 
          theme(text=element_text(size=20))

ggscatter(rFull, x = "simTau", y = "rTau", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Simulated \u03c4", ylab = "Recovered \u03c4",
          color="#d55e00",cor.coef.size = 7) +
          theme_classic() + 
          theme(text=element_text(size=20))

#Parm recovery for Omega model
ggscatter(rOmega, x = "simOmega", y = "rOmega", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Simulated \u03c9", ylab = "Recovered \u03c9",
          color="#cc79a7",cor.coef.size = 7) +
  theme_classic() + 
  theme(text=element_text(size=20))



#Model Comparison
library("Rfast")
#For genFull
loo.genFull_fitFull <- loo(genFull_fitFull)
loo.genFull_fitTheta <- loo(genFull_fitTheta)
loo.genFull_fitOmega <- loo(genFull_fitOmega)
#For genTheta
loo.genTheta_fitFull <- loo(genTheta_fitFull)
loo.genTheta_fitTheta <- loo(genTheta_fitTheta)
loo.genTheta_fitOmega <- loo(genTheta_fitOmega)
#For genOmega
loo.genOmega_fitFull <- loo(genOmega_fitFull)
loo.genOmega_fitTheta <- loo(genOmega_fitTheta)
loo.genOmega_fitOmega <- loo(genOmega_fitOmega)

fitMat <- matrix(NA,3,3)

fits.genFull <- matrix(0,nSub,3) 
fits.genFull[,1]<-loo.genFull_fitFull$pointwise[,4]
fits.genFull[,2]<-loo.genFull_fitTheta$pointwise[,4]
fits.genFull[,3]<-loo.genFull_fitOmega$pointwise[,4]
fitMat[1,] <- table(rowMins(fits.genFull))

fits.genTheta <- matrix(0,nSub,3) 
fits.genTheta[,1]<-loo.genTheta_fitFull$pointwise[,4]
fits.genTheta[,2]<-loo.genTheta_fitTheta$pointwise[,4]
fits.genTheta[,3]<-loo.genTheta_fitOmega$pointwise[,4]
fitMat[2,] <- table(rowMins(fits.genTheta))

fits.genOmega <- matrix(0,nSub,3) 
fits.genOmega[,1]<-loo.genOmega_fitFull$pointwise[,4]
fits.genOmega[,2]<-loo.genOmega_fitTheta$pointwise[,4]
fits.genOmega[,3]<-loo.genOmega_fitOmega$pointwise[,4]
fitMat[3,] <- table(rowMins(fits.genOmega))

fitMat
confusionMat <- fitMat/50

invMat <- matrix(NA,3,3)

for (i in 1:nrow(fitMat)){
  for (j in 1:ncol(fitMat)){
    invMat[i,j] <- fitMat[i,j]/sum(fitMat[,j])
  }
} 


#DIC Plots
Models <- c("Full","Omega","Theta")
simData <- c("Full","Omega","Theta")

logList <- list()
nLLList <- list()
DIClist <- list()

fitMat <- matrix(NA,length(simData),length(Models))
for (gen in 1:length(simData)){
  fitList <- matrix(0,nSub,length(Models))
  for (fit in 1:length(Models)){
    
    fitObj <- get(paste0("gen",simData[gen],"_fit",Models[fit]))
    ex.fitObj <- rstan::extract(fitObj)
    
    deviance <-  -2*colMeans(ex.fitObj$log_lik)
    varDeviance <- apply(ex.fitObj$log_lik, MARGIN = 2, FUN=function(i) var(-2*i))
    fitList[,fit] = deviance + (0.5*varDeviance)
  }
  fitMat[gen,] <- tabulate(rowMins(fitList))
}

confusionMat <- fitMat/50

#invMat <- matrix(NA,3,3)
#
#for (i in 1:nrow(fitMat)){
#  for (j in 1:ncol(fitMat)){
#    invMat[i,j] <- fitMat[i,j]/sum(fitMat[,j])
#  }
#} 


plotTable <- function(mat,name){
df <- data.frame(round(mat,2))
rownames(df) <- c("Full","Omega","Theta")
colnames(df) <- c("Full","Omega","Theta")

df2 <- df %>% rownames_to_column() %>% gather(colname, value, -rowname) 
head(df2)

matFig <- ggplot(df2, aes(x = rowname, y = colname, fill = value)) +
  geom_tile(color = "black") +
  geom_text(aes(label = round(value,3)), color = "black", size = 6) +
  labs(title=name,x="Recovered", y = "Simulated",size=15) +
  scale_fill_distiller() +
  #scale_fill_viridis_c(option="magma",guide = "none") +
  theme(axis.text=element_text(size=15),
        plot.title = element_text(size = 20, hjust = 0.5),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        legend.position = 'none')+
  coord_fixed()

matFig
}

plotTable(confusionMat, "Confusion Matrix")
plotTable(invMat, "Inversion Matrix")


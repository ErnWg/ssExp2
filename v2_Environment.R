#This version creates an experimental environment more aligned with the norbury task
rm(list=ls())

packages <- c("plyr", 'dplyr','tidyverse',  "ggplot2", 'pander', 'emmeans','sjPlot', 'lmerTest', 'data.table', "ggbeeswarm",
              'lmerTest',  'brms', 'BayesFactor', "gridExtra", "lsr",'ggridges','cowplot','entropy','zoo','Hmisc', 
              'DescTools', 'magrittr', 'purrr', 'forcats', 'modelr', 'tidybayes', 'rstan' , 'loo', 'R.matlab','ggpubr','gridExtra')
invisible(lapply(packages, require, character.only = TRUE))

sigma <- 2.5
R.true <- c(0,10)
S.true <- c(0.15,0.5,0.85)
nP1 <- 30
nP2 <- 30
nT <- nP1 + nP2
nB <- 6

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
theta <- runif(nSub,-5,5)
omega <- runif(nSub,-50,50)
tau <- runif(nSub,0.1,5)
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

#Frequency of options
choice.Full <- simEx.Full$chosenObj[1,,,]
optimal.Freq <- prop.table(table(choice.Full))

#From full model params, index high low theta/omega
lowTheta <- param.simFull$theta < median(param.simFull$theta)
highTheta <- param.simFull$theta > median(param.simFull$theta)
lowOmega <- param.simFull$omega < median(param.simFull$omega)
highOmega <- param.simFull$omega > median(param.simFull$omega)

lowTheta.Freq <- prop.table(table(choice.Full[lowTheta,,]))
highTheta.Freq <- prop.table(table(choice.Full[highTheta,,]))

lowOmega.Freq <- prop.table(table(choice.Full[lowOmega,,]))
highOmega.Freq <- prop.table(table(choice.Full[highOmega,,]))

plotFreq <- function(table,name){
  mat <- matrix(NA, nrow = 3, ncol = 2)
  #row if stim, col is reward
  mat[3,1] <- table[1]
  mat[3,2] <- table[2]
  mat[2,1] <- table[3]
  mat[2,2] <- table[4]
  mat[1,1] <- table[5]
  mat[1,2] <- table[6]
  df <- data.frame(mat)
  rownames(df) <- c("0.85","0.5","0.15")
  colnames(df) <- c("Low","High")
  
  df2 <- df %>% rownames_to_column() %>% gather(colname, value, -rowname) 
  head(df2)
  
  plotFig <- ggplot(df2, aes(x = rowname, y = colname, fill = value)) +
    geom_tile(color = "black") +
    geom_text(aes(label = round(value,3)), color = "green", size = 4) +
    labs(title=name,x="Stimulation Probabilities", y = "Reward Level") +
    scale_fill_viridis_c(option = "magma", guide = "none") +
    coord_fixed()

  return(plotFig)
}

plotFreq(optimal.Freq,"Optimal")
plotFreq(lowOmega.Freq,"Low Omega")
plotFreq(highOmega.Freq,"High Omega")
plotFreq(lowTheta.Freq,"Low Theta")
plotFreq(highTheta.Freq,"High Theta")

#Function for learning curves

plotLearnCurve <- function(data,name){
  LearnCurves <- data.frame(data)
  LearnCurves$Trial <- 1:nrow(data)

  LearnCurves.long <- gather(LearnCurves, Bandit, Reward, X1:X6)
    learningPlot <- ggplot(LearnCurves.long,aes(y=Reward,x=Trial,color=Bandit)) +
    labs(title=name,x="Trial Number", y = name) +
    geom_rect(aes(xmin=(nP1+.5), xmax=(nP1+nP2), ymin=-Inf, ymax=Inf),fill="grey",alpha=.1) +
    geom_line(linetype = "solid")+
    geom_point(size=1) +
    scale_fill_viridis_c(option ="viridis")+
    #scale_colour_manual(values = c("#cc79a7","#0070c0","#019e73"),name="Bandits",labels = c("Left","Down","Right")) +
    theme(legend.position= c(1,1), legend.justification = c(1,1),text=element_text(size=15)) + theme(legend.position = 'none') +
    annotate("text", x = c(nP1/2,nP2/2+nP1), y = c(1.1*max(data),1.1*max(data)), label = c("Phase 1", "Phase 2 \n (Sensory Stimulation)"), 
            color="Black", size=3 , angle=0, fontface="bold") +
    theme_classic() 

  learningPlot
  return(learningPlot)
}

#Plot Learning Curves
plotLearnCurve(simEx.Full$Rmu_ts[1,1,3,,],"Estimated Reward")
plotLearnCurve(simEx.Full$Rsig_ts[1,1,3,,],"Estimated Reward Uncertainty")
plotLearnCurve(simEx.Full$Smu_ts[1,1,3,,],"Estimated Stimulation Probability")
plotLearnCurve(simEx.Full$Ssig_ts[1,1,3,,],"Estimated Stimulation Uncertainty")


#Plot average score received

scores <- simEx.Full$Score[1,,,]
apply(scores,1,sum)/4

scores.Omega <- simEx.Omega$Score[1,,,]
apply(scores.Omega,1,sum)/4


rm(list=ls())

packages <- c("plyr", 'dplyr','tidyr',  "ggplot2", 'pander', 'emmeans','sjPlot', 'lmerTest', 'data.table', "ggbeeswarm",
              'lmerTest',  'brms', 'BayesFactor', "gridExtra", "lsr",'ggridges','cowplot','entropy','zoo','Hmisc', 
              'DescTools', 'magrittr', 'purrr', 'forcats', 'modelr', 'tidybayes', 'rstan' , 'loo', 'R.matlab','ggpubr')
invisible(lapply(packages, require, character.only = TRUE))


### Plot Experimental environment 
#Parameters
lowR <- 0
highR <- 10
nP1 <- 100
nP2 <- 1
sigma <- 2
nT <- nP1 + nP2

#set.seed(29101996)
lowReward <- data.frame(b1 = round(rnorm(nT,lowR,sigma)),
                        b2 = round(rnorm(nT,lowR,sigma)),
                        b3 = round(rnorm(nT,highR,sigma)),
                        Trial = 1:nT)
highReward <- data.frame(b1 = round(rnorm(nT,lowR,sigma)),
                         b2 = round(rnorm(nT,highR,sigma)),
                         b3 = round(rnorm(nT,highR,sigma)),
                         Trial = 1:nT)
lowReward.long <- gather(lowReward, Bandit, Reward, b1:b3)
highReward.long <- gather(highReward, Bandit, Reward, b1:b3)

exEnv <- ggplot(highReward.long,aes(y=Reward,x=Trial,color=Bandit)) +
  labs(title="High Reward Environment",x="Trial Number", y = "Unshifted Rewards") +
  geom_rect(aes(xmin=(nP1+.5), xmax=(nP1+nP2), ymin=-Inf, ymax=Inf),fill="grey",alpha=.1) +
  geom_line(linetype = "solid")+
  geom_point(size=2) +
  scale_colour_manual(values = c("#cc79a7","#0070c0","#019e73"),name="Bandits",labels = c("Left","Down","Right")) +
  theme(legend.position= c(1,1), legend.justification = c(1,1),text=element_text(size=15)) + theme(legend.position = 'none') +
  annotate("text",x = c(nP1/2,nP2/2+nP1), y = c(max(highReward.long$Reward) +4,max(highReward.long$Reward) +4), label = c("Phase 1 \n (Complete Reward Feedback)", "Phase 2 \n (Sensory + No Reward Feedback)") , 
           color="Black", size=3 , angle=0, fontface="bold") +
  theme_classic() 

exEnv


### Plot Learning Curves for Reward Distribution according to Kalman Filter Model
initU0 <- sqrt(sigma^2 * 20)
initR0 <- (lowR + highR)/2
funcKF <- function(initR0,initU0,sigma,obs,nT_P1,nT_P2){
  Rmu <- rep(initR0,nT_P1)
  Rsig <- rep(initU0,nT_P1)
  for(p1 in 1:nT_P1){
    delta = obs[p1] - Rmu[p1] 
    kgain = Rsig[p1]^2 / (Rsig[p1]^2 + sigma^2)
    
    Rmu[p1+1] = Rmu[p1] + kgain * delta
    Rsig[p1+1] = sqrt((1-kgain)*Rsig[p1]^2)
  }
  Rmu_ts <- c(Rmu,rep(Rmu[nT_P1+1],nT_P2-1))
  Rsig_ts <- c(Rsig,rep(Rsig[nT_P1+1],nT_P2-1))
  
  return(list(Rmu_ts,Rsig_ts))
}
  
rewardCurves <- data.frame(b1 = funcKF(initR0,initU0,sigma,highReward$b1,nP1,nP2)[[1]],
                           b2 = funcKF(initR0,initU0,sigma,highReward$b2,nP1,nP2)[[1]],
                           b3 = funcKF(initR0,initU0,sigma,highReward$b3,nP1,nP2)[[1]],
                           Trial = 1:(nP1+nP2))
rewardCurves.long <- gather(rewardCurves, Bandit, Reward, b1:b3)

learningPlot <- ggplot(rewardCurves.long,aes(y=Reward,x=Trial,color=Bandit)) +
  labs(title="High Reward Environment",x="Trial Number", y = "Estimated Rewards") +
  geom_rect(aes(xmin=(nP1+.5), xmax=(nP1+nP2), ymin=-Inf, ymax=Inf),fill="grey",alpha=.1) +
  geom_line(linetype = "solid")+
  geom_point(size=2) +
  scale_colour_manual(values = c("#cc79a7","#0070c0","#019e73"),name="Bandits",labels = c("Left","Down","Right")) +
  theme(legend.position= c(1,1), legend.justification = c(1,1),text=element_text(size=15)) + theme(legend.position = 'none') +
  annotate("text", x = c(nP1/2,nP2/2+nP1), y = c(10,10), label = c("Phase 1 \n (Complete Reward Feedback)", "Phase 2 \n (Sensory + No Reward Feedback)") , 
           color="Black", size=3 , angle=0, fontface="bold") +
  theme_classic() 

learningPlot

uncertaintyCurves <- data.frame(b1 = funcKF(initR0,initU0,sigma,highReward$b1,nP1,nP2)[[2]],
                                Trial = 1:(nP1+nP2))
uncertaintyPlot <- ggplot(uncertaintyCurves,aes(y=b1,x=Trial)) +
  labs(title="Uncertainty Curve",x="Trial Number", y = "Estimation Uncertainty") +
  geom_rect(aes(xmin=(nP1+.5), xmax=(nP1+nP2), ymin=-Inf, ymax=Inf),fill="grey",alpha=.1) +
  geom_line(linetype = "solid")+
  geom_point(size=2) +
  theme(legend.position=F)+
  annotate("text", x = c(nP1/2,nP2/2+nP1), y = c(16,16), label = c("Phase 1 \n (Complete Reward Feedback)", "Phase 2 \n (Sensory + No Reward Feedback)") , 
           color="Black", size=3 , angle=0, fontface="bold") +
  theme_classic() 
 
uncertaintyPlot


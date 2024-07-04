#List of reusable functions

######## Create reward environments ################
createEnv <- function(nSubjects,nRep,nTrials,df){
  df.full <- do.call("rbind",replicate(nRep,df,simplify=F))
  nEnv <- length(unique(df$Environment))
  nRounds <- nEnv * nRep
  nC <- length(unique(df$Choice))
  rewards <- array(0, dim=c(nSubjects,nRounds,nTrials,nC))
  stim <- array(0, dim=c(nSubjects,nRounds,nTrials,nC))
  Rmean <- aperm(array(unlist(df.full$Rmean),c(nC,nRounds)))
  Rnoise <- aperm(array(unlist(df.full$Rvar),c(nC,nRounds)))
  stimLearn <- aperm(array(unlist(df.full$stimLearn),c(nC,nRounds)))
  
  for (s in 1:nSubjects){
    ls.rewards <- list()
    ls.stim <- list()
    for (i in 1:nrow(df.full)){
      ls.rewards <- c(ls.rewards,rnorm(nTrials,df.full$Rmean[i],sqrt(df.full$Rvar[i]))) #vector of size nTrials
      ls.stim <- c(ls.stim,(runif(nTrials) <= df.full$prob[i]))
    }
    #Need to reshape such that [nRounds,nTrials,nChoices]
    mat.rewards <- aperm(array(unlist(ls.rewards),c(nTrials,nC,nRounds)),c(3,1,2)) #must unlist!
    mat.stim <-aperm(array(unlist(ls.stim),c(nTrials,nC,nRounds)),c(3,1,2))
    rewards[s,,,] <- mat.rewards
    stim[s,,,] <- mat.stim
  }
  return(list(rewards,stim,Rmean,Rnoise,stimLearn[,1]))
}

######## Create reward envinronment for Exp 2 #############
# Adds a random shift in rewards to maintain differences 

createEnv2 <- function(nSubjects,nRep,nTrials,df){
  df.full <- do.call("rbind",replicate(nRep,df,simplify=F))
  nEnv <- length(unique(df$Environment))
  nRounds <- nEnv * nRep
  nC <- length(unique(df$Choice))
  rewards <- array(0, dim=c(nSubjects,nRounds,nTrials,nC))
  stim <- array(0, dim=c(nSubjects,nRounds,nTrials,nC))
  Rmean <- aperm(array(unlist(df.full$Rmean),c(nC,nRounds)))
  Rnoise <- aperm(array(unlist(df.full$Rvar),c(nC,nRounds)))
  stimLearn <- aperm(array(unlist(df.full$stimLearn),c(nC,nRounds)))
  
  for (s in 1:nSubjects){
    ls.rewards <- list()
    ls.stim <- list()
    for (i in 1:nrow(df.full)){
      rd.rewards <- rnorm(nTrials,df.full$Rmean[i],sqrt(df.full$Rvar[i])) + round(rnorm(1,0,1))
      ls.rewards <- c(ls.rewards,rd.rewards) #vector of size nTrials
      ls.stim <- c(ls.stim,(runif(nTrials) <= df.full$prob[i]))
    }
    #Need to reshape such that [nRounds,nTrials,nChoices]
    mat.rewards <- aperm(array(unlist(ls.rewards),c(nTrials,nC,nRounds)),c(3,1,2)) #must unlist!
    mat.stim <-aperm(array(unlist(ls.stim),c(nTrials,nC,nRounds)),c(3,1,2))
    rewards[s,,,] <- mat.rewards
    stim[s,,,] <- mat.stim
  }
  return(list(rewards,stim,Rmean,Rnoise,stimLearn[,1]))
}
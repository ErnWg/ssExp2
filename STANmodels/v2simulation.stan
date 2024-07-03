data {
  //Parameters for experiment 
  int<lower=0> nSub;//number of subjects
  int<lower=0> nRd;//number of rounds 
  int<lower=0> nP1;//number of phase 1 trials 
  int<lower=0> nP2;//number of phase 2 trials
  int<lower=0> nT;//total number of trials nP1 + nP2
  int<lower=0> nB;//number of objects
  real<lower=0> noise;// sigma of reward distribution 
  int<lower=0> opt1[nSub,nRd,nT];
  int<lower=0> opt2[nSub,nRd,nT];
  int<lower=0> stim[nSub,nRd,nT,nB];
  real rewards[nSub,nRd,nT,nB];
  
  //Parameters for model init
  real<lower=0> initR0;

  //Parameters for learning Rule 
  //Kalman filter is optimal so there is only the learning rate for sensory probabilities 
  real<lower=0> dA[nSub]; //learning rate for alpha counts
  real<lower=0> dB[nSub]; //learning rate for beta counts
  
  //Parameters for Choice Rule 
  //
  real theta[nSub]; //sensory value driven
  real omega[nSub]; //sensory uncertainty driven
  real tau[nSub]; //inverse temp for sensory learning
  
}

transformed data {
  vector[nB] initS0;
  vector[nB] noiseVec;
  real initU0;
  
  initS0 = rep_vector(1,nB);
  initU0 = sqrt(noise^2 * 20);
  noiseVec = rep_vector(noise,nB);
}

parameters {
}

transformed parameters{
}

model {
}

generated quantities {
  // For PPC
  real y_pred[nSub,nRd,nT];
  //For behaviour
  int chosenObj[nSub,nRd,nT]; //chosen object on every trial -- mapping from y_pred to object identity
  real Score[nSub,nRd,nT]; //points received
  int stimDeliver[nSub,nRd,nT];
  
  
  //For learning curves
  real Rmu_ts[nSub,nRd,nT,nB];
  real Rsig_ts[nSub,nRd,nT,nB];
  real Smu_ts[nSub,nRd,nT,nB];
  real Ssig_ts[nSub,nRd,nT,nB];
  
  
  // Set all PPs to -1
  for (s in 1:nSub){
    for (rd in 1:nRd){
      for (t in 1:nT){
        y_pred[s,rd,t] = -1;
      }
    }
  }
  
   {// local section 
    for (s in 1:nSub){
      //Initialise Values
      vector[2] Q; // Q[1] = utility of opt1, Q[2] = utility of opt2
      vector[nB] A; //alpha counts
      vector[nB] B; //beta counts
      vector[nB] Rmu; //Expected reward mean 
      vector[nB] Rsig; //Exoected reward uncertainty
      vector[nB] Smu; //Expected stim probability
      vector[nB] Ssig; //Expect stim uncertainty
      vector[nB] feedback;
      vector[nB] kgain;
      vector[nB] wB; //object weight. 1 = present at trial t
      
      
      for (rd in 1:nRd){
      //Initialise values
        //Reward Learning 
        Rmu = rep_vector(initR0, nB);
        Rsig = rep_vector(initU0, nB);
        //Sensory learning 
        A = initS0;
        B = initS0;
    
        for (t in 1:nT){
          
          feedback = rep_vector(0,nB);
          wB = rep_vector(0,nB);
          //Compute features
          for (b in 1:nB){
            //Reward learning curves 
            Rmu_ts[s,rd,t,b] = Rmu[b];
            Rsig_ts[s,rd,t,b] = Rsig[b];
            //Sensory probabilities learning curves
            Smu[b] = A[b]/(A[b] + B[b]);
            Ssig[b] = sqrt((A[b] * B[b]) / (((A[b]+B[b])^2)*(A[b]+B[b]+1)));
            Smu_ts[s,rd,t,b] = Smu[b];
            Ssig_ts[s,rd,t,b] = Ssig[b];
          }
          
          //Model phase 1
          if (t <= nP1){
            Q[1] = Rmu[opt1[s,rd,t]];
            Q[2] = Rmu[opt2[s,rd,t]];
            //generate posterior prediction for PPC 
            y_pred[s,rd,t] = categorical_rng(softmax(tau[s] * Q));
            
            //for behaviour
            if (y_pred[s,rd,t] == 1){chosenObj[s,rd,t] = opt1[s,rd,t];}
            else if (y_pred[s,rd,t] == 2){chosenObj[s,rd,t] = opt2[s,rd,t];} 
            
            Score[s,rd,t] = rewards[s,rd,t,chosenObj[s,rd,t]]; 
            
            //Update reward
            feedback[opt1[s,rd,t]] = rewards[s,rd,t,opt1[s,rd,t]]; //Update feedback vector
            feedback[opt2[s,rd,t]] = rewards[s,rd,t,opt2[s,rd,t]];
            wB[opt1[s,rd,t]] = 1; //Update Object weight
            wB[opt2[s,rd,t]] = 1;
          
            kgain = wB .* ((Rsig .*Rsig) ./ ((Rsig .*Rsig) + (noiseVec .*noiseVec)));//kalman gain 
            Rmu = Rmu + kgain .* (feedback - Rmu);
            Rsig = sqrt((1-kgain) .* (Rsig .* Rsig));
            
            stimDeliver[s,rd,t] = -1;
          }
          
          // //Model phase 2
          else if (t > nP1){
            Q[1] = Rmu[opt1[s,rd,t]] + theta[s]*Smu[opt1[s,rd,t]] + omega[s]*Ssig[opt1[s,rd,t]];
            Q[2] = Rmu[opt2[s,rd,t]] + theta[s]*Smu[opt2[s,rd,t]] + omega[s]*Ssig[opt2[s,rd,t]];

            //generate posterior prediction for PPC
            y_pred[s,rd,t] = categorical_rng(softmax(tau[s] * Q));

            //for behaviour
            if (y_pred[s,rd,t] == 1){chosenObj[s,rd,t] = opt1[s,rd,t];}
            else if (y_pred[s,rd,t] == 2){chosenObj[s,rd,t] = opt2[s,rd,t];}

            Score[s,rd,t] = rewards[s,rd,t,chosenObj[s,rd,t]];

            //Update reward
            feedback[opt1[s,rd,t]] = rewards[s,rd,t,opt1[s,rd,t]]; //Update feedback vector
            feedback[opt2[s,rd,t]] = rewards[s,rd,t,opt2[s,rd,t]];
            wB[opt1[s,rd,t]] = 1; //Update Object weight
            wB[opt2[s,rd,t]] = 1;

            kgain = wB .* ((Rsig .*Rsig) ./ ((Rsig .*Rsig) + (noiseVec .*noiseVec)));//kalman gain
            Rmu = Rmu + kgain .* (feedback - Rmu);
            Rsig = sqrt((1-kgain) .* (Rsig .* Rsig));

            //Update beta distributions

            if (stim[s,rd,t,chosenObj[s,rd,t]] == 1){
              A[chosenObj[s,rd,t]] += dA[s];
              stimDeliver[s,rd,t] = 1;
              }
            else {
              B[chosenObj[s,rd,t]] += dB[s];
              stimDeliver[s,rd,t] = 0;
              }
          }
        }
      }
    }
  }
}

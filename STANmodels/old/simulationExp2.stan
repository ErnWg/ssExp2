data {
  //Parameters for environment 
  int<lower=0> nSub; //number of subjects
  int<lower=0> nRd; //number of rounds
  int<lower=0> nT; //number of Trials in each round 
  int<lower=0> nC;//number of choices; bandits = 4
  int<lower=0> stim[nSub,nRd,nT,nC]; //stimulation delivered? 1 = yes; 0 = no
  real R[nSub,nRd,nT,nC]; //reward payoffs these should be normalised!
  real<lower=0> initS0; //initialistaion of beta distribution for sensory stimulation. We fix prior psuedocounts to 1
  
  //Parameters for learning Rule 
  //Kalman filter is optimal so there is only the learning rate for sensory probabilities 
  real<lower=0> dA[nSub]; //learning rate for alpha counts
  real<lower=0> dB[nSub]; //learning rate for beta counts
  
  //Parameters for Choice Rule 
  //
  real theta[nSub]; //sensory value driven
  real omega[nSub]; //sensory uncertainty driven
  real tauS[nSub]; //inverse temp for sensory learning
  real sticky[nSub]; //sticky bonus

}


parameters {
}

transformed parameters{
}

model {
}

generated quantities {
  // For PPC
  int y_pred[nSub,nRd,nT];//choices
  real rewards[nSub,nRd,nT];//rewards received
  int stimulation[nSub,nRd,nT]; //stimulation delivered?
  
  //For learning curves
  real Rmu_ts[nSub,nRd,nT,nC];
  real Rsig_ts[nSub,nRd,nT,nC];
  real Smu_ts[nSub,nRd,nT,nC];
  real Ssig_ts[nSub,nRd,nT,nC];

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
      vector[nC] Qstim;
  
      //Beta Distribution
      vector[nC] A; //alpha counts
      vector[nC] B; //beta counts
      vector[nC] Smu; //estimated probability of shock
      vector[nC] Ssig; //estimated uncertainty of shock
      //For stickiness
      vector[nC] prev_choice; //keep check of previous choice
      vector[nC] stickyb; // stickiness bonus
  
      for (rd in 1:nRd){
      //Initialise values
        //Sensory learning 
        A = rep_vector(initS0, nC);
        B = rep_vector(initS0, nC);
        prev_choice = rep_vector(0, nC); //initialise as all zeros
    
        for (t in 1:nT){
          for (c in 1:nC){
            //Reward learning curves (NO LEARNING!!!!)
            Rmu_ts[s,rd,t,c] = R[s,rd,t,c];
            Rsig_ts[s,rd,t,c] = 0;
            //Sensory probabilities learning curves
            Smu[c] = A[c]/(A[c] + B[c]);
            Ssig[c] = sqrt((A[c] * B[c]) / (((A[c]+B[c])^2)*(A[c]+B[c]+1)));
            Smu_ts[s,rd,t,c] = Smu[c];
            Ssig_ts[s,rd,t,c] = Ssig[c];
            
            //compute stickiness bonus
            stickyb = prev_choice * sticky[s];
            
            Qstim[c] = R[s,rd,t,c] + theta[s]*Smu[c] + omega[s]*Ssig[c];
          }
          
            
          y_pred[s,rd,t] = categorical_rng(softmax(tauS[s]*Qstim));
          
          rewards[s,rd,t] = R[s,rd,t,y_pred[s,rd,t]];
          stimulation[s,rd,t] = stim[s,rd,t,y_pred[s,rd,t]];
          //Update sensory learning
          if (stim[s,rd,t,y_pred[s,rd,t]] == 1) A[y_pred[s,rd,t]] += dA[s];
          else B[y_pred[s,rd,t]] += dB[s];
          //update previous choice
          prev_choice = rep_vector(0, nC); 
          prev_choice[y_pred[s,rd,t]] = 1; 
        }
      }
    }
  }
}

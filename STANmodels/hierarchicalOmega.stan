data {
  int<lower=0> nSub; //number of subjects
  int<lower=0> nRd; //number of rounds
  int<lower=0> nP1; //number of phase 1 trials
  int<lower=0> nP2; //number of phase 2 trials 
  int<lower=0> nT; //number of trials each round
  int<lower=0> nB; //number of bandits
  int<lower=0> opt1[nSub,nRd,nT]; //option 1 presented
  int<lower=0> opt2[nSub,nRd,nT]; //option 2 presented
  int<lower=0> choice[nSub,nRd,nT]; // 1 or 2 for left and right options
  //int<lower=0> stim[nSub,nRd,nT];
  //Predictions from learning models 
  real Rmu[nSub,nRd,nT,nB]; //Fixed payoffs for each
  //real<lower=0> Rsig[nSub,nRd,nT,nB];
  real<lower=0> Smu[nSub,nRd,nT,nB];
  real<lower=0> Ssig[nSub,nRd,nT,nB];
  
}

transformed data {
}

parameters {
  //Group level parameters
  //Hyper Priors for sensory learning
  //real mu_theta;
  //real<lower=0> sigma_theta;
  real mu_omega;
  real<lower=0> sigma_omega;
  real<lower=0> mu_tau;
  real<lower=0> sigma_tau;
  //real<lower=0> mu_dA;
  //real<lower=0> sigma_dA;
  //real<lower=0> mu_dB;
  //real<lower=0> sigma_dB;
  
  ////////////////////////////////////////////////////////
  //*Note that Tau here reflects inverse temperature
  //Subject level parameters
  //For sensory learning 
  //vector[nSub] theta_pr; //sensory driven choice
  vector[nSub] omega_pr; //sensory uncertainty coeff
  vector[nSub] tau_pr;
  //vector[nSub] dA_pr;
  //vector[nSub] dB_pr;
}

transformed parameters{

  //vector[nSub] theta;
  vector[nSub] omega;
  vector<lower=0,upper=20>[nSub] tau;
  //vector<lower=0,upper=1>[nSub] dA;
  //vector<lower=0,upper=1>[nSub] dB;
  
  for (s in 1:nSub){
    //theta[s] = mu_theta + sigma_theta * theta_pr[s];
    omega[s] = mu_omega + sigma_omega * omega_pr[s];
    tau[s] = Phi_approx(mu_tau + sigma_tau * tau_pr[s])*20;
    //dA[s] = Phi_approx(mu_dA + sigma_dA * dA_pr[s]);
    //dB[s] = Phi_approx(mu_dB + sigma_dB * dB_pr[s]);
  }
}

model {
  //Hyperparameters

  //Hyperparameters for sensory learning
  //mu_theta ~ normal(0,1);
  //sigma_theta ~ normal(0,1);
  mu_omega ~ normal(0,1);
  sigma_omega ~ normal(0,1);
  mu_tau ~ normal(0,1);
  sigma_tau ~ normal(0,1);
  //mu_dA ~ normal(0,1);
  //sigma_dA ~ normal(0,1);
  //mu_dB ~ normal(0,1);
  //sigma_dB ~ normal(0,1);
  
  
  //individual parameters
  //theta_pr ~ normal(0,1);
  omega_pr ~ normal(0,1);
  tau_pr ~ normal(0,1); //half normal
  //dA_pr ~ normal(0,1);
  //dB_pr ~ normal(0,1);
  
  for (s in 1:nSub){
    //Define variables
      vector[2] Q;

      for (rd in 1:nRd){

        for (t in 1:nT){
          if (choice[s,rd,t] != 0){
            
            if (t<=nP1){
              Q[1] = Rmu[s,rd,t,opt1[s,rd,t]]; 
              Q[2] = Rmu[s,rd,t,opt2[s,rd,t]];
              target +=  categorical_logit_lpmf(choice[s, rd, t] | tau[s] * Q);
            }
            
            else if (t>nP1){
              Q[1] = Rmu[s,rd,t,opt1[s,rd,t]] + omega[s]*Ssig[s,rd,t,opt1[s,rd,t]]; 
              Q[2] = Rmu[s,rd,t,opt2[s,rd,t]] + omega[s]*Ssig[s,rd,t,opt2[s,rd,t]];
              target +=  categorical_logit_lpmf(choice[s, rd, t] | tau[s] * Q);
            }
        }
      }
    }
  }
}

generated quantities {
  
  // Log-likelihood calculation
  real log_lik[nSub];
  
  //For group level parameters
  // For PPC
 // real y_pred[nSub, nT];

  // Set all PPs to -1
  //for (s in 1:nSub){
   // for (t in 1:nT){
   //   y_pred[s,t] = -1;
  //  }
  //}
  
  { // local section 
    for (s in 1:nSub){
    //Define variables
      vector[2] Q;
      
      log_lik[s] = 0;
      
      for (rd in 1:nRd){
        
        for (t in 1:nT){
          if (choice[s,rd,t] != 0){
            if (t<=nP1){
              Q[1] = Rmu[s,rd,t,opt1[s,rd,t]]; 
              Q[2] = Rmu[s,rd,t,opt2[s,rd,t]];
              log_lik[s] +=  categorical_logit_lpmf(choice[s, rd, t] | tau[s] * Q);
            }
            
            else if (t>nP1){
              Q[1] = Rmu[s,rd,t,opt1[s,rd,t]] + omega[s]*Ssig[s,rd,t,opt1[s,rd,t]]; 
              Q[2] = Rmu[s,rd,t,opt2[s,rd,t]] + omega[s]*Ssig[s,rd,t,opt2[s,rd,t]];
              log_lik[s] +=  categorical_logit_lpmf(choice[s, rd, t] | tau[s] * Q);
            }
          }
        }
      }
    }
  }
}

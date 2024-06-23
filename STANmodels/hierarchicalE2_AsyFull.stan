data {
  int<lower=0> nSub; //number of subjects
  int<lower=0> nRd; //number of rounds
  int<lower=0> nT; //number of trials each round
  int<lower=0> nC; //number of choices
  int<lower=0> choice[nSub,nRd,nT];
  int<lower=0> stim[nSub,nRd,nT];
  //Predictions from learning models 
  real Rmu[nSub,nRd,nT,nC]; //Fixed payoffs for each
  //real<lower=0> Rsig[nSub,nRd,nT,nC];
  //real<lower=0> Smu[nSub,nRd,nT,nC];
  //real<lower=0> Ssig[nSub,nRd,nT,nC];
  
}

transformed data {
  vector[nC] initS0; //initialise alpha & beta values for beta distribution
  
  initS0 = rep_vector(1,nC); //initialise all counts as 1
}

parameters {
  //Group level parameters
  //Hyper Priors for sensory learning
  real mu_theta;
  real<lower=0> sigma_theta;
  real mu_omega;
  real<lower=0> sigma_omega;
  real<lower=0> mu_tauS;
  real<lower=0> sigma_tauS;
  real<lower=0> mu_dA;
  real<lower=0> sigma_dA;
  real<lower=0> mu_dB;
  real<lower=0> sigma_dB;
  
  ////////////////////////////////////////////////////////
  //*Note that Tau here reflects inverse temperature
  //Subject level parameters
  //For sensory learning 
  vector[nSub] theta_pr; //sensory driven choice
  vector[nSub] omega_pr; //sensory uncertainty coeff
  vector[nSub] tauS_pr;
  vector[nSub] dA_pr;
  vector[nSub] dB_pr;
}

transformed parameters{

  vector[nSub] theta;
  vector[nSub] omega;
  vector<lower=0,upper=20>[nSub] tauS;
  vector<lower=0,upper=1>[nSub] dA;
  vector<lower=0,upper=1>[nSub] dB;
  
  for (s in 1:nSub){
    theta[s] = mu_theta + sigma_theta * theta_pr[s];
    omega[s] = mu_omega + sigma_omega * omega_pr[s];
    tauS[s] = Phi_approx(mu_tauS + sigma_tauS * tauS_pr[s])*20;
    dA[s] = Phi_approx(mu_dA + sigma_dA * dA_pr[s]);
    dB[s] = Phi_approx(mu_dB + sigma_dB * dB_pr[s]);
  }
}

model {
  //Hyperparameters

  //Hyperparameters for sensory learning
  mu_theta ~ normal(0,1);
  sigma_theta ~ normal(0,1);
  mu_omega ~ normal(0,1);
  sigma_omega ~ normal(0,1);
  mu_tauS ~ normal(0,1);
  sigma_tauS ~ normal(0,1);
  mu_dA ~ normal(0,1);
  sigma_dA ~ normal(0,1);
  mu_dB ~ normal(0,1);
  sigma_dB ~ normal(0,1);
  
  
  //individual parameters
  theta_pr ~ normal(0,1);
  omega_pr ~ normal(0,1);
  tauS_pr ~ normal(0,1); //half normal
  dA_pr ~ normal(0,1);
  dB_pr ~ normal(0,1);
  
  for (s in 1:nSub){
    //Define variables
      vector[nC] Qstim;
      //Beta Distribution
      vector[nC] A; //alpha counts
      vector[nC] B; //beta counts
      vector[nC] Smu; //estimated probability of shock
      vector[nC] Ssig; //estimated uncertainty of shock
      
      for (rd in 1:nRd){
        //initialise values in each round
        //Sensory learning 
        A = initS0;
        B = initS0;
        //Keep track of prev choice for stickiness
        //prev_choice = rep_vector(0, nC); //initialise as all zeros
        
        for (t in 1:nT){
          if (choice[s,rd,t] != 0){
            for (c in 1:nC){
            //Beta distribution summary stats
            Smu[c] = A[c]/(A[c] + B[c]);
            Ssig[c] = sqrt((A[c] * B[c]) / (((A[c]+B[c])^2)*(A[c]+B[c]+1)));
            Qstim[c] =  Rmu[s,rd,t,c] + theta[s]*Smu[c] + omega[s]*Ssig[c];
            }
            //compute stickiness bonus
            //stickyb = prev_choice * sticky[s];
      
            target +=  categorical_logit_lpmf(choice[s, rd, t] | tauS[s] * Qstim);
            
            //update
            //sensory learning 
            if (stim[s,rd,t] == 1) A[choice[s,rd,t]] += dA[s];
            else B[choice[s,rd,t]] += dB[s];
            //update previous choice
            //prev_choice = rep_vector(0, nC); 
            //prev_choice[choice[s,rd,t]] = 1; 
  
        }
      }
    }
  }
}

generated quantities {
  
  // Log-likelihood calculation
  real log_lik[nSub];
  
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
      vector[nC] Qstim;
      //Beta Distribution
      vector[nC] A; //alpha counts
      vector[nC] B; //beta counts
      vector[nC] Smu; //estimated probability of shock
      vector[nC] Ssig; //estimated uncertainty of shock
      
      log_lik[s] = 0;
      
      for (rd in 1:nRd){
        //initialise values in each round
        //Sensory learning 
        A = initS0;
        B = initS0;
        //Keep track of prev choice for stickiness
        //prev_choice = rep_vector(0, nC); //initialise as all zeros
        
        for (t in 1:nT){
          if (choice[s,rd,t] != 0){
            for (c in 1:nC){
            //Beta distribution summary stats
            Smu[c] = A[c]/(A[c] + B[c]);
            Ssig[c] = sqrt((A[c] * B[c]) / (((A[c]+B[c])^2)*(A[c]+B[c]+1)));
            Qstim[c] =  Rmu[s,rd,t,c] + theta[s]*Smu[c] + omega[s]*Ssig[c];
            }
            //compute stickiness bonus
            //stickyb = prev_choice * sticky[s];
      
            log_lik[s] +=  categorical_logit_lpmf(choice[s, rd, t] | tauS[s] * Qstim);
            
            //update
            //sensory learning 
            if (stim[s,rd,t] == 1) A[choice[s,rd,t]] += dA[s];
            else B[choice[s,rd,t]] += dB[s];
            //update previous choice
            //prev_choice = rep_vector(0, nC); 
            //prev_choice[choice[s,rd,t]] = 1; 
  
          }
        }
      }
    }
  }
}

functions{
  real logistic(real t, real r, real K, real N0){
  
  return(K/(1 + ((K-N0)/N0)*exp(-r*t) ));
  }
}

data {
  int N_ts;              // Total number of observation times t >= 0
  real obs_times[N_ts];  // Times of the observations
  real y_obs[N_ts];      // Data observed
  real<lower=0> N0;      // can't infer r and N0 simultaneously, except by using strong priors
}

parameters{
  real<lower=0> sigma;
  real<lower=0> K;
  real<lower=0> r;
}

model{
  
  r ~ normal(0, 1);  // if r is a parameter, should include a prior for r, otherwise model doesn't converge
  K ~ normal(0, 10);
  sigma ~ cauchy(0, 1);
  
  
  for (t in 1:N_ts){
    y_obs[t] ~ normal(logistic(obs_times[t], r, K, N0), sigma); //Likelihood
    }
}

generated quantities{
  real y_checks[N_ts];
  real y_checks_mean[N_ts];
  
  for (t in 1:N_ts){
    y_checks_mean[t] = logistic(obs_times[t], r, K, N0);
    y_checks[t] = normal_rng(y_checks_mean[t], sigma);
    }
  
}


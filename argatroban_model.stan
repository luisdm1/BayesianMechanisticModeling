functions {
  vector ode_model(vector times, real logcl, real logv, real dose){
    
    real tinf;
    real cl;
    real v;
    real t1;
    real t2;
    vector[num_elements(times)] f1;
    
    tinf = 240.0;
    cl = exp(logcl);
    v = exp(logv);
    for (time_ix in 1:num_elements(times)){
      if (times[time_ix] <= tinf)  {t1=1;}  else {t1=0;}
      t2 = tinf * (1 - t1) + t1*times[time_ix];
      f1[time_ix] = (dose/cl) * (1 - exp(-cl*t2/v)) * exp(-cl*(1 - t1)*(times[time_ix] - tinf)/v);
    }
    return(f1);
  }
}

data{
  int N_subj;              //Number of subjects
  int Tot_obs;             // total of observations
  int N_rec_psubj[N_subj]; // Number of records per subject
  int N_rand_eff;          // Number of random effects  
  int p;                   // number of population level parameters of \beta vector
  
  int indiv[Tot_obs];      //vector indiv ID for all subjects
  real dose[Tot_obs];      //vector dose records all subjects
  vector[Tot_obs] obs_times;    //vector times of records all subjects
  
  real conc[Tot_obs];      //vector concentration records all subjects
  
  //// Paramters of priors
  // shape parameter lkj_corr_cholesky
  real eta_ljk_shape;                   
  // Population parameters prior
  real beta_loc;        
  real beta_scale;
  // Error process correlation parameters
  real sigma_e_loc;  
  real sigma_e_scale;
  real lambda_loc; real lambda_scale;
}

parameters{
  vector<lower=0>[N_rand_eff] sigma_bi;     //subj sd, vector 2x1 of scales of rand effects prior 
  cholesky_factor_corr[N_rand_eff] L_bi; // //Prior correlations (LKJ prior on cholesky correlation factor)
    matrix[N_rand_eff, N_subj] z_bi;  // matrix of location of the rand effects prior
  
  //Observational model scale and correlation parameters
  real<lower=0> sigma_e;
  vector[p] beta_raw; // non-centered parametrization
  real lambda_raw;
}

transformed parameters{
  
  row_vector[N_subj] logcl;  //log-Cl param vector all subjects
  row_vector[N_subj] logv;   //log-V param vector all subjects
  vector[Tot_obs] m_ij;      // Mean (Non-linear) function
  matrix[N_rand_eff, N_subj] bi;  // matrix subject random intercept
  vector[Tot_obs] Sigma_Y_ij;     // measurement error process
  real<lower=0> lambda;
  vector[p] beta; //population parameters
  
  beta = beta_loc + beta_scale * beta_raw;
  
  lambda = exp(lambda_scale*lambda_raw + lambda_loc); 
  // random effects
  bi = diag_pre_multiply(sigma_bi, L_bi)*z_bi; // subj random effects
  
  
  logcl = beta[1] + bi[1,:];
  logv =  beta[2] + bi[2,:];
  { int obs_ix;
    obs_ix = 1;
    
    //Evaluate ode model (mean function) at times 1,...,n_i for each patient
    // with their corresponding dose and ODE parameters 
    for (subj_ix in 1:N_subj){
      m_ij[obs_ix : (obs_ix + N_rec_psubj[subj_ix] - 1)] = ode_model(segment(obs_times, obs_ix, N_rec_psubj[subj_ix]), logcl[subj_ix], logv[subj_ix], dose[obs_ix]); 
      obs_ix = obs_ix + N_rec_psubj[subj_ix];    }
  }
    
    { vector[Tot_obs] m_ij_squared;
      for (obs_ix in 1:Tot_obs){
        m_ij_squared[obs_ix] = m_ij[obs_ix]^(lambda);}
      Sigma_Y_ij = sigma_e * m_ij_squared;
    }
    
}
  
  model{
    // Priors random effects
    // eta=2.0 => Off-diag near zero => No prior information about correlation 
    // between rand intercepts for logCl and logV
    L_bi ~ lkj_corr_cholesky(eta_ljk_shape);  
    to_vector(z_bi) ~ normal(0, 1);
    sigma_bi ~ cauchy(0,2.5);
    
    beta_raw ~ normal(0, 1);    //Prior population parameters
    sigma_e ~ cauchy(sigma_e_loc, sigma_e_scale);  // Prior observational model
    lambda_raw ~ normal(0, 1);
    
    conc ~ normal(m_ij,Sigma_Y_ij);
  }
  
generated quantities{
  cov_matrix[N_rand_eff] G;  
  matrix[N_rand_eff, N_rand_eff] G_scales_cor;
  row_vector[N_subj] median_cl; 
  row_vector[N_subj] median_v;
  real cv_cl; real cv_v;
  
  //For posterior checks and posterior predictive dist:
  real conc_check [Tot_obs];  
  real conc_check_mean [Tot_obs];
  
  
  G = diag_matrix(sigma_bi)*(L_bi * L_bi') *diag_matrix(sigma_bi);
  G_scales_cor[1, 1] = sqrt(G[1, 1]);
  G_scales_cor[2, 2] = sqrt(G[2, 2]);
  G_scales_cor[1, 2] = G[1, 2]/ (G_scales_cor[1, 1]*G_scales_cor[2,2]);
  G_scales_cor[2, 1] = G_scales_cor[1, 2];
  
  //Median clearance and volume params for all the patients
  median_cl = exp(logcl)*1000;
  median_v = exp(logv)*1000;
  
  // coefficient of variation of clearance and volume for all patients
  cv_cl = sqrt(G[1,1])*100;
  cv_v  = sqrt(G[2,2])*100;
  
  for (obs_ix in 1:Tot_obs){
    conc_check[obs_ix] = normal_rng(m_ij[obs_ix], Sigma_Y_ij[obs_ix]);
    conc_check_mean[obs_ix] = m_ij[obs_ix];
  }
        
}
  

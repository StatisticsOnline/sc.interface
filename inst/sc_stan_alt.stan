functions {
  
  matrix expand_sum_to_zero(matrix M) {
    int R = rows(M);
    int C = cols(M) + 1;
    matrix[R, C] E;
    E[:,1:(C-1)] = M;
    for(r in 1:R) {
      E[r,C] = -1 * sum(M[r,:]);
    }
    return(E);
  }
  
  vector ar_process(vector innovations,
                    real autocor,
                    real end_scale, 
                    int integrate_process) {
    
    int T = num_elements(innovations);
    vector[T] process;
    vector[T] ar_process;
    real scale_inno = sqrt(1 - (autocor^2));
    
    ar_process[1] = scale_inno * innovations[1];
    for(t in 2:T) {
      ar_process[t] = (autocor * ar_process[t-1]) 
                      + (scale_inno * innovations[t]);
    }
    
    return end_scale * ar_process;
  }

}

data {
  real<lower=0> t_df;
  
  // Flags for model settings
  // ---------------------------------------------------------------------------
  int<lower=0, upper=1> include_intercepts;
  int<lower=0, upper=1> include_time_coefs;
  int<lower=0, upper=1> include_unit_coefs;
  int<lower=0, upper=1> zap;
  // ---------------------------------------------------------------------------
  
  // Sizes of quantities.
  // ---------------------------------------------------------------------------
  // Number of units
  int<lower=1> N_units; 
  // Number of time periods observed.
  int<lower=2> T_times;
  // Number of latent factors.
  int<lower=1> K_latent; 
  // Number of observed covariates.
  int<lower=0> L_covars; 
  // ---------------------------------------------------------------------------
  
  // Observed variables.
  // ---------------------------------------------------------------------------
  // Observed potential outcomes of response variable.
  matrix[T_times, N_units] Y_obs; 
  
  // Covariates varying by both unit and time.
  array[N_units] matrix[T_times, L_covars] X; 
  
  // The first time at which each unit is treated. Values > T_times indicate
  // a unit is not treated.
  array[N_units] int treated_time; 
  // ---------------------------------------------------------------------------

  // Hyperparameters for priors.
  // ---------------------------------------------------------------------------
  // Prior scale hyperparameter for the causal effects.
  real<lower=0> causal_effects_prior_scale; 
  
  // Prior scale hyperparameter for per-unit random effects.
  real<lower=0> unit_intercept_prior_scale;
  
  // Prior scale hyperparameter for SD of reg. coefficients across units.
  vector<lower=0>[include_unit_coefs ? L_covars : 0] unit_coefs_sd_prior_scale;
  
  // Standard deviation of regression coefficients across time.
  vector<lower=0>[L_covars] time_coefs_sd;
  // ---------------------------------------------------------------------------
  
  // Upper and lower bounds for [0,1] bounded parameters.
  // ---------------------------------------------------------------------------
  // Bounds for autoregressive parameters.
  real<lower=0, upper=1> phi_latent_lb;
  real<lower=0, upper=2> phi_latent_ub;
  vector<lower=0, upper=1>[L_covars ? 1 : 0] phi_zeta_lb;
  vector<lower=0, upper=1>[L_covars ? 1 : 0] phi_zeta_ub;
  
  // Bounds for proportion of variance explained by latent process.
  real<lower=0, upper=1> lambda_ub;
  real<lower=0, upper=1> lambda_lb;
  // ---------------------------------------------------------------------------
}

transformed data {

  matrix[T_times, N_units] Y_obs_scaled;
  real Y_obs_overall_mean = mean(Y_obs);
  vector[N_units] Y_obs_means;
  vector[N_units] Y_obs_sds;
  vector[N_units] Y_obs_rms;

  for(n in 1:N_units) {
    Y_obs_sds[n] = sd(Y_obs[:,n]);
    Y_obs_rms[n] = sqrt(mean(Y_obs[:,n] .^ 2));
    Y_obs_means[n] = mean(Y_obs[:,n]);
  }

  real max_Y_sd =  max(Y_obs_sds);
  real max_Y_rms = max(Y_obs_rms);
  for(n in 1:N_units) {
    Y_obs_scaled[:,n] = (Y_obs[:,n] - Y_obs_overall_mean) / max_Y_sd; //
  }

  int num_treated_units = 0;
  for(n in 1:N_units) {
    if(treated_time[n] <= T_times) {
      num_treated_units += 1;
    }
  }

  array[num_treated_units] int treated_units;
  {
    int ti = 1;
    for(i in 1:N_units) {
      if(treated_time[i] <= T_times) {
        treated_units[ti] = i;
        ti += 1;
      }
    }
  }

  
  // Number of observed time periods for which each unit is treated.
  array[N_units] int num_treated_times;
  for(n in 1:N_units) {
    num_treated_times[n] = 1 + T_times - treated_time[n];
  }
  
  // For each unit n, number of total observed (period, unit) combinations for
  // which unit < i was treated.
  array[N_units] int num_treated_before = append_array(
    { 0 }, 
    cumulative_sum(num_treated_times[1:(N_units - 1)])
  );
  
  // Total number of observed (period, unit) combinations for which unit is 
  // treated in period.
  int total_treated = sum(num_treated_times);
  
  // Pre-computed zero intercept vector if include_intercept == 0.
  vector[include_intercepts ? 0 : N_units] zero_intercept;
  zero_intercept = rep_vector(0, N_units * (1 - include_intercepts));
  
  // Total number of nonzero loadings in loadings matrix
  int num_loadings = (N_units - K_latent) * K_latent + 
                       ((K_latent * (K_latent + 1)) %/% 2);
  
  // Indices for unpacking lower flattened random effects into loadings matrix
  array[K_latent + 1] int eff_ii;
  
  eff_ii[1] = 0;
  
  for(k in 1:K_latent) {
    eff_ii[k + 1] = eff_ii[k] + (N_units - 1) - (k - 1);
  }
  
  // Bounds for temporal autocorrelation of regression parameters
  vector[2] pz_bounds;
  if(L_covars) {
    pz_bounds[1] = phi_zeta_lb[1];
    pz_bounds[2] = phi_zeta_ub[1];
  } else {
    pz_bounds = [0, 1]';
  }
  
  // Earliest treatment time
  int earliest_treated = min(treated_time);
}

parameters {
  
  // Latent factor component.  
  // ---------------------------------------------------------------------------
  vector[num_loadings - K_latent] random_effects_flat_lower;
  vector<lower = 0>[K_latent] random_effects_flat_upper;
  // ---------------------------------------------------------------------------
  
  // Matrix of error terms for autoregressive latent factor process.
  // These are transformed into the latent process, factors.
  vector[T_times] factors_0_first;
  matrix[T_times, K_latent - 1] factors_0_rest;
  
  // Lag-1 autocorrelation for latent factor process.
  vector<lower=phi_latent_lb, upper=phi_latent_ub>[K_latent] factors_autocor;
  // ---------------------------------------------------------------------------
  
  // ---------------------------------------------------------------------------
  // Error terms for temporally varying regression coefficients.
  // These are transformed into autoregressive coefficients.
  matrix[T_times, L_covars] time_coefs_0;
  
  // Lag-1 autocorrelation between time-varying coefficients.
  vector<lower=pz_bounds[1], upper=pz_bounds[2]>[L_covars] time_coefs_autocor;
    
  // Scale-normalized unit-varying coefficients (ommitting the last).
  // These are constrained to sum to 0 in order to avoid nonidentification
  // with the temporally-varying coefficients.
  matrix[include_unit_coefs ? L_covars : 0, N_units - 1] unit_coefs_0;
  
  // Scale-normalized standard deviation of unit-varying coefficients,
  // describing how (dis)similar coefficients are across units.
  // vector<lower=0>[include_unit_coefs ? L_covars : 0] unit_coefs_sd_0;
  // ---------------------------------------------------------------------------
  
  // Scale-normalized unit-specific intercepts.
  vector[include_intercepts ? N_units : 0] unit_intercept_0;
  
  // Fraction of variance explained by latent factor process.
  real<lower=lambda_lb, upper=lambda_ub> frac_var_latent;

  // Scale-normalized causal effects.
  vector[total_treated] causal_effects_0;
}

transformed parameters {

  // Latent factor component.
  // ---------------------------------------------------------------------------

  
  matrix[N_units, K_latent] factor_loadings;
  
  for(k in 1:K_latent) {
    factor_loadings[1:(k-1), k] = rep_vector(0, k-1);
    factor_loadings[k, k] = random_effects_flat_upper[k];
    factor_loadings[(k+1):N_units, k] = 
      random_effects_flat_lower[(eff_ii[k] + 1):(eff_ii[k+1])];
  }
  
  matrix[N_units, K_latent] factor_loadings_unscaled;
  
  matrix[K_latent, N_units] factor_loadings_t = factor_loadings';
  // ---------------------------------------------------------------------------
  
  // Autoregressive latent factor process, constructed with stationary
  // SD = 1 and lag-1 autocorrelation = factors_autocor.
  matrix[T_times, K_latent] factors_0;
  factors_0 = append_col(factors_0_first, factors_0_rest);

  matrix[T_times, K_latent] factors;
  
  for(k in 1:K_latent) {
    factors[:, k] = ar_process(factors_0[:, k], factors_autocor[k], 1, 0);
  }
  
  // Per-unit latent mean process.
  matrix[T_times, N_units] latent_means = factors * factor_loadings_t;

  matrix[T_times, N_units] latent_component;
  for(n in 1:N_units) {
    latent_component[:, n] = (sqrt(frac_var_latent) * latent_means[:,n]);
  }
  // ---------------------------------------------------------------------------
  
  // Regression process component.
  // ---------------------------------------------------------------------------
  // Autoregressive time-varying coefficients.
  // These are constructed to have stationary SD = time_coefs_sd, and lag-1
  // autocorrelation = time_coefs_autocor.
  matrix[T_times, L_covars] time_coefs;

  for(l in 1:L_covars) {
    time_coefs[:, l] = ar_process(time_coefs_0[:, l],
                                  time_coefs_autocor[l],
                                  time_coefs_sd[l] / max_Y_sd, 0);
  }
  
  // Unit-varying coefficients.
  matrix[include_unit_coefs ? L_covars : 0, N_units] unit_coefs;
  
  // Scale-normalized unit-varying coefficients, constrained to sum to zero.
  matrix[include_unit_coefs ? L_covars : 0, N_units] unit_coefs_expanded;
  
  if(include_unit_coefs) {
    unit_coefs_expanded = expand_sum_to_zero(unit_coefs_0);
  
    unit_coefs = diag_pre_multiply(unit_coefs_sd_prior_scale / max_Y_sd,
                                   unit_coefs_expanded);
  } 
  
  matrix[T_times, N_units] regression_component;
  
  for(n in 1:N_units) {
    if(include_unit_coefs) {
      regression_component[,n] = rows_dot_product(X[n], time_coefs)
                                   + (X[n] * unit_coefs[:,n]);
    } else {
      regression_component[,n] = rows_dot_product(X[n], time_coefs);
    }
  }
  
  // ---------------------------------------------------------------------------
  
  // Unit-specific intercepts.
  vector[N_units] unit_intercept;
  if(include_intercepts) {
    unit_intercept = unit_intercept_prior_scale * unit_intercept_0; 
  } else {
    unit_intercept = zero_intercept;
  }
  
  // Causal effects.                         
  vector[total_treated] causal_effects = causal_effects_prior_scale 
                                           * causal_effects_0;

  // Vector of causal effects for unit n (identically 0 if unit n untreated).
  matrix[T_times, N_units] delta;
  for(n in 1:N_units) {
    delta[:,n] = rep_vector(0, T_times);
    delta[(treated_time[n]):T_times, n] = segment(causal_effects, 
                                                  min({num_treated_before[n]+1,
                                                       total_treated}), 
                                                  num_treated_times[n]); 
  }
}

model {
  
  // Per-unit exchangeable error likelihoods.
  // ---------------------------------------------------------------------------
  for(n in 1:N_units) {                                                
    // Vector of spillover effect for unit n
    
    // Mean vector for unit n, combining all effects except exchangeable error.
    vector[T_times] mean_n = unit_intercept[n]
                               + delta[:,n]
                               + regression_component[:,n]
                               + latent_component[:,n];
    
    // Standard deviation of exchangeable error for unit n.
    real sd_n = sqrt(1 - frac_var_latent);
    
    // Error likelihood for unit n.
    if(t_df > 0) {
      Y_obs_scaled[:,n] ~ student_t(t_df, mean_n, sd_n);
    } else {
      Y_obs_scaled[:,n] ~ normal(mean_n, sd_n);
    }
    
  }
  // ---------------------------------------------------------------------------
  
  // (Roughly) unit priors for scale-normalized quantities
  // ---------------------------------------------------------------------------
  if(t_df > 0) {
    to_vector(factors_0[1, :]) ~ student_t(t_df, 0, 1);
    to_vector(factors_0[2:T_times, :]) ~ student_t(t_df, 0, 1);
  } else {
    to_vector(factors_0[1, :]) ~ normal(0, 1);
    to_vector(factors_0[2:T_times, :]) ~ normal(0, 1);
  }
  
  // This prior ensures unit normal marginal priors on the components of
  // unit_coefs_expanded[l,:]. No Jacobian adjustment is needed since the
  // transformation is linear.
  if(include_unit_coefs) {
    to_vector(unit_coefs_expanded) ~ normal(0, inv(sqrt(1 - inv(N_units))));
  }
  
  to_vector(time_coefs_0) ~ normal(0, 1);
  
  unit_intercept_0 ~ normal(0, 1);
  
  causal_effects_0 ~ normal(0, 1);

  // ---------------------------------------------------------------------------
  
  for(n in 1:N_units) {
    if(n <= K_latent) {
      factor_loadings[n, 1:n] ~ normal(0, 2 / sqrt(n));
    } else {
      factor_loadings[n, 1:K_latent] ~ normal(0, 2 / sqrt(K_latent));
    }
  }
  
}

generated quantities {
  
  matrix[N_units, N_units] latent_cors = tcrossprod(factor_loadings);
  
  // Means for each unit and period.
  matrix[T_times, N_units] latent_trends;
  for(n in 1:N_units) {
    latent_trends[:,n] = unit_intercept[n]
                           + regression_component[,n]
                           + (sqrt(frac_var_latent) 
                                            * latent_means[:,n]);
    latent_trends[:,n] = (max_Y_sd * latent_trends[:,n]) + Y_obs_overall_mean;
  }

  matrix[T_times, num_treated_units] treated_control;
  treated_control = latent_trends[:, treated_units];

  vector[total_treated] causal_effects_scaled = max_Y_sd 
                                                * causal_effects;

  matrix[T_times, num_treated_units] effs = max_Y_sd * delta[:, treated_units];
                  
}

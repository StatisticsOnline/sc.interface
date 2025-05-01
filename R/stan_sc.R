stan_model <- cmdstanr::cmdstan_model(
  stan_file = system.file("sc_stan_alt.stan", package = "sc.interface"),
  stanc_options = list("O1")
)

first_one <- function(xs) {
  return(1 + sum(1 - xs))
}

format_data_to_stan <- function(
  data, response, time, unit, covars, treated_index, scale = FALSE
) {
  if (!is.data.frame(data)) {
    stop("Data must be provided in long format as a data frame.")
  }

  tt <- data |>
    dplyr::select(unit, treated_index) |>
    dplyr::group_by(dplyr::across(unit)) |>
    dplyr::summarize(
      dplyr::across(treated_index, first_one)
    )
  
  tt <- rev(tt[, treated_index, drop = TRUE])

  response_matrix <- as.matrix(
    data |>
      dplyr::select(response, time, unit) |>
      tidyr::pivot_wider(
        values_from = response,
        names_from = unit
      ) |>
      dplyr::select(!time) |>
      dplyr::select(
        sort(tidyselect::peek_vars(), decreasing = TRUE)
      )
  )

  unit_names <- colnames(response_matrix)

  response_matrix_scaled <- response_matrix - mean(response_matrix)
  response_matrix_scaled <- response_matrix_scaled / max(apply(response_matrix, 2, sd))

  num_covars <- length(covars)
  covar_array <- array(dim = c(ncol(response_matrix), nrow(response_matrix), num_covars))
  if (num_covars > 0) {
    for(c in 1:num_covars) {
      covar_mat <- t(as.matrix(
        data |>
          dplyr::select(covars[c], time, unit) |>
          tidyr::pivot_wider(
            values_from = covars[c],
            names_from = unit
          ) |>
          dplyr::select(!time) |>
          dplyr::select(
            sort(tidyselect::peek_vars(), decreasing = TRUE)
          )
      ))
      covar_array[, , c] <- covar_mat
    }
  }

  return(list(
    stan_data = list(
      Y_obs = response_matrix,
      X = covar_array,
      treated_time = tt,
      Ysc = response_matrix_scaled
    ),
    meta_data = list(
      unit_names = unit_names,
      treated_indices = which(tt < max(tt))
    )
  ))
}

fit_stan_model <- function(data_long,
                           num_latent,
                           response,
                           time,
                           unit,
                           covars,
                           treated_index,
                           thin = NULL,
                           output_dir = NULL,
                           sampler_options = NULL,
                           noreg = FALSE,
                           prev_fit = NULL,
                           include_intercepts = TRUE,
                           include_unit_coefs = TRUE,
                           integrate_factors = FALSE,
                           include_spillover = TRUE,
                           include_regression = TRUE,
                           sf_mean = 0.5,
                           sf_kappa = 1) {

  formatted_data <- format_data_to_stan(
    data_long, response, time, unit, covars, treated_index,
    scale = FALSE
  )

  meta_data <- formatted_data$meta_data

  data <- formatted_data$stan_data
  data <- c(data, list(
    T_times = nrow(data$Y_obs),
    N_units = ncol(data$Y_obs),
    L_covars = length(covars),
    Q_edges = 0
  ))

  data$include_intercepts <- as.numeric(include_intercepts)
  data$include_unit_coefs <- as.numeric(include_unit_coefs)
  data$include_spillover <- as.numeric(include_spillover)
  data$include_time_coefs <- 1
  data$integrated_factors <- as.numeric(integrate_factors)

  num_covars <- length(covars)
  avg_response <- colMeans(data$Y_obs)
  avg_covars <- as.data.frame(apply(data$X, c(1, 3), mean))

  if (include_regression) {
    ols_fit <- lm(avg_response ~ ., data = avg_covars)
    ols_coefs <- ols_fit$coefficients[2:(length(covars)+1)]
    coefs_scale_est <- 2 * max(abs(ols_coefs))

    if (!include_unit_coefs) {
      data$unit_coefs_sd_prior_scale <- numeric(length = 0)
    } else {
      data$unit_coefs_sd_prior_scale <- rep(coefs_scale_est, num_covars)
    }
  } else {
    coefs_scale_est <- NULL
  }

  data <- c(data, list(
    causal_effects_prior_scale = 2, # 0.4
    unit_intercept_prior_scale = 1,
    phi_latent_lb = 0, # 0.9
    phi_latent_ub = 0.99999, #0.99999
    phi_zeta_lb = 0,
    phi_zeta_ub = 1,
    overall_sd_prior_scale = 1,
    lambda_ub = 0.999, #0.999
    lambda_lb = 0.9,
    time_coefs_sd = rep(coefs_scale_est, num_covars),
    spillover_effects_prior_scale = rep(0, ncol(data$Y_obs)),
    K_latent = num_latent,
    zap = FALSE,
    T_pos = 0
  ))

  if(!include_regression) {
    data$X <- array(dim = c(data$N_units, data$T_times, 0))
    data$L_covars <- 0
    data$unit_coefs_sd_prior_scale <- numeric(length = 0)
    data$time_coefs_sd <- numeric(length = 0)
    data$phi_zeta_lb <- numeric(length = 0)
    data$phi_zeta_ub <- numeric(length = 0)
  }

  if (data$Q_edges > 0) {
    data$sf_mean <- sf_mean
    data$sf_kappa <- sf_kappa 
  } else {
    data$sf_mean <- numeric(length = 0)
    data$sf_kappa <- numeric(length = 0)
    data$edges_1 <- numeric(length = 0)
    data$edges_2 <- numeric(length = 0)
    data$spatial_scale_norm <- numeric(length = 0)
    data$spatial_scale <- numeric(length = 0)
    data$mean_scale <- numeric(length = 0)
  }

  if (noreg) {
    data$X <- array(dim = c(data$N_units, data$T_times, 0))
    data$L_covars <- 0
    data$unit_coefs_sd_prior_scale <- numeric(length = 0)
    data$time_coefs_sd <- numeric(length = 0)
    data$phi_zeta_lb <- numeric(length = 0)
    data$phi_zeta_ub <- numeric(length = 0)
  }

  f1 <- data$Ysc[,1]
  e1 <- c(f1[1], (f1[2:(data$T_times)] - 0.997*f1[1:(data$T_times - 1)])) / 0.07
  print(round(e1,3))
  e1_init <- e1 + rnorm(data$T_times, 0, 0.0001)
  e1_list <- list(
    factors_0_first = e1_init,
    # factors_autocor_0 = rep(0.99 * log(data$phi_latent_ub / (1 - data$phi_latent_ub)), data$K_latent),
    frac_var_latent = 0.98
  )
  e1_init_list <- list(e1_list, e1_list, e1_list, e1_list)

  data$t_df <- 30

  synth_fit <- stan_model$sample(
    data = data,
    iter_warmup = sampler_options$warm,
    iter_sampling = sampler_options$iter,
    adapt_delta = sampler_options$ad,
    max_treedepth = sampler_options$mt,
    init = e1_init_list,
    thin = thin,
    parallel_chains = floor(0.8 * parallel::detectCores()),
    output_dir = output_dir
  )

  meta_data$stan_data <- data

  synth_fit$cmdstan_diagnose()

  delta_causal <- posterior::as_draws_matrix(
    synth_fit$draws(variables = "causal_effects_scaled")
  )
  delta_causal_qs <- apply(
    delta_causal, 2, function(x) quantile(x, c(0.025, 0.25, 0.75, 0.975))
  )
  print("Estimated Causal Effects (Posterior Quartiles)")
  print(delta_causal_qs)

  return(list(
    posterior = synth_fit,
    meta = meta_data
  ))
}
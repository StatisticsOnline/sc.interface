model <- cmdstanr::cmdstan_model(
  stan_file = system.file("R", "sc_stan_alt.stan", package = "sc.interface"),
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
    Y_obs = response_matrix,
    X = covar_array,
    treated_time = rev(tt[, treated_index, drop = TRUE]),
    Ysc = response_matrix_scaled
  ))
}

fit_model <- function(data_long,
                      num_latent,
                      response,
                      time,
                      unit,
                      covars,
                      treated_index,
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

  data <- format_data_to_stan(
    data_long, response, time, unit, covars, treated_index,
    scale = FALSE
  )

  # print(data)
  # return()

  # treated_index <- numeric(length = length(treated_times))
  # for(i in 1:length(treated_times)) {
  #   treated_index[i] <- 1 + sum(treated_times[i] > unique(data_long[, time]))
  # }

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

  ols_fit <- lm(avg_response ~ ., data = avg_covars)
  ols_coefs <- ols_fit$coefficients[2:(length(covars)+1)]
  coefs_scale_est <- 2 * max(abs(ols_coefs))

  if (!include_unit_coefs) {
    data$unit_coefs_sd_prior_scale <- numeric(length = 0)
  } else {
    data$unit_coefs_sd_prior_scale <- rep(coefs_scale_est, num_covars)
  }

  # data <- c(data, list(
  #   causal_effects_prior_scale = sd(data$Y_obs),
  #   unit_intercept_prior_scale = 2 * max(apply(data$Y_obs, 2, sd)),
  #   phi_latent_lb = 0.9,
  #   phi_latent_ub = 0.9999,
  #   phi_zeta_lb = 0,
  #   phi_zeta_ub = 1,
  #   overall_sd_prior_scale = 2 * max(apply(data$Y_obs, 2, sd)),
  #   lambda_ub = 0.99, #0.99
  #   lambda_lb = 0,
  #   time_coefs_sd = rep(coefs_scale_est, num_covars),
  #   spillover_effects_prior_scale = rep(0, ncol(data$Y_obs)),
  #   K_latent = num_latent,
  #   zap = TRUE,
  #   T_pos = 0
  # ))

  data <- c(data, list(
    causal_effects_prior_scale = 0.4,
    unit_intercept_prior_scale = 0.2,
    phi_latent_lb = 0.9,
    phi_latent_ub = 0.99999, #0.99999
    phi_zeta_lb = 0,
    phi_zeta_ub = 1,
    overall_sd_prior_scale = 1,
    lambda_ub = 0.999, #0.99
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
  # e1 <- rep(1, data$T_times)
  print(round(e1,3))
  e1_init <- e1 # + rnorm(data$T_times, 0, 0.01)
  e1_list <- list(
    factors_0_first = e1_init,
    factors_autocor_0 = rep(0.95 * log(data$lambda_ub / (1 - data$lambda_ub)), data$K_latent),
    frac_var_latent = 0.998
  )
  e1_init_list <- list(e1_list, e1_list, e1_list, e1_list)

  data$frac_var_latent <- 0.999
  data$t_df <- 30

  print(apply(data$Ysc, 2, sd))

  synth_fit <- model$sample(
    data = data,
    iter_warmup = sampler_options$warm,
    iter_sampling = sampler_options$iter,
    adapt_delta = sampler_options$ad,
    max_treedepth = sampler_options$mt,
    init = e1_init_list,
    parallel_chains = floor(0.8 * parallel::detectCores()),
    # output_dir = "./cmdstan_files/"
  )
  synth_fit$cmdstan_diagnose()

  delta_causal <- posterior::as_draws_matrix(
    synth_fit$draws(variables = "causal_effects_scaled")
  )
  delta_causal_qs <- apply(
    delta_causal, 2, function(x) quantile(x, c(0.025, 0.25, 0.75, 0.975))
  )
  print("Estimated Causal Effects (Posterior Quartiles)")
  print(delta_causal_qs)

  return(synth_fit)
}
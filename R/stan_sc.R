# Compile Stan model
stan_model <- cmdstanr::cmdstan_model(
  stan_file = system.file("sc_stan.stan", package = "sc.interface"),
  stanc_options = list("O1")
)

# Utility function
first_one <- function(xs) {
  return(1 + sum(1 - xs))
}

#' Take data in long (data frame) format, and convert to list of
#' data variables for Stan model
#'
#' @param data A data frame containing data for fitting causal
#' latent factor model
#' @param response A string indicating the name of the column in
#' data containing the response/outcome variable (e.g. GDP)
#' @param time A string indicating the name of the column in
#' data containing the time variable (e.g. year)
#' @param unit A string indicating the name of the column in
#' data containing the unit variable (e.g. country)
#' @param covars A vector of strings indicating the columns
#' in data containing the covariate variables (e.g. inflation)
#' @param treated_index A string indicating the column in data
#' containing the treatement indicator (0 or 1)
#' @returns A list with components stan_data (which contains
#' data variables X, Y_obs, and treated_time formatted for
#' the Stan model) and meta_data (which contains the names of
#' the units from the data variable as well as a vector of
#' which unit's were treated for at least one time, by index).
format_data_to_stan <- function(
  data, response, time, unit, covars, treated_index,
  prop_treated, wishart_cov, wishart_df
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
  response_matrix_scaled <-
    response_matrix_scaled / max(apply(response_matrix, 2, sd))

  order_by_units <- function(names) {
    return(order(factor(names, levels = as.factor(unit_names))))
  }
  prop_treated <- prop_treated[order_by_units(names(prop_treated))]
  cov_names <- order_by_units(colnames(wishart_cov))
  wishart_cov <- wishart_cov[cov_names, cov_names, drop = FALSE]

  num_covars <- length(covars)
  covar_array <-
    array(dim = c(ncol(response_matrix), nrow(response_matrix), num_covars))
  if (num_covars > 0) {
    for (c in 1:num_covars) {
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

  if (!all(tt[tt < max(tt)] == min(tt))) {
    stop("Stan model requires all treated units to be treated at same time.")
  }

  return(list(
    stan_data = list(
      Y_obs = response_matrix,
      X = covar_array,
      treated_time = min(tt),
      treated_units = which(tt < max(tt)),
      N_treated = length(which(tt < max(tt))),
      Ysc = response_matrix_scaled,
      prop_treated = prop_treated,
      cov_eff_re_df = wishart_df,
      cov_eff_re_S = wishart_cov
    ),
    meta_data = list(
      unit_names = unit_names,
      treated_indices = which(tt < max(tt))
    )
  ))
}

#' Take data in long (data frame) format, and convert to list of
#' data variables for Stan model
#'
#' @param data_long A data frame containing data for fitting causal
#' latent factor model
#' @param num_latent An integer >= 1 indicating the number of latent
#' factors to include in the model
#' @param response A string indicating the name of the column in
#' data containing the response/outcome variable (e.g. GDP)
#' @param time A string indicating the name of the column in
#' data containing the time variable (e.g. year)
#' @param unit A string indicating the name of the column in
#' data containing the unit variable (e.g. country)
#' @param covars A vector of strings indicating the columns
#' in data containing the covariate variables (e.g. inflation)
#' @param treated_index A string indicating the column in data
#' containing the treatement indicator (0 or 1)
#' @param thin A number giving the amount of thinning to apply
#' to Stan's posterior draws. The default value (NULL) uses
#' Stan's default thinning (0)
#' @param output_dir Directory to store Stan output files.
#' The default value (NULL) uses Stan's default, which stores
#' output files in a temporary directory. If saving the ScStanFit
#' object (e.g. in an RData file), it is strongly recommended to
#' set this to a non-null value. Otherwise, methods called on the
#' ScStanFit object may attempt to read from the Stan output files,
#' which will usually be automatically deleted after some time if
#' this parameter is NULL.
#' @param sampler_options A list of options to pass to Stan's
#' sampler, including adapt_delta and max_treedepth, which may be
#' increased if divergent transitions or treedepth exceedences are
#' encountered, respectively. See Stan documentation for details.
#' @param include_intercepts A boolean indicating whether to include a
#' random intercept for each unit in the model.
#' @param include_unit_coefs A boolean indicating whether to include
#' a unit-varying components in the regression coefficients. (This parameter
#' is only used if include_regression is TRUE and covars has length > 0.)
#' @param include_regression A boolean indicating whether to include the
#' regression component in the model.
#' @returns An ScStanFit object containing posterior draws along and
#' basic summary and plotting functions.
#' @export
fit_stan_model <- function(data_long,
                           num_latent,
                           response,
                           time,
                           unit,
                           covars,
                           treated_index,
                           prop_treated,
                           wishart_df,
                           wishart_cov,
                           thin = NULL,
                           output_dir = NULL,
                           sampler_options = NULL,
                           include_intercepts = TRUE,
                           include_unit_coefs = TRUE,
                           include_regression = TRUE,
                           hierarchical_delta = TRUE,
                           include_ce_random_effects = FALSE) {

  formatted_data <- format_data_to_stan(
    data_long, response, time, unit, covars, treated_index,
    prop_treated, wishart_cov, wishart_df
  )

  meta_data <- formatted_data$meta_data

  data <- formatted_data$stan_data
  data <- c(data, list(
    T_times = nrow(data$Y_obs),
    N_units = ncol(data$Y_obs),
    L_covars = length(covars)
  ))

  data$include_intercepts <- as.numeric(include_intercepts)
  data$include_unit_coefs <- as.numeric(include_unit_coefs)
  data$include_time_coefs <- 1

  num_covars <- length(covars)
  avg_covars <- as.data.frame(apply(data$X, c(1, 3), mean))

  if (include_regression) {
    avg_response <- colMeans(data$Y_obs)
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
    causal_effects_prior_scale = 2,
    unit_intercept_prior_scale = 1,
    phi_latent_lb = 0,
    phi_latent_ub = 0.99999,
    phi_zeta_lb = 0,
    phi_zeta_ub = 1,
    overall_sd_prior_scale = 1,
    lambda_ub = 0.999,
    lambda_lb = 0,
    time_coefs_sd = rep(coefs_scale_est, num_covars),
    spillover_effects_prior_scale = rep(0, ncol(data$Y_obs)),
    K_latent = num_latent,
    zap = FALSE,
    T_pos = 0,
    hierarchical_delta = hierarchical_delta,
    include_ce_random_effects = include_ce_random_effects
  ))

  if (!include_regression) {
    data$X <- array(dim = c(data$N_units, data$T_times, 0))
    data$L_covars <- 0
    data$unit_coefs_sd_prior_scale <- numeric(length = 0)
    data$time_coefs_sd <- numeric(length = 0)
    data$phi_zeta_lb <- numeric(length = 0)
    data$phi_zeta_ub <- numeric(length = 0)
  }

  if (!include_regression) {
    data$X <- array(dim = c(data$N_units, data$T_times, 0))
    data$L_covars <- 0
    data$unit_coefs_sd_prior_scale <- numeric(length = 0)
    data$time_coefs_sd <- numeric(length = 0)
    data$phi_zeta_lb <- numeric(length = 0)
    data$phi_zeta_ub <- numeric(length = 0)
  }

  f1 <- data$Ysc[, 1]
  e1 <-
    c(f1[1], (f1[2:(data$T_times)] - 0.997 * f1[1:(data$T_times - 1)])) / 0.07
  e1_init <- e1 + rnorm(data$T_times, 0, 0.0001)
  e1_list <- list(
    factors_0_first = e1_init
  )
  e1_init_list <- list(e1_list, e1_list, e1_list, e1_list)

  data$t_df <- 0

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
    delta_causal, 2, function(x) quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975))
  )
  writeLines("Estimated Causal Effects (Posterior Quartiles)")
  print(delta_causal_qs)

  return(list(
    posterior = synth_fit,
    meta = meta_data
  ))
}
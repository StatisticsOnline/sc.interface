get_wg_long <- function() {
  data <- haven::read_dta(
    system.file("R", "repgermany.dta", package = "sc.interface")
  )

  # Average covariates over pre-treatment times
  data_pre <- data |>
    dplyr::group_by(country) |>
    dplyr::summarize(dplyr::across(
      c(trade, infrate, industry, schooling, invest80),
      ~ mean(.x[data$year < 1990 & data$year >= 1980], na.rm = TRUE)
    ), dplyr::across(c(gdp, year))) |>
    dplyr::arrange(dplyr::desc(country)) |>
    dplyr::mutate(D = as.numeric(year >= 1990 & country == "West Germany")) |>
    dplyr::ungroup()

  data_pre$time <- data_pre$year

  return(as.data.frame(data_pre))
}

#' @export 
run_wg_example <- function(latent = 5, inc_reg = TRUE) {
  wg_data_long <- get_wg_long()

  # gs_fit <- gsynth::gsynth(
  #   Y = "gdp",
  #   X = c(), 
  #   D = "D",
  #   index = c("country", "year"),
  #   data = wg_data_long,
  #   r = latent,
  #   force = "time"
  # )

  # bp_fit <- bpCausal::bpCausal(
  #   data = wg_data_long,
  #   index = c("country", "time"),
  #   Yname = "gdp",
  #   Dname = "D",
  #   Xname = c("trade", "infrate", "industry", "schooling", "invest80"),
  #   Aname = c("trade", "infrate", "industry", "schooling", "invest80"),
  #   Zname = c(),
  #   xlasso = 1,
  #   zlasso = 1,
  #   alasso = 1,
  #   flasso = 1,
  #   r = latent,
  #   burn = 5000,
  #   niter = 15000
  # )

  # wg_data_stan <- get_wg_stan(latent = latent)
  # st_fit <- fit_model(
  #   wg_data_stan,
  #   sampler_options = list(warm = 1000, iter = 2000, ad = 0.8, mt = 12),
  #   include_spillover = FALSE,
  #   integrate_factors = FALSE,
  #   include_intercepts = FALSE,
  #   include_unit_coefs = FALSE
  # )
  # max_year <- max(wg_data_long$year) + 1
  # num_countries <- length(unique(wg_data_long$country))

  st_fit <- fit_model(
    wg_data_long,
    treated_index = "D",
    num_latent = latent,
    response = "gdp",
    time = "year",
    unit = "country",
    covars = c("trade"), # c("trade", "infrate", "industry", "schooling", "invest80"),
    sampler_options = list(warm = 1000, iter = 3000, ad = 0.8, mt = 11),
    include_spillover = FALSE,
    integrate_factors = FALSE,
    include_intercepts = FALSE,
    include_unit_coefs = FALSE,
    include_regression = inc_reg
  )

  return(list(
    st_fit = st_fit
    # gs_fit = gs_fit,
    # bp_fit = bp_fit
  ))  
}

## DEPRECATED: Stan interface now expects long data
## in the same format as gsynth and bpCausal 

get_wg_stan <- function(
  latent, scale_X = FALSE
) {

  data <- haven::read_dta(
    system.file("R", "repgermany.dta", package = "sc.interface")
  )
  data_pre <- data |>
    subset(year < 1990 & year >= 1980) |>
    dplyr::group_by(country) |>
    dplyr::summarize(dplyr::across(
      c(trade, infrate, industry, schooling, invest80),
      ~ mean(.x, na.rm = TRUE)
    )) |>
    dplyr::arrange(dplyr::desc(country)) |>
    dplyr::select(!country)

  if (scale_X) {
    data_pre_matrix <- scale(as.matrix(data_pre))
  } else {
    data_pre_matrix <- as.matrix(data_pre)
  }

  years <- unique(data$year)
  T_times = length(years)
  N_units = length(unique(data$country))
  L_covars = ncol(data_pre)

  size_data <- list(
    T_times = T_times,
    N_units = N_units,
    L_covars = L_covars,
    Q_edges = 0
  )

  X <- array(dim = c(N_units, T_times, L_covars))
  for (t in 1:T_times) {
    X[, t, ] <- data_pre_matrix
  }

  Y_obs <- as.matrix(data |>
                     dplyr::select(gdp, year, country) |>
                     tidyr::pivot_wider(
                       values_from = gdp,
                       names_from = country
                     ) |>
                     dplyr::select(!year) |>
                     dplyr::select(
                      sort(tidyselect::peek_vars(), decreasing = TRUE)))

  treated_time <- sum(1990 >= unique(data$year))

  variable_data <- list(
    X = X,
    Y_obs = Y_obs,
    treated_time = c(treated_time, rep(T_times + 1, N_units - 1)),
    time_values = years
  )

  ols_fit <- lm(colMeans(Y_obs) ~ ., data = as.data.frame(data_pre_matrix))
  ols_coefs <- ols_fit$coefficients[2:(size_data$L_covars+1)]
  coefs_scale_est <- 2 * max(abs(ols_coefs))

  ## Values of spatial_scale and mean_scale chosen so that GM of marginal
  ## variances is approximately 1 and so that the GM of correlations between
  ## treated unit and neighbors is approximately 0.55.
  ## Run test_and_check_prior in order to test different hyperparameter choices.
  hyperparameter_data <- list(
    causal_effects_prior_scale = max(apply(Y_obs, 2, sd)),
    unit_intercept_prior_scale = 2 * max(apply(Y_obs, 2, sd)),
    phi_latent_lb = 0.9, #0.99,
    phi_latent_ub = 0.9999, #0.98,
    phi_zeta_lb = 0, #0.99,
    phi_zeta_ub = 1, #0.98
    overall_sd_prior_scale = 2 * max(apply(Y_obs, 2, sd)),
    lambda_ub = 0.99,
    lambda_lb = 0.99,
    time_coefs_sd = rep(coefs_scale_est, L_covars),
    unit_coefs_sd_prior_scale = rep(coefs_scale_est, L_covars),
    spillover_effects_prior_scale = rep(0, ncol(Y_obs)),
    K_latent = latent,
    T_pos = 0
  )

  hyperparameter_data <- c(hyperparameter_data, list(
    edges_1 = numeric(length = 0),
    edges_2 = numeric(length = 0),
    spatial_scale_norm = numeric(length = 0)
  ))

  hyperparameter_data <- c(hyperparameter_data, list(
    spatial_scale = numeric(length = 0),
    mean_scale = numeric(length = 0)
  ))

  model_data <- c(size_data, variable_data, hyperparameter_data)
  return(model_data)
}
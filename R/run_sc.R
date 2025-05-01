#' Generate synthetic control estimates from specified backends
#' 
#' `run_sc()` estimates causal effects for treated units in
#' panel data using three variants of the synthetic control
#' idea: (i) the generalized synthetic control estimator 
#' (implemented by the gsynth package), (ii) the DF-LFM
#' Bayesian latent factor model, and (iii) a custom variant of
#' the DF-LFM model implemented in Stan.
#' 
#' @param data The input panel data in long format, which 
#'  must include, at minimum, columns for the response variable, 
#'  unit variable, time variable, and treatement indicator variable.
#' @param response Name of column containing the response variable.
#' @param panel_index Vector of two strings containing (i) the name
#'  of the column containing the unit variable, and (ii) the name
#'  of the column containing the time variable.
#' @param covars Vector of names of columns corresponding to regression
#'  covariate variables which are adjusted for.
#' @param num_latent The number of latent factors to fit (used for both
#'  bpCausal and Stan models). For bpCausal, some number of factors less
#'  than num_latent may be selected by the model's specification search.
#' @param bpc_warmup Number of warmup / burn-in iterations to run for
#'  bpCausal's sampler.
#' @param bpc_iter Number of post-warmup iterations to run for bpCausal's
#'  sampler.
#' @param stan_warmup Number of warmup / burn-in iterations to run for
#'  Stan's sampler.
#' @param stan_iter Number of post-warmup iterations to run for Stan's
#'  sampler.
#' @param stan_ad Adapt-delta hyperparameter for Stan's dHMC sampler.
#' @param stan_mt Max-treedepth hyperparameter for Stan's dHMC sampler.
#' @param include_covars Vector of backends for which covariate variables
#'  should be adjusted for.
#' @param backends Vector of backends for which to perform inference.
#' 
#' @returns List with entries `gs_fit`, `bp_fit`, and `st_fit` containing
#'  inference output from gsynth, bpCausal, and Stan methods respectively,
#'  or `NULL` if corresponding method not included in `backends`.
#' 
#' @export
run_sc <- function(
  data,
  response,
  panel_index,
  covars,
  num_latent,
  treated_index,
  bpc_warmup = 5000,
  bpc_iter = 15000,
  stan_warmup = 1000,
  stan_iter = 2000,
  stan_ad = 0.8,
  stan_mt = 11,
  include_covars = c("gs", "bp", "st"),
  backends = c("gs", "bp", "st")
) {

  data[, panel_index[1]] <- as.factor(data[, panel_index[1]])

  do_st <- "st" %in% backends
  do_gs <- "gs" %in% backends
  do_bp <- "bp" %in% backends

  if (do_gs) {
    gs_inf <- gsynth::gsynth(
      Y = response,
      X = covars,
      D = treated_index,
      index = panel_index,
      data = data,
      r = num_latent,
      force = "time",
      se = TRUE,
      CV = TRUE,
      inference = "parametric",
      nboots = 2000
    )
  } else {
    gs_inf <- NULL
  }

  if (do_bp) {

    bp_inf <- bpCausal::bpCausal(
      data = data,
      index = panel_index,
      Yname = response,
      Dname = treated_index,
      Xname = covars,
      Aname = covars,
      Zname = c(),
      xlasso = 1,
      zlasso = 1,
      alasso = 1,
      flasso = 1,
      r = num_latent,
      burn = bpc_warmup,
      niter = bpc_iter
    )
  } else {
    bp_inf <- NULL
  }

  if (do_st) {
    st_inf <- fit_stan_model(
      data,
      treated_index = treated_index,
      num_latent = num_latent,
      response = response,
      time = panel_index[2],
      unit = panel_index[1],
      covars = covars,
      sampler_options = list(
        warm = stan_warmup, 
        iter = stan_iter, 
        ad = stan_ad, 
        mt = stan_mt
      ),
      include_spillover = FALSE,
      integrate_factors = FALSE,
      include_intercepts = TRUE,
      include_unit_coefs = FALSE,
      include_regression = ("st" %in% include_covars),
      output_dir = "./"
    )
  } else {
    st_inf <- NULL
  }

  return(list(
    st_fit = st_inf,
    gs_fit = gs_inf,
    bp_fit = bp_inf
  ))  
}
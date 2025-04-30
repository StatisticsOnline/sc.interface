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
      include_regression = ("st" %in% include_covars)
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
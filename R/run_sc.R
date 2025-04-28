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
  stan_mt = 12,
  include_covars = c("gs", "bp", "st"),
  backends = c("gs", "bp", "st")
) {
  gs_fit <- gsynth::gsynth(
    Y = response,
    X = covars,
    D = treated_index,
    index = panel_index,
    data = data,
    r = num_latent,
    force = "time"
  )

  bp_fit <- bpCausal::bpCausal(
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

  st_fit <- fit_model(
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
    include_intercepts = FALSE,
    include_unit_coefs = FALSE
  )

  # print(wg_data_long[,c("infrate", "country", "year")])
  # print(wg_data_stan$X[,,1])
  # print(st_fit$X[,,1])
  return(list(
    st_fit = st_fit,
    gs_fit = gs_fit,
    bp_fit = bp_fit
  ))  
}
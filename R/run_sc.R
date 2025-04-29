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
  # data[data$country == "Australia", "D"] <- data[data$country == "West Germany", "D"]

  do_st <- "st" %in% backends
  do_gs <- "gs" %in% backends
  do_bp <- "bp" %in% backends

  if (do_gs) {
    gs_fit <- gsynth::gsynth(
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

    gs_inf <- list(
      Y.ct = gs_fit$Y.ct,
      Eff = gs_fit$est.ind[, 1, ],
      SE = gs_fit$est.ind[, 2, ],
      out = gs_fit
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
      include_intercepts = FALSE,
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

summarize_gs <- function(fit, stats = c("y.ct", "eff", "err"), format = "long") {

  stat_funs <- list(
    "y.ct" = function(fit) fit$Y.ct,
    "eff" = function(fit) abind::adrop(fit$est.ind[, 1, , drop = FALSE], drop = 2),
    "err" = function(fit) abind::adrop(fit$est.ind[, 2, , drop = FALSE], drop = 2)
  )

  statfs <- stat_funs[stats]
  stat_sums <- lapply(statfs, function(f) f(fit))
  sum <- do.call(function(...) abind::abind(..., along = 3), stat_sums)
  dimnames(sum)[[1]] <- seq(1, dim(sum)[1])

  if (format == "long") {
    sum <- plyr::adply(sum, c(1, 2), .id = c("unit", "time_index"))
  }

  return(sum)
}

summarize_bp <- function(fit, stats = c("y.ct", "eff", "err"), format = "long") {
  stat_funs <- list(
    y.ct = function(x, i) mean(x),
    err = function(x, i) sd(x),
    eff = function(x, i) as.numeric(fit$yo_t)[i] - mean(x)
  )

  apply_and_format <- function(stat, stat_name, wide = FALSE) {
    stats_long <- sapply(seq(1, nrow(fit$yct)), function(i) stat(fit$yct[i, ], i))
    stats_df <- data.frame(
      stats_long, fit$raw.id.tr, fit$time.tr
    )
    if(format == "long") {
      colnames(stats_df) <-  c(stat_name, "unit", "time_index")
      return(stats_df)
    } else if (format == "wide") {
      colnames(stats_df) <-  c("stat", "unit", "time_index")
      stats_mat <- as.matrix(stats_df |> tidyr::pivot_wider(
        names_from = unit,
        values_from = stat
      ) |> 
      dplyr::select(!time_index))
      return(stats_mat)
    } else {
      stop("Format misspecified, must be either 'wide' or 'long'.")
    }
  }

  sum <- NULL
  si <- 1

  for(stat in stats) {
    stat_sum <- apply_and_format(stat_funs[[stat]], stat, wide = (format == "wide"))
    if(si > 1 && format == "long") {
      sum <- merge(sum, stat_sum, sort = FALSE)
    } else if (si > 1 && format == "wide") {
      sum <- abind::abind(sum, stat_sum, along = 3)
    } else {
      sum <- stat_sum
    }
    si <- si + 1
  }

  if (format == "wide") {
    dimnames(sum)[[1]] <- seq(1, dim(sum)[1])
    dimnames(sum)[[3]] <- stats
  }

  return(sum)
}

summarize_stan <- function(fit, stats = c("y.ct", "eff", "err"), format = "long") {
  stan_draws <- posterior::merge_chains(fit$posterior$draws(
    c("effs", "treated_control")
  ))
  effs <- posterior::extract_variable_array(stan_draws, "effs")
  tc <- posterior::extract_variable_array(stan_draws, "treated_control")

  stat_funs <- list(
    "y.ct" = function() apply(tc, c(3, 2), mean),
    "eff" = function() apply(effs, c(3, 2), mean),
    "err" = function() apply(effs, c(3, 2), sd)
  )

  statfs <- stat_funs[stats]
  stat_sums <- lapply(statfs, function(f) f())
  sum <- do.call(function(...) abind::abind(..., along = 3), stat_sums)
  dimnames(sum)[[1]] <- seq(1, dim(sum)[1])
  dimnames(sum)[[2]] <- fit$meta$unit_names[fit$meta$treated_indices]

  if (format == "long") {
    sum <- plyr::adply(sum, c(1, 2), .id = c("time_index", "unit"))
  }

  return(sum)
}
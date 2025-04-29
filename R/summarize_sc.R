#' @export 
summarize <- function(fits, stats = c("y.ct", "eff", "err"), format = "long", level = NULL) {
  if (!is.null(level) && format != "long") {
    stop("Intervals only available in long format.")
  }
  
  fnames <- names(fits)
  sums <- lapply(fnames, function(fname) {
    fit <- fits[[fname]]
    if(is.null(fit)) {
      return(NULL);
    } else {
      if(fname == "st_fit") {
        sum <- summarize_stan(fit, stats = stats, format = format)
        if(!is.null(level)) {
          sum <- merge(sum, intervals_stan(fit, level = level))
        }
        return(sum)
      } else if (fname == "gs_fit") {
        sum <- summarize_gs(fit, stats = stats, format = format)
        if(!is.null(level)) {
          sum <- merge(sum, intervals_gs(fit, level = level))
        }
        return(sum)
      } else if (fname == "bp_fit") {
        sum <- summarize_bp(fit, stats = stats, format = format)
        if(!is.null(level)) {
          sum <- merge(sum, intervals_bp(fit, level = level))
        }
        return(sum)
      }
    }
  })
  names(sums) <- fnames

  if (format == "long") {
    if(!is.null(sums[["st_fit"]])) {
      sums[["st_fit"]]$method <- "stan"
    }
    if(!is.null(sums[["bp_fit"]])) {
      sums[["bp_fit"]]$method <- "bpcausal"
    }
    if(!is.null(sums[["gs_fit"]])) {
      sums[["gs_fit"]]$method <- "gsynth"
    }
    sums <- do.call(function(...) rbind(..., make.row.names = FALSE), sums)
  }
  return(sums)
}

summarize_gs <- function(
  fit,
  stats = c("y.ct", "eff", "err"),
  format = "long") {

  stat_funs <- list(
    "y.ct" = function(fit) fit$Y.ct,
    "eff" = function(fit) abind::adrop(fit$est.ind[, 1, , drop = FALSE], drop = 2),
    "err" = function(fit) abind::adrop(fit$est.ind[, 2, , drop = FALSE], drop = 2)
  )

  statfs <- stat_funs[stats]
  stat_sums <- lapply(statfs, function(f) f(fit))
  sum <- do.call(function(...) abind::abind(..., along = 3), stat_sums)
  dimnames(sum)[[1]] <- seq(1, dim(sum)[1])

  cns <- c("time_index", "unit", dimnames(sum)[[3]])

  if (format == "long") {
    sum <- plyr::adply(sum, c(1, 2))
    colnames(sum) <- cns
    sum$time_index <- as.numeric(sum$time_index)
  }

  return(sum)
}

intervals_gs <- function(fit, level) {
  qlb <- (1 - level) / 2
  qub <- 1 - qlb

  qln <- paste0(round(100 * qlb, 1), "%")
  qun <- paste0(round(100 * qub, 1), "%")

  ests <- abind::adrop(fit$est.ind[, 1, , drop = FALSE], drop = 2)
  ses <- abind::adrop(fit$est.ind[, 2, , drop = FALSE], drop = 2)

  data <- abind::abind(ests, ses, along = 3)
  dimnames(data)[[1]] <- seq(1, dim(data)[1])
  dimnames(data)[[3]] <- c("est", "se")

  data_long <- plyr::adply(data, c(1, 2), .id = c("time_index", "unit")) |>
    dplyr::mutate(
      {{qln}} := est + qnorm(qlb) * se,
      {{qun}} := est + qnorm(qub) * se
    ) |> 
    dplyr::select(!c(est, se))

  return(data_long)
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

intervals_bp <- function(fit, level) {
  qlb <- (1 - level) / 2
  qub <- 1 - qlb

  qln <- paste0(round(100 * qlb, 1), "%")
  qun <- paste0(round(100 * qub, 1), "%")

  qs <- apply(
    -1 * (fit$yct - as.numeric(fit$yo_t)), 1, 
    function(y) quantile(y, c(qlb, qub))
  )
  qs_df <- data.frame(
    qs[1, ], qs[2, ], fit$raw.id.tr, fit$time.tr
  )
  colnames(qs_df) <- c(qln, qun, "unit", "time_index")

  return(qs_df)
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

  cns <- c("time_index", "unit", dimnames(sum)[[3]])

  if (format == "long") {
    sum <- plyr::adply(sum, c(1, 2))
    colnames(sum) <- cns
    sum$time_index <- as.numeric(sum$time_index)
  }

  return(sum)
}

intervals_stan <- function(fit, level) {
  qlb <- (1 - level) / 2
  qub <- 1 - qlb

  stan_draws <- posterior::merge_chains(fit$posterior$draws("effs"))
  effs <- posterior::extract_variable_array(stan_draws, "effs")
  qs_wide <- apply(
    effs, c(3, 2),
    function(eff) quantile(eff, c(qlb, qub))
  )
  dimnames(qs_wide)[[3]] <- fit$meta$unit_names[fit$meta$treated_indices]
  qs <- plyr::adply(qs_wide, c(2, 3), .id = c("time_index", "unit"))


  return(qs)
}
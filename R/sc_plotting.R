plot_panel <- function(data, 
                       response,
                       unit,
                       time,
                       fit_data = NULL, 
                       chain = NULL,
                       fit_data2 = NULL, 
                       labels = NULL,
                       only_treated = FALSE,
                       which_units = NULL) {
  
  if(is.null(labels)) {
    labels <- list(x = "Time", y = "Effect")
  }

  plot_df <- data.frame(
    Y = data[, response],
    unit = data[, unit],
    times = data[, time]
  )
  
  if(!is.null(fit_data)) {
    if (is.null(chain)) {
      draws <- posterior::merge_chains(fit_data$draws())
      chain <- 1
    } else {
      draws <- fit_data$draws()
    }
    
    latent_means <- posterior::extract_variable_array(
      draws, variable = "latent_trends"
    )[, chain, , ]
    expected_latent_means <- colMeans(latent_means, dims = 1)
    plot_df$latent_trends = as.numeric(expected_latent_means)
  }
  
  if(!is.null(fit_data2)) {
    draws2 <- posterior::merge_chains(fit_data2$draws())
    latent_means2 <- posterior::extract_variable_array(draws2, variable = "latent_trends")[,1,,]
    expected_latent_means2 <- colMeans(latent_means2, dims = 1)
    plot_df$latent_trends2 = as.numeric(expected_latent_means2)
  }
  
  if(!is.null(which_units)) {
    plot_df <- plot_df %>% 
      filter(unit %in% unit_labels[which_units]) 
  } else {
    # which_units <- seq(1, data$N_units)
  }
  
  # tt <- inf_data$treated_time
  # tt[tt > inf_data$T_times] <- NA
  # tt <- tt - 1
  # # tt <- inf_data$time_values[tt]
  
  if(only_treated) {
    treated_units <- which(!is.na(tt))
    plot_df <- filter(plot_df, unit %in% treated_units)
    tt <- tt[treated_units]
  }

  data_plot <- ggplot2::ggplot(plot_df) +
    ggplot2::geom_line(ggplot2::aes(x = times, y = Y), linewidth = 1) +
    # ggplot2::geom_vline(
    #   data = data.frame(
    #     unit = which_units,
    #     tt = tt
    #   ),
    #   ggplot2::aes(xintercept = tt)
    # ) +
    ggplot2::xlab(labels$x) +
    ggplot2::ylab(labels$y) +
    ggplot2::theme_bw()
  
  if(!only_treated) {
    data_plot <- data_plot + 
      ggplot2::facet_wrap(ggplot2::vars(unit), ncol = 5, nrow = 4, scales = "free_y") +
      ggplot2::theme(strip.background = ggplot2::element_rect(fill="white"))
  }
  
  if(!is.null(fit_data)) {
    data_plot <- data_plot +
      ggplot2::geom_line(ggplot2::aes(x = times, y = latent_trends), linewidth = 1, color = 'blue')
  }
  
  if(!is.null(fit_data2)) {
    data_plot <- data_plot +
      ggplot2::geom_line(ggplot2::aes(x = times, y = latent_trends2), color = 'red')
  }
  
  return(data_plot)
}

plot_treated <- function(inf_data, 
                         fits = NULL,
                         fit_labels = NULL,
                         labels = NULL) {
  
  if(is.null(labels)) {
    labels <- list(x = "Year", y = "Effect")
  }
  
  t_index <- which(inf_data$treated_time < nrow(inf_data$Y_obs))
  
  if(length(t_index) > 1) {
    print("Too many treated units!")
    return()
  }
  
  plot_df <- data.frame(
    Observed = inf_data$Y_obs[, t_index], 
    times = inf_data$time_values
  )
  
  for(fi in 1:length(fits)) {
    fit_data <- fits[[fi]]
    draws <- posterior::merge_chains(fit_data$draws())
    latent_means <- posterior::extract_variable_array(draws, variable = "latent_trends")[,1,,t_index]
    expected_latent_means <- colMeans(latent_means)
    plot_df[fit_labels[fi]] <- as.numeric(expected_latent_means)
  }
  
  plot_df <- plot_df %>%
    dplyr::pivot_longer(
      !times,
      values_to = "Y",
      names_to = "Source"
    )
  
  tt <- inf_data$treated_time[t_index] - 1
  
  data_plot <- ggplot2::ggplot(plot_df) +
    ggplot2::geom_line(ggplot2::aes(x = times, y = Y, color = Source), linewidth = 1) +
    ggplot2::geom_vline(
      data = data.frame(
        tt = inf_data$time_values[tt]
      ),
      ggplot2::aes(xintercept = tt)
    ) +
    ggplot2::xlab(labels$x) +
    ggplot2::ylab(labels$y) +
    ggplot2::scale_color_manual(values = c("black", "blue", "red", "darkgreen", "gold")) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title=element_blank())
  
  return(data_plot)
}

plot_effect <- function(
    fit_data, cred_level1 = 0.5, cred_level2 = 0.9, time_values = NULL, labels = NULL
  ) {
  
  delta <- as_draws_matrix(fit_data$draws(variables = "causal_effects"))
  delta_mean <- colMeans(delta)
  
  if(is.null(labels)) {
    labels <- list(x = "Year", y = "Effect")
  }
  
  if(is.null(time_values)) {
    time_values <- seq(1, length(delta_mean))
  } else {
    time_values <- tail(time_values, length(delta_mean))
  }
  
  delta_qs1 <- apply(
    delta, 2, 
    function(x) quantile(x, c(0.5 * (1 - cred_level1), (0.5 * (1 + cred_level1))))
  )
  
  delta_qs2 <- apply(
    delta, 2, 
    function(x) quantile(x, c(0.5 * (1 - cred_level2), (0.5 * (1 + cred_level2))))
  )
  
  delta_df <- data.frame(
    delta_m = delta_mean,
    delta_ql = delta_qs1[1,],
    delta_qu = delta_qs1[2,],
    delta_qll = delta_qs2[1,],
    delta_quu = delta_qs2[2,],
    time = time_values
  )
  delta_plot <- ggplot(delta_df) +
    geom_ribbon(
      ggplot2::aes(x = time, ymin = delta_qll, ymax = delta_quu),
      fill = 'grey', alpha = 0.5) +
    geom_ribbon(
      ggplot2::aes(x = time, ymin = delta_ql, ymax = delta_qu),
      fill = 'grey', alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_line(ggplot2::aes(x = time, y = delta_m), linewidth = 1) +
    xlab(labels$x) +
    ylab(labels$y) +
    theme_bw()
  return(delta_plot)
}

plot_regression <- function(fit_data) {
  draws <- posterior::merge_chains(fit_data$draws())
  time_coefs <- extract_variable_array(draws, variable = "time_coefs")[,1,,]
  tc_qs <- apply(time_coefs, c(2,3), function(x) quantile(x, c(0.025, 0.5, 0.975)))
  
  tc_df <- adply(tc_qs, c(2,3))
  tc_df <- plyr::rename(tc_df, c(X1 = "Times", X2 = "Covs", "2.5%" = "cil", "97.5%" = "ciu", "50%" = "med"))
  
  tc_df$Times <- as.numeric(tc_df$Times)
  tc_df$Covs <- as.numeric(tc_df$Covs)
  
  tcplot <- ggplot(tc_df) +
    geom_line(ggplot2::aes(x = Times, y = med), linewidth = 1) +
    geom_ribbon(ggplot2::aes(x = Times, ymin = cil, ymax = ciu), alpha = 0.5) +
    facet_wrap(ggplot2::vars(Covs), nrow = 1) +
    theme_bw()
  
  return(tcplot)
}

plot_breakdown <- function(fit_data, unit = 1, time) {
  lc <- as.numeric(fit_data$draws(paste0("latent_component[",time,",",unit,"]")))
  rc <- as.numeric(fit_data$draws(paste0("regression_component[",time,",",unit,"]")))
  bd_plot <- ggplot(data.frame(lc = lc, rc = rc)) +
    geom_point(ggplot2::aes(x = lc, y = rc)) +
    xlab("Latent Component") +
    ylab("Regression Component") +
    theme_bw()
  return(bd_plot)
}

plot_breakdown_series <- function(fit_data, unit = 1, time_values) {
  draws <- posterior::merge_chains(fit_data$draws())
  
  rc <- extract_variable_array(draws, variable = "regression_component")[,1,,unit]
  rcm <- colMeans(rc)
  lc <- extract_variable_array(draws, variable = "latent_component")[,1,,unit]
  lcm <- colMeans(lc)
  
  n <- length(time_values)
  comp <- c(rep("Latent Component", n), rep("Regression Component", n))
             
  bd_plot <- ggplot(data.frame(
                      est = c(lcm, rcm), 
                      comp = comp, 
                      time = rep(time_values,2))) +
    geom_line(ggplot2::aes(x = time, y = est, color = comp), linewidth = 1) +
    scale_color_manual(values = c("blue", "red")) + 
    xlab("Year") +
    ylab("Estimated Effect") +
    theme_bw() +
    theme(legend.title=element_blank())
  
  return(bd_plot)
}

plot_latent <- function(fit_data) {
  draws <- posterior::merge_chains(fit_data$draws())
  factors <- extract_variable_array(draws, variable = "factors")[,1,,]
  f_qs <- apply(factors, c(2,3), function(x) quantile(x, c(0.025, 0.5, 0.975)))
  
  f_df <- adply(f_qs, c(2,3))
  f_df <- plyr::rename(f_df, c(X1 = "Times", X2 = "Facs", "2.5%" = "cil", "97.5%" = "ciu", "50%" = "med"))
  
  f_df$Times <- as.numeric(f_df$Times)
  f_df$Covs <- as.numeric(f_df$Facs)
  
  fplot <- ggplot(f_df) +
    geom_line(ggplot2::aes(x = Times, y = med), linewidth = 1) +
    geom_ribbon(ggplot2::aes(x = Times, ymin = cil, ymax = ciu), alpha = 0.5) +
    facet_wrap(ggplot2::vars(Covs), nrow = 1) +
    theme_bw()
  
  return(fplot)
  
}

plot_loadings <- function(fit_data) {
  draws <- posterior::merge_chains(fit_data$draws())
  factors <- extract_variable_array(draws, variable = "factor_loadings")[,1,,]
  f_qs <- apply(factors, c(2,3), function(x) quantile(x, c(0.025, 0.5, 0.975)))
  
  f_df <- adply(f_qs, c(2,3))
  f_df <- plyr::rename(f_df, c(X1 = "Times", X2 = "Facs", "2.5%" = "cil", "97.5%" = "ciu", "50%" = "med"))
  
  f_df$Times <- as.numeric(f_df$Times)
  f_df$Covs <- as.numeric(f_df$Facs)
  
  fplot <- ggplot(f_df) +
    geom_line(ggplot2::aes(x = Times, y = med), linewidth = 1) +
    geom_ribbon(ggplot2::aes(x = Times, ymin = cil, ymax = ciu), alpha = 0.5) +
    facet_wrap(ggplot2::vars(Covs), nrow = 1) +
    theme_bw()
  
  return(fplot)
  
}

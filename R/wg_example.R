get_wg_long <- function() {
  data <- haven::read_dta(
    system.file("extdata", "repgermany.dta", package = "sc.interface")
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

#' Run backends on West German reunification data
#' 
#' Runs all three backends on the classic West German reunification example.
#' Synthetic control estimates are constructed for the effect of reunification
#' on the West German GDP. By default, estimates are loaded in from pre-computed
#' data, but analysis can be optionally rerun.
#' 
#' @param refit Should the synthetic control estimates be recomputed from each
#'  available backend.
#' 
#' @return List of inference outputs from each of the three available backend
#'  methods. 
#' @export 
run_wg_example <- function(refit = FALSE) {

  if (refit) {
    wg_data_long <- get_wg_long()
    wg_fits <- run_sc(
      wg_data_long, 
      response = "gdp",
      panel_index = c("country", "year"),
      covars = c(),
      num_latent = 10,
      treated_index = "D",
      stan_iter = 4000,
      include_covars = c()
    )
  } else {
    load(system.file("data", "wg_gs_fit.RData", package = "sc.interface"))
    load(system.file("data", "wg_bp_fit.RData", package = "sc.interface"))
    load(system.file("data", "wg_stan_fit_meta.RData", package = "sc.interface"))
    stan_files = c(
      system.file("extdata", "wg_stan_c1.csv", package = "sc.interface"),
      system.file("extdata", "wg_stan_c2.csv", package = "sc.interface"),
      system.file("extdata", "wg_stan_c3.csv", package = "sc.interface"),
      system.file("extdata", "wg_stan_c4.csv", package = "sc.interface")
    )
    wg_stan_fit_post <- cmdstanr::as_cmdstan_fit(stan_files, check_diagnostics = FALSE)
    wg_stan_fit_gq <- stan_model$generate_quantities(
      wg_stan_fit_post, data = wg_stan_fit_meta$stan_data
    )
    wg_fits <- list(
      st_fit = list(
        posterior = wg_stan_fit_gq,
        meta = wg_stan_fit_meta
      ),
      bp_fit = wg_bp_fit,
      gs_fit = wg_gs_fit
    )
  }

  return(wg_fits)  
}
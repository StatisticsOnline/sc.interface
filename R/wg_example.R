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
param_names <- c(
  "random_effects_flat_lower",
  "random_effects_flat_upper",
  "factors_0",
  "factors_autocor_0",
  # "time_coefs",
  # "time_coefs_autocor",
  # "overall_sd_0",
  # "frac_var_latent",
  "causal_effects_0",
  "spillover_effects_0"
)

diagnose_stan_sc <- function(stanfit) {

  writeLines("Checking Stan synthetic control fit...\n\n")
  
  post_samples <- stanfit$draws(variables = param_names, format = "draws_matrix")
  param_names <- colnames(post_samples)

  sample_sds <- apply(post_samples, 2, sd)
  bad_sd_index <- which(sample_sds > 2 | sample_sds < 0.2)

  writeLines("Parameters on scales very different from 1 can slow sampling.\n")
  writeLines("The following parameters had standard deviations very different from 1:\n")
  for (bsi in bad_sd_index) {
    pm <- round(mean(post_samples[, bsi]), 4)
    writeLines(paste0(param_names[bsi], ": ", round(sample_sds[bsi], 3), " (", pm, ")\n"))
  }
}

plot_autocor <- function(stanfit, i = 1) {
  draws <- stanfit$draws(
    c("factors", "factors_0", "factors_autocor", "factors_autocor_0"),
    format = "draws_df"
  )
  f01 <- draws[[paste0("factors_0[1,", i, "]")]]
  f02 <- draws[[paste0("factors_0[2,", i, "]")]]
  fac <- draws[[paste0("factors_autocor_0[", i ,"]")]]
  plot(f01 ~ fac, pch = 20)

  print(cor(f01, f02))
  print(cor(f01[fac > quantile(fac, 0.8)], f02[fac > quantile(fac, 0.8)]))

}

# max_rhat <- function(stanfit) {
#   posterior::rhat()
# }

check_factor_overlap <- function(stanfit) {
  factors_flat <- stanfit$draws("factors", format = "draws_matrix")
  factors <- posterior::extract_variable_array(
    posterior::merge_chains(factors_flat),
    variable = "factors"
  )[, 1, , ]

  cors <- function(x) {
    xcor <- cor(x)
    cors <- xcor[upper.tri(xcor)]
    return(cors)
  }

  fcors <- apply(factors, 1, cors)
  print(sort(abs(rowMeans(abs(fcors)))))
  
}

check_mixing <- function(stanfit) {
  draws <- stanfit$draws()
  diags <- posterior::summarise_draws(draws, c("rhat", "ess_bulk"))
  return(c(
    max(as.numeric(diags$rhat), na.rm = TRUE),
    min(as.numeric(diags$ess_bulk), na.rm = TRUE)
  ))
}
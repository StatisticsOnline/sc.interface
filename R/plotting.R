#' @export 
plot_treated <- function(fits, level = 0.95, start = 1, compare = NULL, ylim = NULL) {

  lb_name <- paste0(50 * (1 -level), "%")
  ub_name <- paste0(100 - (50 * (1 - level)), "%")
  plot_stat <- "eff"

  if (!is.null(compare)) {
    methods <- c("stan", "bpcausal", "gsynth")
    comp_colors <- character(length = 3)
    comp_ind <- which(methods == compare)
    comp_colors[comp_ind] <- "blue"
    comp_colors[-c(comp_ind)] <- c("grey65", "wheat2")
    names(comp_colors) <- methods
  }

  plot_df <- summarize(fits, stats = plot_stat, level = level)
  plot_df$res <- plot_df[[plot_stat]]
  plot_df <- plot_df |> dplyr::filter(time_index >= start)

  tplot <- ggplot2::ggplot(data = plot_df) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data[[lb_name]], x = time_index, color = method),
      linetype = "dashed"
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = .data[[ub_name]], x = time_index, color = method),
      linetype = "dashed"
    ) +
    ggplot2::geom_line(ggplot2::aes(y = res, x = time_index, color = method)) +
    ggplot2::xlab("Time") +
    ggplot2::ylab("Effect") +
    ggplot2::theme_bw()

  if(!is.null(ylim)) {
    tplot <- tplot + ggplot2::ylim(ylim)
  }

  if(is.null(compare)) {
    tplot <- tplot + ggplot2::scale_color_brewer(palette = "Dark2")
  } else {
    tplot <- tplot + ggplot2::scale_color_manual(values = comp_colors)
    tplot <- tplot + ggplot2::geom_ribbon(
      data = dplyr::filter(plot_df, method == compare),
      ggplot2::aes(x = time_index, ymax = .data[[ub_name]], ymin = .data[[lb_name]]),
      alpha = 0.04, fill = "blue"
    )
  }

  print(tplot)
}
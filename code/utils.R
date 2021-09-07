# utililty fxns

run_SS_plots <- function(dir) {
  r4ss::SS_plots(r4ss::SS_output(dir = dir, verbose  = FALSE, printstats = F,
                                 hidewarn = TRUE), verbose = FALSE)
}

#' create relative error plots for the paramters specified
#' 
#' @param paramters A vector of parameter names, all of which should end with 
#'  _re
#' @param scalar_dat_wide The scalar dataframe created by get_results_all in the
#'  wide format.
#' @return A ggplot object.
plot_rel_error <- function(parameters, scalar_dat_wide) {
  error <-  scalar_dat_wide[, c("ID", "scenario", parameters)]
  error$scenario_fac <- factor(error$scenario, levels = unique(error$scenario))
  error_tidy <- error[,c("ID", "scenario_fac", "scenario", parameters)]
  error_tidy <- error_tidy %>% 
    gather("Parameter", "Relative_Error", (3 + seq_along(parameters)))
  
  boxplot <- ss3sim::plot_boxplot(error_tidy,
                          x = "scenario",
                          y = "Relative_Error",
                          horiz = "Parameter") +
    theme_classic(base_size = 15) +
    geom_hline(yintercept = 0, color = "blue")
  invisible(boxplot)
}

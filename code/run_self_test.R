# This is the central script to run analyses.

# Load packages, scripts, set options ----

#devtools::install_github("ss3sim/ss3sim@4d84586")
#devtools::install_github("r4ss/r4ss@2d7044e")
library(ss3sim)
library(dplyr)
library(tidyr)
library(ggplot2)

# create folders (if dont exist) ----
out <- "results_self_test"
dir.create(out)

# create the df scenarios ----
df <- setup_scenarios_defaults()
df[, "user_recdevs"] <- "matrix(0, nrow = 200, ncol = 20)"
df[, "co.par_name"] <- "SR_sigmaR"
df[, "co.par_int"] <- 0.001
df[, "ce.par_name"] <- "SR_sigmaR"
df[, "ce.par_int"] <- 0.001
df[, "ce.forecast_num"] <- 0
df[, "si.years.2"] <- "seq(60, 100, by = 2)"
df[, "si.sds_obs.2"] <- 0.001
df[, "scenarios"] <- c("self_test_1")
df[, "bias_adjust"] <- TRUE
df[, "hess_always"] <- TRUE

# Run simulations --------------------------------------------------------------
#run ss3sim in a different directory, but make sure to reset the wd on exit.
# made this into a function mostly so that changing back the wd is not forgotten!

# ' use run ss3sim when changing wd
#' @param results_wd The working directory to run the ss3sim analysis in
#' @noRd
run_analysis <- function(n_iter, simdf, results_wd) {
  wd <- getwd()
  setwd(results_wd)
  on.exit(setwd(wd), add = TRUE)
  scen_name <- run_ss3sim(iterations = n_iter, simdf = simdf)
}

scen_name <- run_analysis(n_iter = 1:20, 
                 simdf = df,
                 results_wd = out)

# Get results ------------------------------------------------------------------
res <- get_results_all(directory = out)
res$scalar <- convert_to_wide(res$scalar)
res$ts <- convert_to_wide(res$ts)

# make plots ------

plot_path <- file.path(out, "plots")
dir.create(plot_path)
scalar_dat <- calculate_re(res$scalar, add = TRUE)
ts_dat <- calculate_re(res$ts, add = TRUE)

# get error for each scenario
growth_error <- scalar_dat[, c("ID", "scenario","VonBert_K_Fem_GP_1_re", 
                               "L_at_Amin_Fem_GP_1_re", "L_at_Amax_Fem_GP_1_re")]
growth_error$scenario_fac <- factor(growth_error$scenario, levels = unique(growth_error$scenario))
growth_error_tidy <- growth_error %>% 
  select(ID, scenario_fac, scenario, VonBert_K_Fem_GP_1_re, 
         L_at_Amin_Fem_GP_1_re, L_at_Amax_Fem_GP_1_re) %>% 
  gather("Parameter", "Relative_Error", 4:6)
# may want to fiter out any scalar_dat$max_grad with a value higher than some
# pre defined quantity? Need to define convergence somehow.
scalar_dat$max_grad # for now, just print this, but don't save.

# Create plots -----------------------------------------------------------------
# Make boxplots of relative error in growth parameters, for ones that were
# estimated.

g_boxplot <- plot_boxplot(growth_error_tidy,
                                 x = "scenario",
                                 y = "Relative_Error",
                                 horiz = "Parameter")
g_boxplot + theme_classic()
ggsave(file.path(plot_path, "growth_re.png"), 
       width = 7, height = 5, units = "in")

# try some other plots
SSB_plot <- plot_lines(ts_dat, y = "SpawnBio_re", vert = "scenario")
SSB_plot + theme_classic()
ggsave(file.path(plot_path, "SSB_re.png"), width = 7, height = 5,
       units = "in")

# example of a cum mean plot. May want to look at other params, too.
cummean_plot <- ss3sim::plot_cummean(data = scalar_dat, 
                  var = "VonBert_K_Fem_GP_1_em", group = "scenario")
cummean_plot$plot + theme_classic()
ggsave(file.path(plot_path, "cummean.png"), width = 7, height = 5,
       units = "in")

#TODO: diagnose convergence problems: looks like there are some problems with 
# models not finishing. Figure out why - maybe try  turning off bias adjust to
# see if this fixes any issues? Maybe estimating too many parameters is an 
# issue?




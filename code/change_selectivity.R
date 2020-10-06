# This is the central script to run analyses.

# Load packages, scripts, set options ----

#devtools::install_github("ss3sim/ss3sim@4d84586")
#devtools::install_github("r4ss/r4ss@2d7044e")
library(ss3sim)
library(dplyr)
library(tidyr)
library(ggplot2)

# create folders (if dont exist) ----
out <- "results_self_test_2"
dir.create(out)

# create the df scenarios ----
df <- setup_scenarios_defaults()
df[, "co.par_name"] <- "SR_sigmaR"
df[, "co.par_int"] <- 0.001
df[, "ce.par_name"] <- "c('SR_sigmaR','L_at_Amin_Fem_GP_1', 'L_at_Amax_Fem_GP_1')"
df[, "ce.par_int"] <- "c(0.001, 20, 132)"
df[, "ce.par_phase"] <- "rep(-1, 3)"
df[, "ce.forecast_num"] <- 0
df[, "si.years.2"] <- "seq(60, 100, by = 2)"
df[, "si.sds_obs.2"] <- 0.001
df[, "scenarios"] <- c("self_test_1")
df[, "bias_adjust"] <- FALSE
df[, "hess_always"] <- FALSE


# Note: EM runs slow in this one (model took 22 mins to run)
df <- rbind(df, df)
df[2, "co.par_name"] <- "c('SizeSel_P4_Fishery(1)', 'SizeSel_P6_Fishery(1)')"
df[2, "co.par_int"] <- "c(15, 999)"
df[2, "ce.par_name"] <- "c('SizeSel_P4_Fishery(1)', 'SizeSel_P6_Fishery(1)')"
df[2, "ce.par_int"] <- "c(5.100, -5)"
df[, "ce.par_phase"] <- "rep(-1, 2)"
df[2, "scenarios"] <- "om_log_em_dome"

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

scen_name <- run_analysis(n_iter = 1,
                          simdf = df,
                          results_wd = out)



# This is the central script to run analyses.

# Load packages, scripts, set options ----

#devtools::install_github("ss3sim/ss3sim@6a31af3")
#devtools::install_github("r4ss/r4ss@d4430de")
library(ss3sim)
library(dplyr)
library(tidyr)
library(ggplot2)

# create folders (if dont exist) ----
out <- "results_piner"
dir.create(out)

# create the df scenarios ----
df <- setup_scenarios_defaults()
df[,"cf.fvals.1"] <- "rnorm(75, 0.2, 0.08)"
df[, "co.par_name"] <- "c('NatM_p_1_Fem_GP_1','SR_BH_steep','SR_sigmaR','SR_LN(R0)','L_at_Amin_Fem_GP_1','L_at_Amax_Fem_GP_1','VonBert_K_Fem_GP_1', 'CV_young_Fem_GP_1','CV_old_Fem_GP_1')"
# is von bert k correct?
df[, "co.par_int"] <- "c(0.3, 0.75, 0.6, 9, 3, 50, 0.3/1.65,0.1,0.1)"
df[, "ce.par_name"] <- "c('NatM_p_1_Fem_GP_1','SR_BH_steep','SR_sigmaR','SR_LN(R0)','L_at_Amin_Fem_GP_1','L_at_Amax_Fem_GP_1','VonBert_K_Fem_GP_1','CV_young_Fem_GP_1','CV_old_Fem_GP_1')"
df[, "ce.par_int"] <- "c(0.3, 0.75, 0.6, 9, 3, 50, 0.3/1.65, 0.1, 0.1)"
df[, "ce.par_phase"] <- "c(-1, -1, -1, 1, 2, 2, 2, -1, -1)"

df[, grep("sl", names(df), value = TRUE)] <- NULL
df[, grep("sa\\.[[:alpha:]]*\\.2", names(df), value = TRUE)] <- NULL
df[,"sc.years.1"] <- 50
df[,"sc.Nsamp_lengths.1"]  <- 50
df[,"sc.Nsamp_ages.1"] <- 25
df[, "scenarios"] <- c("piner_1")
df[, "bias_adjust"] <- FALSE
df[, "hess_always"] <- FALSE

# Run simulations --------------------------------------------------------------
#run ss3sim in a different directory, but make sure to reset the wd on exit.
# made this into a function mostly so that changing back the wd is not forgotten!

# ' use run ss3sim when changing wd
#' @param results_wd The working directory to run the ss3sim analysis in
#' @noRd
run_analysis <- function(iter_vec, simdf, results_wd) {
  wd <- getwd()
  setwd(results_wd)
  on.exit(setwd(wd), add = TRUE)
  scen_name <- run_ss3sim(iterations = iter_vec, simdf = simdf)
}

scen_name <- run_analysis(iter_vec = 1:10,
                          simdf = df,
                          results_wd = out)

# check out scenario to see if it has converged
output <- r4ss::SS_output(file.path(out, scen_name, "3", "em"), verbose = FALSE)
r4ss::SS_plots(output, verbose = FALSE)

# Notes
# Needed to change which params were estimated and which weren't in order to get convergence
# too strict of sampling on the index caused problems? Maybe lead to not enough
# contrast in data?
# params on bounds: this happens sometimes.
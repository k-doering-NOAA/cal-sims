# This is the central script to run analyses.

# Load packages, scripts, set options ------------------------------------------
#library(devtools)
#install_github("ss3sim/ss3sim@e898580")
#install_github("r4ss/r4ss@2d7044e")
library(ss3sim)

# create folders (if dont exist) -----------------------------------------------
dir.create("results_run_sims")

# Change pars here -------------------------------------------------------------
# Any fixed values that need to be changed should be dones so here, so that they
# are not lost futher down in the script.
# put output in this folder
out <- normalizePath("results_run_sims")

# create the df scenarios ----
# make a default data frame
df_orig <- setup_scenarios_defaults()
df <- df_orig
for(i in 1:2) {
  df <- rbind(df, df_orig)
}

# edit the data frame
# use the same f for each one. Want a 2 way trip, but these values are somewhat arbitrary currently.
df[, "cf.fvals.1"] # use this default

# use the same index for each 1
df$si.years.2 <- "seq(50, 100, by = 2)"
df$si.sds_obs.2 <- "0.4"

# length comp
df$sl.Nsamp.1 <- NULL
df$sl.years.1 <- NULL
df$sl.Nsamp.2 <- 100
df$sl.years.2 <- "seq(30, 100, by = 5)"

# age comp
df$sa.Nsamp.2 <- c("10","10", "10")
df$sa.years.2 <- c("40", "seq(40, 100, by = 5)", 
                   "seq(50, 100, by = 5)")

# cal comp : change across cases
df$sc.Nsamp_lengths.2 <- "20"
df$sc.Nsamp_ages.2 <- "10"
df$sc.years.2 <- c("30", "c(30, 35)", "c(30, 35, 40, 45)")

# scenario options
df[, "bias_adjust"] <- FALSE
df[, "hess_always"] <- FALSE
df[, "scenarios"] <- c("scen_0", "scen_1", "scen_2")

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
  run_ss3sim(iterations = n_iter, simdf = simdf)
}

run_analysis(n_iter = 1, 
             simdf = df,
             results_wd = out)

# Get results ------------------------------------------------------------------
get_results_all(directory = out)

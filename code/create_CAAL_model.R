# Investigate adding CAAL little by little to see how much is needed

# Load packages, scripts, set options ----

#devtools::install_github("ss3sim/ss3sim@5f10816")
#devtools::install_github("r4ss/r4ss@4e2502c")

library(ss3sim)
library(r4ss)
library(dplyr)
library(tidyr)
library(ggplot2)

# create folders (if dont exist) ----
out <- "create_CAAL"
dir.create(out)

df <- setup_scenarios_defaults()
set.seed(345)
# remove the years where CAAL data will be used
df[, "sl.years.1"] <- "c(26:94, 96:99)"
df[, "sa.years.1"] <- "c(26:94, 96:99)"
#add the conditionals for the fishery for years 95 and 100
df[,"sc.years.1"] <- "c(95, 100)"
df[,"sc.Nsamp_lengths.1"] <- "50"
df[, "sc.Nsamp_ages.1"] <- "30"
# specify the scenario name
df[, "scenarios"] <- c("create_CAAL_model")


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

scen_name <- run_analysis(iter_vec = 1,
                          simdf = df,
                          results_wd = out)

# check that model congverged

em_output <- r4ss::SS_output(file.path(out, "create_CAAL_model", 1, "em"))
r4ss::SS_plots(em_output)

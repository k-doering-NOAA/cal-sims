
library(ss3sim)
library(dplyr)
library(tidyr)
library(ggplot2)

# create folders (if dont exist) ----
out <- "results_na_vals_df"
dir.create(out)

# create the df scenarios ----
df <- setup_scenarios_defaults()
df <- rbind(df, df)

df[2, "sl.years.2"] <- NA # b/c can't put as null
df[2, "sl.Nsamp.2"] <- NA # b/c can't put as null
df[, "scenarios"] <- c("ctl", "no_lencomps_2")
df[, "bias_adjust"] <- FALSE
df[, "hess_always"] <- FALSE

df

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

# Running scenarios and iterations sequentially.
# Running SS in ctl/1/om
# Running SS in ctl/1/em
# Running SS in no_lencomps_2/1/om
# Error in if (any(unlist(params$years) < dat_list$styr) | any(unlist(params$years) >  : 
#                                                              missing value where TRUE/FALSE needed 


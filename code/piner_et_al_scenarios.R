# This is the central script to run analyses.

# Load packages, scripts, set options ----

devtools::install_github("ss3sim/ss3sim@bb3ebc0")
devtools::install_github("r4ss/r4ss@99c3f03")
library(ss3sim)
library(r4ss)
library(dplyr)
library(tidyr)
library(ggplot2)

# create folders (if dont exist) ----
out <- "results_piner"
dir.create(out)
mod_files_path <- "mod_files"
dir.create(mod_files_path)

# create the model to use -----
om_init_path <- system.file("extdata", "models", "cod-om", package = "ss3sim")
r4ss::copy_SS_inputs(dir.old = om_init_path, dir.new = file.path(mod_files_path, "cod-om-new-bins"))

dat <- SS_readdat(file.path(mod_files_path, "cod-om-new-bins", "codOM.dat"),
                  verbose = FALSE)
dat$minimum_size <- 1
dat$maximum_size <- 71
dat$lbin_vector <- seq(20, 71, by = 3)
SS_writedat(dat, outfile = file.path(mod_files_path, "cod-om-new-bins", "codOM.dat"), 
            overwrite = TRUE, verbose = FALSE)

# create the df scenarios ----
df <- setup_scenarios_defaults()
set.seed(123)
df[,"cf.years.1"] <- "26:100"
# for now, use same values across iterations.
df[,"cf.fvals.1"] <- "rep(0.2, 75)"
df[, "co.par_name"] <- "c('NatM_p_1_Fem_GP_1','SR_BH_steep','SR_sigmaR','SR_LN(R0)','L_at_Amin_Fem_GP_1','L_at_Amax_Fem_GP_1','VonBert_K_Fem_GP_1', 'CV_young_Fem_GP_1','CV_old_Fem_GP_1','SizeSel_P5_Fishery(1)', 'SizeSel_P6_Fishery(1)')"
# is von bert k correct?
df[, "co.par_int"] <- "c(0.3, 0.75, 0.6, 9, 3, 50, 0.3/1.65,0.1,0.1, 4.99, 4.99)"
df[, "ce.par_name"] <- "c('NatM_p_1_Fem_GP_1','SR_BH_steep','SR_sigmaR','SR_LN(R0)','L_at_Amin_Fem_GP_1','L_at_Amax_Fem_GP_1','VonBert_K_Fem_GP_1','CV_young_Fem_GP_1','CV_old_Fem_GP_1','SizeSel_P5_Fishery(1)', 'SizeSel_P6_Fishery(1)', 'SizeSel_P1_Fishery(1)', 'SizeSel_P2_Fishery(1)', 'SizeSel_P3_Fishery(1)', 'SizeSel_P4_Fishery(1)')"
df[, "ce.par_int"] <- "c(0.3, 0.75, 0.6, 9, 3, 50, 0.3/1.65, 0.1, 0.1, 4.99, 4.99, 50.8, -3, 5.1, 15)"
df[, "ce.par_phase"] <- "c(-1, -1, -1, 1, 2, 2, 2, -1, -1, -1, -1, -1, -1, -1, -1)"

df[, grep("sl", names(df), value = TRUE)] <- NULL
df[, grep("sa\\.[[:alpha:]]*\\.2", names(df), value = TRUE)] <- NULL
df[,"sc.years.1"] <- 100
df[,"sc.Nsamp_lengths.1"]  <- 250
df[,"sc.Nsamp_ages.1"] <- 250
df[, "scenarios"] <- c("piner_250")
df[, "bias_adjust"] <- FALSE
df[, "hess_always"] <- FALSE
om_path <- normalizePath(file.path('mod_files', 'cod-om-new-bins'), winslash = "/")
df[,"om"] <- om_path
df[,"em"] <- system.file("extdata", "models", "cod-em", package = "ss3sim")

df <- rbind(df,df)
# add the second scenario
df[2, "scenarios"] <- "piner_4000"
df[2,"sc.Nsamp_lengths.1"] <- 4000 
df[2, "sc.Nsamp_ages.1"] <- 4000

# piner et al sample size is 250, 500, 1000, 2000, and 4000 individ-uals). Note these
# were aonly taken for 1 year o f data.

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

# Notes
# Needed to change which params were estimated and which weren't in order to get convergences
# too strict of sampling on the index caused problems? Maybe lead to not enough
# contrast in data?
# params on bounds: this happens sometimes.

# to do next: look at 10 iterations and what is going on.

# Get results ------------------------------------------------------------------
res <- get_results_all(directory = out)

res <- vector(mode = "list", length = 2)
res[[1]] <- read.csv( file.path(out, "ss3sim_scalar.csv"))
res[[2]] <- read.csv(file.path(out, "ss3sim_ts.csv"))
names(res) <-  c("scalar", "ts")
              


# Check EM convergence -----
# make sure aren't giant
res$scalar[res$scalar$model_run == "em","max_grad"] > 1
# look at params on bounds
# on iteration has a big gradient, but SSB doesn't look that off.

#some params are stuck low, which may be an issue.
# looks like it is true that these params aren't on bounds, but are fairly low -
# I don't think this should be an issue, tho?
res$scalar[res$scalar$model_run == "em",
   c("iteration", "scenario", "params_on_bound", "params_stuck_low",
     "params_stuck_high")]

# Look at SSB - how off is it? - tracking pretty well, so these sims aren't so bad.
ggplot(dat = res$ts, aes(x = year, y = SpawnBio)) +
  geom_line(aes(color = model_run)) +
  facet_grid(rows = vars(scenario), cols = vars(iteration))

# make an ss output plot for reference to get a feel fror how these sims ran....

# make plots ------
plot_path <- file.path(out, "plots")
dir.create(plot_path)
scalar_dat_wide <- convert_to_wide(res$scalar)
ts_dat_wide <- convert_to_wide(res$ts)
scalar_dat_wide <- calculate_re(scalar_dat_wide, add = TRUE)
ts_dat_wide <- calculate_re(ts_dat_wide, add = TRUE)

# get error for each scenario ----
growth_error <- scalar_dat_wide[, c("ID", "scenario","VonBert_K_Fem_GP_1_re", 
                               "L_at_Amin_Fem_GP_1_re", "L_at_Amax_Fem_GP_1_re")]

growth_error$scenario_fac <- factor(growth_error$scenario, levels = unique(growth_error$scenario))
growth_error_tidy <- growth_error %>% 
  select(ID, scenario_fac, scenario, VonBert_K_Fem_GP_1_re, 
         L_at_Amin_Fem_GP_1_re, L_at_Amax_Fem_GP_1_re) %>% 
  gather("Parameter", "Relative_Error", 4:6)

g_boxplot <- plot_boxplot(growth_error_tidy,
                          x = "scenario",
                          y = "Relative_Error",
                          horiz = "Parameter")
g_boxplot + theme_classic(base_size = 15)


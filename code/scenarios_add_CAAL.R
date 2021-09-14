# Investigate adding CAAL little by little to see how much is needed

# Load packages, scripts, set options ----

#devtools::install_github("ss3sim/ss3sim@eeb0b33")
#devtools::install_github("r4ss/r4ss@e588b87")
library(ss3sim)
library(r4ss)
library(dplyr)
library(tidyr)
library(ggplot2)

source("code/utils.R")

# create folders (if dont exist) ----
out <- "results_add_CAAL"
dir.create(out)
mod_files_path <- "mod_files"
dir.create(mod_files_path)

# create the om model to use -----
# om_init_path <- system.file("extdata", "models", "cod-om", package = "ss3sim")
# r4ss::copy_SS_inputs(dir.old = om_init_path, dir.new = file.path(mod_files_path, "cod-om-new-bins"))
# 
# dat <- SS_readdat(file.path(mod_files_path, "cod-om-new-bins", "codOM.dat"),
#                   verbose = FALSE)
# dat$minimum_size <- 1
# dat$maximum_size <- 71
# dat$lbin_vector <- seq(11, 71, by = 1)
# dat$N_lbins <- length(dat$lbin_vector)
# dat$lencomp <- dat$lencomp[,1:6]
# new_lencomp <- as.data.frame(matrix(data = 1, nrow = nrow(dat$lencomp), ncol = dat$N_lbins))
# dat$lencomp <- cbind(dat$lencomp, new_lencomp)
# colnames(dat$lencomp)[7:ncol(dat$lencomp)] <- paste0("l", dat$lbin_vector)
# SS_writedat(dat, outfile = file.path(mod_files_path, "cod-om-new-bins", "codOM.dat"), 
#             overwrite = TRUE, verbose = FALSE)
# # modify selectivity pattern for om
# om_ctl <- SS_readctl(file.path(mod_files_path, "cod-om-new-bins", "codOM.ctl"), 
#                      verbose = FALSE, use_datlist = TRUE, datlist = dat)
# om_ctl[["size_selex_types"]]["Fishery", "Pattern"] <- 0 # change to constant
# par_rows_to_rm <- grep("Fishery", rownames(om_ctl[["size_selex_parms"]]))
# om_ctl[["size_selex_parms"]] <- om_ctl[["size_selex_parms"]][-par_rows_to_rm, ]
# SS_writectl(om_ctl, file.path(mod_files_path, "cod-om-new-bins", "codOM.ctl"), 
#             verbose = FALSE, overwrite = TRUE)

# create the em model to use ----

# em_init_path <- system.file("extdata", "models", "cod-em", package = "ss3sim")
# r4ss::copy_SS_inputs(dir.old = em_init_path, dir.new = file.path(mod_files_path, "cod-em-constant-sel"))
# em_ctl <- SS_readctl(file.path(mod_files_path, "cod-em-constant-sel", "codEM.ctl"), 
#            verbose = FALSE, use_datlist = TRUE, datlist = dat)
# # edit next lines.
# em_ctl[["size_selex_types"]]["Fishery", "Pattern"] <- 0 # change to constant
# em_par_rows_to_rm <- grep("Fishery", rownames(em_ctl[["size_selex_parms"]]))
# em_ctl[["size_selex_parms"]] <- em_ctl[["size_selex_parms"]][-em_par_rows_to_rm, ]
# SS_writectl(em_ctl, file.path(mod_files_path, "cod-em-constant-sel", "codEM.ctl"), 
#             verbose = FALSE, overwrite = TRUE)

# create the df scenarios ----

# TODO: add in a good time series of lengths that can estimate selectivity well
# if the lengths are informing growth, make the data more realistic by decreasing sample size
# then, add in the CAAL scenarios and see what this does to growth

# At first, just set up 1 scenario 

df <- setup_scenarios_defaults()
set.seed(123)
df[,"cf.years.1"] <- "26:100"
# for now, use same values across iterations.
df[,"cf.fvals.1"] <- "rep(0.1, 75)"
df[, "co.par_name"] <- "c('NatM_p_1_Fem_GP_1','SR_BH_steep','SR_sigmaR','SR_LN(R0)','L_at_Amin_Fem_GP_1','L_at_Amax_Fem_GP_1','VonBert_K_Fem_GP_1', 'CV_young_Fem_GP_1','CV_old_Fem_GP_1')"
# is von bert k correct?
df[, "co.par_int"] <- "c(0.3, 0.75, 0.6, 9, 3, 50, 0.3/1.65,0.1,0.1)"
df[, "ce.par_name"] <- "c('NatM_p_1_Fem_GP_1','SR_BH_steep','SR_sigmaR','SR_LN(R0)','L_at_Amin_Fem_GP_1','L_at_Amax_Fem_GP_1','VonBert_K_Fem_GP_1','CV_young_Fem_GP_1','CV_old_Fem_GP_1', 'LnQ_base_Survey(2)')"
df[, "ce.par_int"] <- "c(0.3, 0.75, 0.6, 9, 3, 50, 0.3/1.65, 0.1, 0.1, 0)"
df[, "ce.par_phase"] <- "c(-1, -1, -1, 1, 2, 2, 2, -1, -1, 1)"
# get rid ofall length and age samples
df[, grep("sl", names(df), value = TRUE)] <- NULL
df[, grep("sa", names(df), value = TRUE)] <- NULL


df[,"cb.bin_vector"] <- "seq(11, 71, by = 1)"
df[, grep("sa\\.[[:alpha:]]*\\.2", names(df), value = TRUE)] <- NULL
df[, "hess_always"] <- FALSE
om_path <- normalizePath(file.path('mod_files', 'cod-om-new-bins'), winslash = "/")
df[,"om_dir"] <- om_path
df[,"em_dir"] <- normalizePath(file.path('mod_files', 'cod-em-constant-sel'), winslash = "/")
# leave the population bin width as in the model for now:
df[, "cb.pop_binwidth"] <- 1
df[, "cb.pop_minimum_size"] <- 1
df[, "cb.pop_maximum_size"] <- 71
df[, "cb.lbin_method"] <- 2


# setup the 4 scenarios
df <- rbind(df,df)
df <- rbind(df, df) # 4 scenarios
# add in conditional age at length samples for the survey
df[,"sc.Nsamp_lengths.2"]  <- 250
df[,"sc.Nsamp_ages.2"] <- 250
# first 
df[1,"sc.years.2"] <- 100
# for the second scenario
df[2,"sc.years.2"] <- "c(98, 100)"
# for the third
df[3,"sc.years.2"] <- "c(96, 98, 100)"
# for the fourth
df[4,"sc.years.2"] <- "c(94, 96, 98, 100)"
df[, "scenarios"] <- c("cal_1_yr", "cal_2_yr", "cal_3_yr", "cal_4_yr")


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

scen_name <- run_analysis(iter_vec = 1:100,
                          simdf = df,
                          results_wd = out)

# What happens when CAAL is input as conditionals? ----
# I think this can be added by:
# 1. copy over the EM files from a scenario with multiple years of CAAL
# 2. Transform the CAAL data into the input for marginals
# 3. Run the EM and use data weighting
# This can be written outside of ss3sim, but maybe incorporating it in would be 
# helpful
# may want to also look at having only 1 yr of CAAL retained and the rest 
# marginals, converting all to marginals, etc.
# with multiple fleets, it would be useful to look at what happens when only 
# 1 fleet is transfered to marginals; can also look at many fleet models to 
# understand this.

scen <- "cal_2_yr"
new_scen <- paste0(scen, "marg")
new_loc <- file.path(out, new_scen)
dir.create(new_loc)
scen_folders <- file.path(
  list.dirs(file.path(out, scen), recursive = FALSE, full.names = TRUE))
# copy over om and em files to a new scenario files
lapply(scen_folders, function(f, new_loc) {
  om_files <- list.files(file.path(f, "om"),full.names = TRUE, recursive = FALSE)
  em_files <- list.files(file.path(f, "om"),full.names = TRUE, recursive = FALSE)
  dir.create(file.path(new_loc, basename(f)))
  dir.create(file.path(new_loc, basename(f), "om"))
  dir.create(file.path(new_loc, basename(f), "em"))
  file.copy(from = om_files, to = file.path(new_loc, basename(f), "om", basename(om_files)))
  file.copy(from = em_files, to = file.path(new_loc, basename(f), "em", basename(om_files)))
}, new_loc = new_loc)
# tranform CAAL data into input for marginals
scen_folders <- list.dirs(new_loc, recursive = FALSE, full.names = TRUE)
ss_path <- ss3sim::get_bin()

lapply(scen_folders, function(s) {
  tmp_dat <- r4ss::SS_readdat(file.path(s, "em"))
  # change the CAAL data;
  # add up all columns for the same flt/ yr and change lbin hi and lbin lo to -1s
  # write back out file
  # run em
})
# tmp_dat <- r4ss::SS_readdat(file.path(new_loc, "1", "em")) # replace the 1

# summarize this iteration

# What happens when there is time varying growth? ----
# Add in time varying growth to OM, scenarios with and without it in the EM
# Look at what happens when these data ar put in as conditionals.


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
res$scalar[res$scalar$model_run == "em",c("scenario", "max_grad")]

# check relative SSB
calc_SSB <- res$ts %>%
  select(iteration, scenario, year, model_run, SpawnBio)
OM_vals <- calc_SSB %>% 
  filter(model_run == "om") %>% 
  rename(SpawnBio_OM = SpawnBio)
EM_vals <- calc_SSB %>% 
  filter(model_run == "em") %>% 
  rename(SpawnBio_EM = SpawnBio)
bind_vals <- full_join(EM_vals, OM_vals, by = c("iteration", "scenario", "year")) %>% 
  mutate(SSB_ratio = SpawnBio_EM/SpawnBio_OM)
filter_SSB <- bind_vals %>% 
  filter(SSB_ratio > 3 | SSB_ratio < 0.33) # standards more lax
not_converged <- unique(filter_SSB[,c("iteration", "scenario")])
nrow(not_converged) #23 by these standards, kind of high.


# remove iterations that aren't converged

res_converged <- lapply(res, function(x, not_converged) {
  converged <- dplyr::anti_join(x, not_converged, by = c("iteration", "scenario"))
  converged
}, not_converged = not_converged)

# look at params on bounds
# on iteration has a big gradient, but SSB doesn't look that off.

#some params are stuck low, which may be an issue.
# looks like it is true that these params aren't on bounds, but are fairly low -
# I don't think this should be an issue, tho?
res$scalar[res$scalar$model_run == "em" & !is.na(res$scalar$params_on_bound),
   c("iteration", "scenario", "params_on_bound", "params_stuck_low",
     "params_stuck_high")]

# Look at SSB - how off is it? - tracking pretty well, so these sims aren't so bad.
# ggplot(dat = res$ts, aes(x = year, y = SpawnBio)) +
#   geom_line(aes(color = model_run)) +
#   facet_grid(rows = vars(scenario), cols = vars(iteration))
# # some look way off, so probably want to filter out.


# make an ss output plot for reference to get a feel fror how these sims ran....

# make plots ------
plot_path <- file.path(out, "plots")
dir.create(plot_path)
scalar_dat_wide <- convert_to_wide(res_converged$scalar)
ts_dat_wide <- convert_to_wide(res_converged$ts)
scalar_dat_wide <- calculate_re(scalar_dat_wide, add = TRUE)
ts_dat_wide <- calculate_re(ts_dat_wide, add = TRUE)

# plotting function for relative error

# get error for each scenario ----
g_boxplot <- plot_rel_error(
  parameters = c("VonBert_K_Fem_GP_1_re", "L_at_Amin_Fem_GP_1_re",
                 "L_at_Amax_Fem_GP_1_re"), 
  scalar_dat_wide = scalar_dat_wide)
g_boxplot
ggsave(file.path(plot_path, "rel_error_growth.png"), plot = g_boxplot,
       width = 12, height = 8, units = "in")

# get error in key pop quantities -----
p_boxplot <- plot_rel_error(
  parameters = c("SSB_MSY_re","TotYield_MSY_re", "SSB_Unfished_re"), 
  scalar_dat_wide = scalar_dat_wide)
p_boxplot
ggsave(file.path(plot_path, "rel_error_pop.png"), plot = p_boxplot,
       width = 12, height = 8, units = "in")

# make selectivity plots, eventually. ---

s_boxplot <- plot_rel_error(parameters = c("Size_DblN_ascend_se_Survey_2_re", 
                              "Size_DblN_peak_Survey_2_re"), scalar_dat_wide = scalar_dat_wide)
s_boxplot
ggsave(file.path(plot_path, "rel_error_sel.png"), plot = s_boxplot,
       width = 12, height = 8, units = "in")



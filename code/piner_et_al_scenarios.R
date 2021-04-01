# This is the central script to run analyses.

# Load packages, scripts, set options ----

devtools::install_github("ss3sim/ss3sim@issue_286")
devtools::install_github("r4ss/r4ss@2e5f0d5")
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

# create the om model to use -----
om_init_path <- system.file("extdata", "models", "cod-om", package = "ss3sim")
r4ss::copy_SS_inputs(dir.old = om_init_path, dir.new = file.path(mod_files_path, "cod-om-new-bins"))

dat <- SS_readdat(file.path(mod_files_path, "cod-om-new-bins", "codOM.dat"),
                  verbose = FALSE)
dat$minimum_size <- 1
dat$maximum_size <- 71
dat$lbin_vector <- seq(11, 71, by = 1)
dat$N_lbins <- length(dat$lbin_vector)
dat$lencomp <- dat$lencomp[,1:6]
new_lencomp <- as.data.frame(matrix(data = 1, nrow = nrow(dat$lencomp), ncol = dat$N_lbins))
dat$lencomp <- cbind(dat$lencomp, new_lencomp)
colnames(dat$lencomp)[7:ncol(dat$lencomp)] <- paste0("l", dat$lbin_vector)
SS_writedat(dat, outfile = file.path(mod_files_path, "cod-om-new-bins", "codOM.dat"), 
            overwrite = TRUE, verbose = FALSE)
# modify selectivity pattern for om
om_ctl <- SS_readctl(file.path(mod_files_path, "cod-om-new-bins", "codOM.ctl"), 
                     verbose = FALSE, use_datlist = TRUE, datlist = dat)
om_ctl[["size_selex_types"]]["Fishery", "Pattern"] <- 0 # change to constant
par_rows_to_rm <- grep("Fishery", rownames(om_ctl[["size_selex_parms"]]))
om_ctl[["size_selex_parms"]] <- om_ctl[["size_selex_parms"]][-par_rows_to_rm, ]
SS_writectl(om_ctl, file.path(mod_files_path, "cod-om-new-bins", "codOM.ctl"), 
            verbose = FALSE, overwrite = TRUE)

# create the em model to use ----

em_init_path <- system.file("extdata", "models", "cod-em", package = "ss3sim")
r4ss::copy_SS_inputs(dir.old = em_init_path, dir.new = file.path(mod_files_path, "cod-em-constant-sel"))
em_ctl <- SS_readctl(file.path(mod_files_path, "cod-em-constant-sel", "codEM.ctl"), 
           verbose = FALSE, use_datlist = TRUE, datlist = dat)
# edit next lines.
em_ctl[["size_selex_types"]]["Fishery", "Pattern"] <- 0 # change to constant
em_par_rows_to_rm <- grep("Fishery", rownames(em_ctl[["size_selex_parms"]]))
em_ctl[["size_selex_parms"]] <- em_ctl[["size_selex_parms"]][-em_par_rows_to_rm, ]
SS_writectl(em_ctl, file.path(mod_files_path, "cod-em-constant-sel", "codEM.ctl"), 
            verbose = FALSE, overwrite = TRUE)

# create the df scenarios ----
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
# get rid of fishery length samples
df[, "sl.years.1"] <- NA
df[, "sl.Nsamp.1"] <- NA
df[, "sa.Nsamp.2"] <- NA
df[, "sa.years.2"] <- NA
df[, "sa.Nsamp.2"] <- 50
df[, "sa.years.2"] <- "26:100"

# keep the length samples
#df[, grep("sl", names(df), value = TRUE)] <- NULL

df[,"cb.bin_vector"] <- "seq(11, 71, by = 1)"
df[, grep("sa\\.[[:alpha:]]*\\.2", names(df), value = TRUE)] <- NULL
df[,"sc.years.1"] <- 100
df[,"sc.Nsamp_lengths.1"]  <- 250
df[,"sc.Nsamp_ages.1"] <- 250
df[, "scenarios"] <- c("bins_1_250_cal")
df[, "bias_adjust"] <- FALSE
df[, "hess_always"] <- FALSE
om_path <- normalizePath(file.path('mod_files', 'cod-om-new-bins'), winslash = "/")
df[,"om"] <- om_path
df[,"em"] <- normalizePath(file.path('mod_files', 'cod-em-constant-sel'), winslash = "/")
# leave the population bin width as in the model for now:
df[, "cb.pop_binwidth"] <- 1
df[, "cb.pop_minimum_size"] <- 1
df[, "cb.pop_maximum_size"] <- 71
df[, "cb.lbin_method"] <- 2


df <- rbind(df,df)
df <- rbind(df, df, df, df) # 8 scenarios
# add the second scenario
df[2, "scenarios"] <- "bins_1_4000_cal"
#df[2, "cb.bin_vector"] <- "seq(11, 71, by = 1)" # don't need this line, already set
df[2,"sc.Nsamp_lengths.1"] <- 4000
df[2, "sc.Nsamp_ages.1"] <- 4000
# 3rd scen
df[3, "scenarios"] <- "bins_3_250_cal"
df[3,"cb.bin_vector"] <- "seq(11, 71, by = 3)"
# 4th scen
df[4, "scenarios"] <- "bins_3_4000_cal"
df[4,"cb.bin_vector"] <- "seq(11, 71, by = 3)"
df[4,"sc.Nsamp_lengths.1"] <- 4000
df[4, "sc.Nsamp_ages.1"] <- 4000
# 5th scen
df[5, "scenarios"] <- "bins_6_250_cal"
df[5, "cb.bin_vector"] <- "seq(11, 71, by = 6)"
# 6th scen
df[6, "scenarios"] <- "bins_6_4000_cal"
df[6, "cb.bin_vector"] <- "seq(11, 71, by = 6)"
df[6,"sc.Nsamp_lengths.1"] <- 4000
df[6, "sc.Nsamp_ages.1"] <- 4000
# 7 scen: marginals only
df[7, "scenarios"] <- "bins_1_250_marg"
df[7:8, grep("sc\\.", names(df), value = TRUE)] <- NA
df[7, "sl.Nsamp.1"] <- 250
df[7, "sl.years.1"] <- 100
df[7, "sa.Nsamp.2"] <- 250
df[7, "sa.years.2"] <- 100
# 8 scen
df[8, "scenarios"] <- "bins_1_4000_marg"
df[8, "sl.Nsamp.1"] <- 4000
df[8, "sl.years.1"] <- 100
df[8, "sa.Nsamp.2"] <- 4000
df[8, "sa.years.2"] <- 100

# piner et al sample size is 250, 500, 1000, 2000, and 4000 individ-uals). Note these
# were aonly taken for 1 year o f data.


# Run simulations --------------------------------------------------------------
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
res$scalar[res$scalar$model_run == "em",c("scenario", "max_grad")]
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

# make selectivity plots, eventually.


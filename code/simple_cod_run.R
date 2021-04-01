# Simple cod run that sets up unbiased "deterministic" runs
# turn off process error (from recdevs) and turn down sampling error.

#devtools::install_github("ss3sim/ss3sim@66923d7")

library(ss3sim)
library(ggplot2)
library(dplyr)
library(tidyr)

df <- setup_scenarios_defaults()
df[,"scenarios"] <- "base_run_some_survey_err" # I hope
# for now, use same values across iterations.
df[,"cf.years.1"] <- "26:100"
df[,"cf.fvals.1"] <- "c(seq(0.01, 0.12, length.out = 25), seq(0.12, 0.005, length.out = 50) )"
df[, "user_recdevs"] <- "matrix(0, nrow = 200, ncol = 20)" # turn off recdevs
# Make the sigma R small and fix it in the model
df[, "co.par_name"] <- "SR_sigmaR"
df[, "co.par_int"] <- 0.001
df[, "ce.par_name"] <- "SR_sigmaR"
df[, "ce.par_int"] <- 0.001
df[, "ce.par_phase"] <- -1 
# Turn off bias adjustment and hessian
df[, "bias_adjust"] <- FALSE #not necessary b/c no rec devs used
df[, "hess_always"] <- FALSE 
# setup 1 more scenario that reduces the survey error to be really low. 
# note this will result in higher max gradients for some of the em runs.
#(in the single digits instead of e-05)
df <- rbind(df,df)
df[2,"scenarios"] <- "base_run_low_survey_err"
# add the following 2 line to make the index obs more precise
df[2, "si.sds_obs.2"] <- 0.001
out <- file.path("results_simple_cod")
dir.create(out)

scen_name <- run_analysis(iter_vec = 1:10,
                          simdf = df,
                          results_wd = out)

res <- get_results_all(directory = out)

# run the following if get_results_all has already been created:
# res <- vector(mode = "list", length = 2)
# res[[1]] <- read.csv( file.path(out, "ss3sim_scalar.csv"))
# res[[2]] <- read.csv(file.path(out, "ss3sim_ts.csv"))
# names(res) <-  c("scalar", "ts")

# Check EM convergence -----
# make sure aren't giant.
res$scalar[res$scalar$model_run == "em",c("scenario", "max_grad")]
# look at params on bounds

res$scalar[res$scalar$model_run == "em",
           c("iteration", "scenario", "params_on_bound", "params_stuck_low",
             "params_stuck_high")]

# Look at SSB trajectories in om and em
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

# look at relative error in key quantities -------
MSY_err <- scalar_dat_wide[, c("ID", "scenario","SSB_MSY_re","depletion_re", 
                                 "Totbio_Unfished_re", "VonBert_K_Fem_GP_1_re",
                                "L_at_Amin_Fem_GP_1_re", "L_at_Amax_Fem_GP_1_re")]

# get error for each scenario ----
growth_error <- scalar_dat_wide[, c("ID", "scenario", "SSB_MSY_re", "depletion_re", 
                                    "Totbio_Unfished_re", "VonBert_K_Fem_GP_1_re", 
                                    "L_at_Amin_Fem_GP_1_re", "L_at_Amax_Fem_GP_1_re")]

growth_error$scenario_fac <- factor(growth_error$scenario, levels = unique(growth_error$scenario))
growth_error_tidy <- growth_error %>% 
  select(ID, scenario_fac, scenario, SSB_MSY_re, depletion_re, 
         Totbio_Unfished_re,  VonBert_K_Fem_GP_1_re, 
         L_at_Amin_Fem_GP_1_re, L_at_Amax_Fem_GP_1_re) %>% 
  gather("Parameter", "Relative_Error", 4:9)

g_boxplot <- plot_boxplot(growth_error_tidy,
                          x = "scenario",
                          y = "Relative_Error",
                          horiz = "Parameter")
g_boxplot + 
  theme_classic(base_size = 15) + 
  geom_hline(yintercept = 0, color = "red")

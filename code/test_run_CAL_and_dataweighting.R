# This is the central script to run analyses.


# change this to look at data weighting (this can replace run_sims instead, if
# that is better)
# Load packages, scripts, set options ------------------------------------------
library(devtools)
#install_github("ss3sim/ss3sim@CAL") # use this branch for CAL
library(ss3sim)
library(tidyverse)
library(r4ss)
#source("code/funs/tune_comp.R")

# create folders (if dont exist) -----------------------------------------------
dir.create("casefiles")
dir.create("code")
dir.create("docs")
dir.create("results")

# Change pars here -------------------------------------------------------------
# Any fixed values that need to be changed should be dones so here, so that they
# are not lost futher down in the script.
# put output in this folder
out <- normalizePath("results")

# Set up the casefiles ---------------------------------------------------------
# F. Want a 2 way trip, but these values are somewhat arbitrary currently.
F_case <- c("# description: 2 way trip", 
            "years; 1:100", 
            "fisheries; 1",
            "fvals; c(rep(0,25), seq(0.005, 0.15, length.out =  50), seq(0.15, 0.06, length.out = 25))")
writeLines(F_case, con = file.path("casefiles", "F0-cod.txt"))
# index. Assume a value for every other year for now? Need to check if this is 
# realistic. Only include survey obs
index_case <- c("fleets; c(2)", 
                "years; list(seq(20,100, by = 2))", 
                "sds_obs; list(0.4)") # is this a good default value?
writeLines(index_case, con = file.path("casefiles", "index3-cod.txt"))

# age comp.
age_case <- c("fleets; c(2)", 
              "Nsamp; list(20)",
              "years; list(seq(30, 100, by = 5))") # is this a good default value?
writeLines(age_case, con = file.path("casefiles", "agecomp3-cod.txt"))

# calcomp
calcomp_case <- c("# Note that Nsamp_ages must be less than Nsamp_lengths",
                  "fleets; c(2)",
                  "years; list(seq(30, 45, by = 5))",
                  "Nsamp_lengths; list(30)",
                  "Nsamp_ages; list(20)",
                  "ESS_lengths; list(600)",
                  "ESS_ages; list(500)",
                  "method; 'simple_random' #other option: length_stratified")

writeLines(calcomp_case, file.path("casefiles", "calcomp3-cod.txt"))

# length comp. Assume collected every 5 years for now. Required in CAL years.
# Note need a separate file for each scenario although they are all the same
len_case <- c("fleets; c(1,2)",
              "Nsamp; list(100,100)",  # I think ss3sim treats this as number of fish? and not N tows?
              "years; list(seq(30, 100, by = 5), seq(30, 100, by = 5)) ")
writeLines(len_case, file.path("casefiles", "lcomp3-cod.txt"))


# Set up for simulations -------------------------------------------------------
# use absolute paths so it doesn't matter what the wd is.
scen <- expand_scenarios(cases = list(D = 3, F = 0), 
                         species = "cod")
extdat <- system.file("extdata", package = "ss3sim")
om <- file.path(extdat, "models", "cod-om")
em <- file.path(extdat, "models", "cod-em")
case_folder <- normalizePath("casefiles", mustWork = TRUE) # get the absolute path.
# Run simulations --------------------------------------------------------------
#run ss3sim in a different directory, but make sure to reset the wd on exit.
# made this into a function mostly so that changing back the wd is not forgotten!

# ' use run ss3sim when changing wd
#' @param results_wd The working directory to run the ss3sim analysis in
#' @noRd
run_analysis <- function(n_iter, scen, case_folder, om, em, case_files, results_wd) {
  wd <- getwd()
  setwd(results_wd)
  on.exit(setwd(wd), add = TRUE)
  run_ss3sim(iterations = n_iter,
             scenarios = scen,
             case_folder = case_folder,
             om_dir = om,
             em_dir = em,
             case_files = case_files
  )
}

args <- get_caseargs(folder = case_folder, scenario = scen,
                     case_files = list(F = "F",
                                       D = c("index", "lcomp", "calcomp", "agecomp")),
                    )
# add list component the 
DM_pars <- list(method = "DM", 
                                 niters_weighting = 1,
                                 fleets = 2)
MI_pars <- list(method = "MI", 
                niters_weighting = 1,
                fleets = 2)
Fran_pars <- list(method = "Francis",
                  niters_weighting = 1,
                  fleets = 2)

wd <- getwd()
setwd(out)
ss3sim_base(iterations = 1:10,
            scenarios = "D3-F0-cod-DM",
            f_params = args$F,
            index_params = args$index,
            lcomp_params = args$lcomp,
            agecomp_params = args$agecomp,
            calcomp_params = args$calcomp,
            wtatage_params = NULL,
            mlacomp_params = NULL,
            em_binning_params = NULL,
            estim_params = NULL,
            tv_params = NULL,
            operat_params = NULL,
            om_dir = om,
            em_dir = em,
            retro_params = NULL,
            data_params = NULL,
            weight_comps_params = DM_pars,
            user_recdevs = NULL,
            user_recdevs_warn = TRUE,
            bias_adjust = FALSE,
            hess_always = FALSE,
            print_logfile = TRUE)

ss3sim_base(iterations = 1:10,
            scenarios = "D3-F0-cod-MI",
            f_params = args$F,
            index_params = args$index,
            lcomp_params = args$lcomp,
            agecomp_params = args$agecomp,
            calcomp_params = args$calcomp,
            wtatage_params = NULL,
            mlacomp_params = NULL,
            em_binning_params = NULL,
            estim_params = NULL,
            tv_params = NULL,
            operat_params = NULL,
            om_dir = om,
            em_dir = em,
            retro_params = NULL,
            data_params = NULL,
            weight_comps_params = MI_pars,
            user_recdevs = NULL,
            user_recdevs_warn = TRUE,
            bias_adjust = FALSE,
            hess_always = FALSE,
            print_logfile = TRUE)

ss3sim_base(iterations = 1:10,
            scenarios = "D3-F0-cod-Fran",
            f_params = args$F,
            index_params = args$index,
            lcomp_params = args$lcomp,
            agecomp_params = args$agecomp,
            calcomp_params = args$calcomp,
            wtatage_params = NULL,
            mlacomp_params = NULL,
            em_binning_params = NULL,
            estim_params = NULL,
            tv_params = NULL,
            operat_params = NULL,
            om_dir = om,
            em_dir = em,
            retro_params = NULL,
            data_params = NULL,
            weight_comps_params = Fran_pars,
            user_recdevs = NULL,
            user_recdevs_warn = TRUE,
            bias_adjust = FALSE,
            hess_always = FALSE,
            print_logfile = TRUE)

# No data weighting
ss3sim_base(iterations = 1:10,
            scenarios = "D3-F0-cod-No_DW",
            f_params = args$F,
            index_params = args$index,
            lcomp_params = args$lcomp,
            agecomp_params = args$agecomp,
            calcomp_params = args$calcomp,
            wtatage_params = NULL,
            mlacomp_params = NULL,
            em_binning_params = NULL,
            estim_params = NULL,
            tv_params = NULL,
            operat_params = NULL,
            om_dir = om,
            em_dir = em,
            retro_params = NULL,
            data_params = NULL,
            weight_comps_params = NULL,
            user_recdevs = NULL,
            user_recdevs_warn = TRUE,
            bias_adjust = FALSE,
            hess_always = FALSE,
            print_logfile = TRUE)

setwd(wd)
# Note DM gives an error msg about data file having 2 sections...could not replicate
# from building pkg locally and running. not sure why this is the case.
# Get results ------------------------------------------------------------------

get_results_all(directory = out, 
                user_scenarios = paste0("D3-F0-cod-", c("DM", "MI", "Fran", 
                                                        "No_DW")))
# read in csvs
scal_res <- read.csv(file.path(out, "ss3sim_scalar.csv"), stringsAsFactors = FALSE)
ts_res <- read.csv(file.path(out, "ss3sim_ts.csv"), stringsAsFactors = FALSE)

# convergence diagnostics


plot_path <- file.path(out, "plots", "CAL_and_data_weight")
dir.create(plot_path)

ggplot(scal_res, aes(x = max_grad)) + # Not really that helpful for determining convergence.
  geom_histogram(fill = "#56B4E9", color = "black", bins = 10) +
  facet_wrap(~scenario)+
  theme_classic()
ggsave(file.path(plot_path, "max_gradient_by_scen.png"), width = 12, height = 8, units = "in")

convergence <- select(scal_res, ID, params_on_bound_em, params_stuck_high_em, params_stuck_low_em)
write.csv(convergence, file.path(out, "D3-F0-cod-convergence.csv"))

# Size_DblN_peak_Fishery(1) still on bounds.

# calculate relative error for each growth parameter
growth <- scal_res %>%
           select(matches("L_at_Am|VonBert_K|^ID$|^scenario$|^iteration$")) %>% 
           mutate(L_at_Amin_rel_err = 
                    (L_at_Amin_Fem_GP_1_om - L_at_Amin_Fem_GP_1_em)/
                    L_at_Amin_Fem_GP_1_om) %>% 
           mutate(L_at_Amax_rel_err = 
                   (L_at_Amax_Fem_GP_1_om - L_at_Amax_Fem_GP_1_em)/
                   L_at_Amax_Fem_GP_1_om) %>% 
           mutate(VonBert_K_rel_err = 
                   (VonBert_K_Fem_GP_1_om - VonBert_K_Fem_GP_1_em)/
                   VonBert_K_Fem_GP_1_om) %>% 
           select(matches("^ID$|^scenario$|^iteration$|rel_err")) %>% 
           gather("var", "rel_err", 4:6)
  


#Plot the relative error in the growth parameters
plot_scalar_boxplot(growth,  "scenario", "rel_err", horiz = "var", axes.free = F, ) +
  geom_hline(yintercept = 0, color = "blue")+
  theme_classic()
ggsave(file.path(plot_path, "growth_re_boxplot.png"), width = 6, height = 4, units = "in", )

ggplot(growth, aes(rel_err)) +
  geom_histogram(fill = "#56B4E9", color = "black", bins = 10) +
  geom_vline(xintercept = 0, color = "red")+
  facet_grid(rows = vars(var), cols = vars(scenario))+
  theme_classic()
ggsave(file.path(plot_path, "growth_re_hist.png"), width = 12, height = 8, units = "in")


# stability ----
#uses ggplot, dplyr packages
scal_res <- read.csv(file.path(out, "ss3sim_scalar.csv"), stringsAsFactors = FALSE)
# have enough iterations been done?
vonbert_select <- scal_res %>% 
  select(iteration, scenario, VonBert_K_Fem_GP_1_em) %>% 
  arrange(scenario, iteration) %>% 
  group_by(scenario) %>% 
  mutate(cummean_K = cummean(VonBert_K_Fem_GP_1_em)) %>% 
  ungroup()

ggplot(vonbert_select, aes(x = iteration, y = cummean_K))+
  geom_line(aes(color = scenario))+
  geom_point(aes(color = scenario)) +
  #ylim(0, 0.25)+ # use if want to set manual limits 
  theme_classic()
ggsave(file.path(plot_path, "cummean_K.png"), width = 12, height = 8, units = "in")

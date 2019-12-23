# plot the results of code/run_sims.R
# please restart the r session before using this script.

# Load packages, scripts, set options ------------------------------------------
#library(devtools)
#install_github("ss3sim/ss3sim@CAL")
library(ss3sim)
library(dplyr)
library(tidyr)
library(ggplot2)
options(stringsAsFactors = FALSE)
# specify paths, create folders ------------------------------------------------
plot_path <- file.path("results", "plots")
# path to the summarized output.
sum_path <- file.path("results", "summarized_output")
dir.create(plot_path)

# read in output ---------------------------------------------------------------
scalar_dat <- read.csv(file.path("results", "ss3sim_scalar.csv"))
ts_dat <- read.csv(file.path("results", "ss3sim_ts.csv"))

# Manipulate output ------------------------------------------------------------
# add relative error
scalar_dat <- calculate_re(scalar_dat, add = TRUE)
ts_dat <- calculate_re(ts_dat, add = TRUE)

# get error for each scenario
growth_error <- scalar_dat[, c("ID", "scenario","VonBert_K_Fem_GP_1_re", 
                               "L_at_Amin_Fem_GP_1_re", "L_at_Amax_Fem_GP_1_re")]
growth_error$scenario_fac <- factor(growth_error$scenario, levels = unique(growth_error$scenario))
growth_error_tidy <- growth_error %>% 
                      select(ID, scenario_fac, scenario, VonBert_K_Fem_GP_1_re, 
                             L_at_Amin_Fem_GP_1_re, L_at_Amax_Fem_GP_1_re) %>% 
                      gather("Parameter", "Relative_Error", 4:6)
# may want to fiter out any scalar_dat$max_grad with a value higher than some
# pre defined quantity? Need to define convergence somehow.
scalar_dat$max_grad # for now, just print this, but don't save.

# Create plots -----------------------------------------------------------------
# Make boxplots of relative error in growth parameters, for ones that were
# estimated.

g_boxplot <- plot_scalar_boxplot(growth_error_tidy,
                                 x = "scenario",
                                 y = "Relative_Error",
                                 horiz = "Parameter")
g_boxplot + theme_classic()
ggsave(file.path("results", "plots", "growth_re.png"), 
       width = 7, height = 5, units = "in")

# try some other plots
SSB_plot <- plot_ts_lines(ts_dat, y = "SpawnBio_re", vert = "scenario")
SSB_plot + theme_classic()
ggsave(file.path("results", "plots", "SSB_re.png"), width = 7, height = 5,
       units = "in")

# Is there an easy way to toggle par estimation in the EM? A method in ss3sim 
# to do this, perhaps?

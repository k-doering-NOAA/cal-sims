# This is the central script to run analyses.

# Load packages, scripts, set options ------------------------------------------
#library(devtools)
#install_github("ss3sim/ss3sim@CAL") # use this branch for CAL
library(ss3sim)

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
                "years; list(1:100)", 
                "sds_obs; list(0.4)") # is this a good default value?
writeLines(index_case, con = file.path("casefiles", "index0-cod.txt"))
writeLines(index_case, con = file.path("casefiles", "index1-cod.txt"))
writeLines(index_case, con = file.path("casefiles", "index2-cod.txt"))
# age comp. Don't want for now, so make a null case.
# will need to make different cases depending on how much CAL data is included.

# scenario 0: Assume all marginals every 5 yrs.
agecomp_case <- c("# description: Null case", 
                  "Nsamp; list(10)", 
                  "fleets; c(2)", 
                  "years; list(seq(30, 100, by = 5))", 
                  "cpar; NULL")
writeLines(agecomp_case, con = file.path("casefiles", "agecomp0-cod.txt"))
# scenario 1: 2 yrs CAL data, rest of years marginal age comps.
agecomp_case[4] <- "years; list(seq(40,100, by = 5))" #maybe would be better not at beginning of time series?
writeLines(agecomp_case, con = file.path("casefiles", "agecomp1-cod.txt"))
# Scenario 2: 4 years CAL data, rest of years marginal age comps
agecomp_case[4] <- "years; list(seq(50,100, by = 5))" #maybe would be better not at beginning of time series?
writeLines(agecomp_case, con = file.path("casefiles", "agecomp2-cod.txt"))
# calcomp
# Scenario 0: none needed (or can create a null scenario)
calcomp_case <- c("fleets; NULL",
                  "Nsamp; NULL", 
                  "years; NULL")
writeLines(calcomp_case, file.path("casefiles", "calcomp0-cod.txt"))
# Scenario 1:
calcomp_case <- c("fleets; c(2)", 
                  "Nsamp; list(20)", # what is a realistic value for this?
                  "years; list(c(30, 35))")
writeLines(calcomp_case, file.path("casefiles", "calcomp1-cod.txt"))
# Scenario 2:
calcomp_case[3] <- "years; list(c(30, 35, 40, 45))"
writeLines(calcomp_case, file.path("casefiles", "calcomp2-cod.txt"))

# length comp. Assume collected every 5 years for now. Required in CAL years.
# Note need a separate file for each scenario although they are all the same
len_case <- c("fleets; c(2)",
              "Nsamp; list(100)",  # I think ss3sim treats this as number of fish? and not N tows?
              "years; list(seq(30, 100, by = 5)) ")
writeLines(len_case, file.path("casefiles", "lcomp0-cod.txt"))
writeLines(len_case, file.path("casefiles", "lcomp1-cod.txt"))
writeLines(len_case, file.path("casefiles", "lcomp2-cod.txt"))

# Set up for simulations -------------------------------------------------------
# use absolute paths so it doesn't matter what the wd is.
scen <- expand_scenarios(cases = list(D = 0:2, F = 0), 
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

run_analysis(n_iter = 1:5, scen = scen, case_folder = case_folder, om = om,
             em = em, 
             case_files = list(F = "F",
                               D = c("index", "lcomp", "agecomp", "calcomp")), 
             results_wd = out)

# Get results ------------------------------------------------------------------
get_results_all(directory = out)

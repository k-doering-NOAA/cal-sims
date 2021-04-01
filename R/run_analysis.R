#' run ss3sim in a different directory, but make sure to reset the wd on exit.
#'
#' made this into a function mostly so that changing back the wd is not forgotten!
#' use run ss3sim when changing wd
#' @param iter_vec The desired iterations
#' @param simdf Data frame of specs to run the simulation
#' @param results_wd The working directory to run the ss3sim analysis in
#' @noRd
run_analysis <- function(iter_vec, simdf, results_wd) {
  wd <- getwd()
  setwd(results_wd)
  on.exit(setwd(wd), add = TRUE)
  scen_name <- run_ss3sim(iterations = iter_vec, simdf = simdf)
}

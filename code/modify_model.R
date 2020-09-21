# 

modify_model <- function(input_mod = "cod", life_history = "flatfish", 
                         out_path = getwd()) {
  life_history <- match.arg(arg = life_history,
                            choices = "flatfish",
                            several.ok = FALSE)
  if(input_mod == "cod") {
    input_mod <- system.file("extdata", "models", "cod-om", package = "ss3sim")
  }
  # read in files -----
  start <- r4ss::SS_readstarter(file.path(input_mod, "starter.ss"),
                                verbose = FALSE)
  dat <- r4ss::SS_readdat(file.path(input_mod, start$datfile), verbose = FALSE)
  ctl <- r4ss::SS_readctl(file.path(input_mod, start$ctlfile), verbose = FALSE, 
                          use_datlist = TRUE, datlist = dat)
  #need to modify forecast file?
  # edit
  if(life_history == "flatfish") {
    # edit ctl file ----
    # edit MG par bounds, initial value, est or not
    ctl$MG_parms["L_at_Amin_Fem_GP_1","LO"] <- 2
    ctl$MG_parms["L_at_Amin_Fem_GP_1","HI"] <- 20
    ctl$MG_parms["L_at_Amin_Fem_GP_1","INIT"] <- 12.6666
    ctl$MG_parms["L_at_Amax_Fem_GP_1","LO"] <- 25
    ctl$MG_parms["L_at_Amax_Fem_GP_1","HI"] <- 80
    ctl$MG_parms["L_at_Amax_Fem_GP_1","INIT"] <- 47.4245
    ctl$MG_parms["VonBert_K_Fem_GP_1","LO"] <- 0.01
    ctl$MG_parms["VonBert_K_Fem_GP_1","HI"] <- 2
    ctl$MG_parms["VonBert_K_Fem_GP_1","INIT"] <- 0.347769
    ctl$MG_parms["CV_young_Fem_GP_1", "LO"] <- 0.01
    #ctl$MG_parms["CV_young_Fem_GP_1", "HI"]
    ctl$MG_parms["CV_young_Fem_GP_1", "INIT"] <- 0.2
    ctl$MG_parms["CV_old_Fem_GP_1", "LO"] <- 0.01
    #ctl$MG_parms["CV_old_Fem_GP_1", "HI"]
    ctl$MG_parms["CV_young_Fem_GP_1", "INIT"] <- 0.2
    ctl$MG_parms["Wtlen_1_Fem_GP_1", "INIT"] <- 0.00001
    
    #TODO: edit other lines
    
    
    # edit SR par bounds, initial value, est or not
    # need to edit q or selectivity?
    # edit dat file ----
    # edit pop bins
    dat$binwidth <- 1
    dat$minimum_size <- 2
    dat$maximum_size <- 86
    # edit length bins
    dat$lbin_vector <- seq(from = 4, to = 86, by = 2)
    dat$N_lbins <- length(dat$lbin_vector)
  }
  
  # consider diffs needed between the OM and EM
  # write out files -----
  start$datfile <- "data.ss"
  start$ctlfile <- "control.ss"
  mod_path <- file.path(out_path, paste0(life_history, "-om"))
  dir.create(mod_path)
  r4ss::SS_writectl(ctl, file.path(mod_path, "control.ss"), verbose = FALSE)
  r4ss::SS_writedat(dat, file.path(mod_path, "data.ss"), verbose = FALSE)
  r4ss::SS_writestarter(start, dir = mod_path, verbose = FALSE)
  
  return(TRUE)
}

modify_model()
library(reticulate)
library(plotly)

source("server_files/load_input_helpers.R")
source("server_files/helpers.R")
source("server_files/plot_functions.R")

user        <- Sys.info()['user']
reticulate::use_python(paste0("/Users/",user,"/myenv/bin/python"), required = TRUE)

reticulate::source_python("helpers.py")
reticulate::source_python("moltenprot_shiny.py")

dsf <- DsfFitter()
dsf$load_supr_dsf("/Users/oburastero/Downloads/example_data.supr")

dsf$set_signal(dsf$signals[1])

dsf$estimate_fluo_derivates(5)

dsf$decompose_spectra()


dsf$set_signal(dsf$signals[1])

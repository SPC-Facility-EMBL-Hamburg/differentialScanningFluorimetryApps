library(reticulate)
library(plotly)
library(tidyverse)

source("server_files/load_input_helpers.R")
source("server_files/helpers.R")
source("server_files/plot_functions.R")

user        <- Sys.info()['user']
#reticulate::use_python(paste0("/Users/",user,"/myenv/bin/python"), required = TRUE)
reticulate::use_python("/home/os/myenv/bin/python", required = TRUE)

reticulate::source_python("helpers.py")
reticulate::source_python("main.py")

dsf <- ManyDsfFitters()

dsf$add_experiment('/home/os/thermochemicalDenaturationApp/appFiles/Chemelt/www/panta.xlsx')
dsf$add_experiment('/home/os/thermochemicalDenaturationApp/appFiles/Chemelt/www/MX3005P.txt')

print(dsf$all_signals)

dsf$unify_signals(
    signal_list = c("350nm","Fluorescence"),
    new_signal_name = "Unified_Signal"
)

print(dsf$all_signals)

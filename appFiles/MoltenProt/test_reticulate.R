library(reticulate)
library(plotly)
library(tidyverse)

source("server_files/load_input_helpers.R")
source("server_files/helpers.R")
source("server_files/plot_functions.R")

user        <- Sys.info()['user']
reticulate::use_python(paste0("/Users/",user,"/myenv/bin/python"), required = TRUE)

reticulate::source_python("helpers.py")
reticulate::source_python("main.py")

dsf <- ManyDsfFitters()
dsf$add_experiment('./www/demo.xlsx')

dsf$set_signal(dsf$all_signals[1])

dsf$set_colors(rep('blue',48))

dsf$select_conditions(c(rep(TRUE,4),rep(FALSE,44)))

dsf$estimate_fluo_derivates()
dsf$set_baseline_types(2,2)
dsf$estimate_baselines_parameters()
dsf$equilibrium_two_state(0)

params_all <- dsf$get_experiment_properties('params_all',flatten=TRUE)
params_name <- dsf$get_experiment_properties('params_name',flatten=TRUE)
params_name <- unique(params_name)

fitted_conditions <- dsf$get_experiment_properties('fitted_conditions',flatten=TRUE)

get_sorted_params_table(
params_all,
fitted_conditions,
16,
params_name,
'Tm'
)


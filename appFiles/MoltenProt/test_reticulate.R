library(reticulate)
library(plotly)

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

dsf$select_conditions(c(rep(TRUE,4),rep(FALSE,44)))

dsf$estimate_fluo_derivates()
dsf$set_baseline_types(2,2)
dsf$estimate_baselines_parameters()
dsf$equilibrium_two_state(0)

dfs <- list()
for (exp in dsf$available_experiments) {
  py_obj <- dsf$experiments[[exp]]
  fluo <- py_obj$fluo
  conds <- py_obj$conditions
  temps <- py_obj$temps
  print(conds)
  df <- make_df4plot_row_wise(fluo,conds,temps)
  colnames(df)[colnames(df) == 'temp'] <- 'temperature'
  df$temperature <- df$temperature - 273.15 # Kelvin to Centigrade
  dfs[[length(dfs)+1]] <- df
}

df <- colbind_pad(dfs)

print(head(df))
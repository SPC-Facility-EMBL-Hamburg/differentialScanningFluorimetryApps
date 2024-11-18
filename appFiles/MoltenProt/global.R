packages <- c('reshape2','Cairo','tidyverse','reticulate','pracma','data.table',
              'grid',"plotly","shinyalert","shinydashboard","shinycssloaders",
              "rhandsontable","tableHTML")

invisible(lapply(packages, library, character.only = TRUE))

appName     <- "MoltenProt"
user        <- Sys.info()['user']

reticulate::use_python(paste0("/home/",user,"/myenv/bin/python"), required = TRUE)

# developer path
base_dir <- paste0("/home/",user,"/spc_shiny_servers/differentialScanningFluorimetryApps/appFiles/",appName,"/")

# path for the docker user
if (user == 'shiny') {
  base_dir <- paste0("/home/shiny/",appName,'/')
}

global_chunck_n     <- 16 # should be global_plot_columns * global_plot_rows
global_plot_columns <- 4
global_plot_rows    <- 4

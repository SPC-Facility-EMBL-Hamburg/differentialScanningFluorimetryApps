packages <- c('reshape2','Cairo','tidyverse','reticulate','pracma','data.table',
              'scales','viridis','minpack.lm','broom','grid',"plotly","shinyalert",
              "shinydashboard","shinycssloaders","rhandsontable","tableHTML")

invisible(lapply(packages, library, character.only = TRUE))

appName     <- "FoldAffinity"
user        <- Sys.info()['user']

reticulate::use_python(paste0("/home/",user,"/myenv/bin/python"), required = TRUE)

# developer path
base_dir <- paste0("/home/",user,"/differentialScanningFluorimetryApps/appFiles/",appName,"/")

# path for the docker user
if (user == 'shiny') {
  base_dir <- paste0("/home/shiny/",appName,'/')
}

# Number of elements per plot in the fluorescence fit
# (the plots are plot_columns X plot_rows grids

global_n_rows_conditions_table <- 24
global_chunck_n <- 16    # should be global_plot_columns * global_plot_rows

global_plot_columns <- 4
global_plot_rows    <- 4

# viridis color palette
global_colors_palette_signal <- c(
  "#440154FF", "#46075BFF", "#470E61FF", "#481467FF", "#481A6CFF", "#482071FF",
  "#482576FF", "#472A7AFF", "#462F7EFF", "#453581FF", "#443A83FF", "#424086FF",
  "#404588FF", "#3E4989FF", "#3C4F8AFF", "#3A538BFF", "#38588CFF", "#365C8DFF",
  "#34618DFF", "#32658EFF", "#30698EFF", "#2E6E8EFF", "#2C718EFF", "#2B758EFF",
  "#297A8EFF", "#277E8EFF", "#26828EFF", "#24868EFF", "#238A8DFF", "#218E8DFF",
  "#20928CFF", "#1F968BFF", "#1F9A8AFF", "#1F9F88FF", "#1FA287FF", "#21A685FF",
  "#25AB82FF", "#28AE80FF", "#2DB27DFF", "#33B679FF", "#3ABA76FF", "#40BD72FF",
  "#49C16DFF", "#52C569FF", "#5AC864FF", "#64CB5FFF", "#6ECE58FF", "#77D153FF",
  "#82D34CFF", "#8DD645FF", "#98D83EFF", "#A3DA37FF", "#AFDD2FFF", "#BADE28FF",
  "#C6E021FF", "#D1E21BFF", "#DDE318FF", "#E8E419FF", "#F3E61EFF", "#FDE725FF"
)

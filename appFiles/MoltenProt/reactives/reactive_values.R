reactives <- reactiveValues(
  data_loaded                    = FALSE, # To control the display of plots/tables
  multiple_files                 = FALSE, # To control the display of the unify signals button
  have_plot_data                 = FALSE, # To control the display of plots/tables
  nconditions                    = 0,     # Number of total conditions imported
  report_was_created             = FALSE, # To activate the download report button
  global_max_conditions          = 384,   # To determine the load input table number of columns.
  global_n_rows_conditions_table = 96,    # Number of rows for the load input table
  full_spectra                   = FALSE, # Full spectra instead of single wavelength
  spectra_panel_names            = NULL,
  reportDir                      = NULL,
  model_is_two_state             = NULL,
  model_name                     = NULL,
  min_wl                         = 0,
  max_wl                         = 800,
  data_was_fitted                = FALSE,
  fluo_fit_real                  = NULL,
  fluo_fit_pred                  = NULL
  )

output$data_loaded             <- reactive({
  return(reactives$data_loaded)})

output$report_was_created             <- reactive({
  return(reactives$report_was_created)})

output$full_spectra             <- reactive({
  return(reactives$full_spectra)})

output$data_was_fitted  <- reactive({
  return(reactives$data_was_fitted)})

output$model_is_two_state  <- reactive( { return(reactives$model_is_two_state) } )

output$model_name  <- reactive( { return(reactives$model_name) } )

output$multiple_files <- reactive({
  return(reactives$multiple_files)})

outputOptions(output, "data_loaded"       , suspendWhenHidden = FALSE)
outputOptions(output, "report_was_created", suspendWhenHidden = FALSE)
outputOptions(output, "full_spectra"      , suspendWhenHidden = FALSE)
outputOptions(output, "model_is_two_state", suspendWhenHidden = FALSE)
outputOptions(output, "model_name", suspendWhenHidden = FALSE)
outputOptions(output, "data_was_fitted", suspendWhenHidden = FALSE)
outputOptions(output, "multiple_files"    , suspendWhenHidden = FALSE)

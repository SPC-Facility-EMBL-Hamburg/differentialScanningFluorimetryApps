
output$download_params_table <- downloadHandler(filename = function() {
  "Fitted_Parameters.csv"},content = function(file) {

  params_all <- dsf$get_experiment_properties('params_all',flatten=TRUE)
  params_name <- dsf$get_experiment_properties('params_name',flatten=TRUE)
  params_name <- unique(params_name)

  fitted_conditions <- dsf$get_experiment_properties('fitted_conditions',flatten=TRUE)

  df <- get_sorted_params_table(
    params_all,
    fitted_conditions,
    global_chunck_n,
    params_name,
    input$sort_table_parameter
  )

    write.csv(df,file,row.names = F,quote = F)
  })

output$download_params_errors_table <- downloadHandler(filename = function() { 
    "Fitted_Parameters_Relative_Errors.csv"},content = function(file) {

    errors_percentage_all <- dsf$get_experiment_properties('errors_percentage_all',flatten=TRUE)
  fitted_conditions <- dsf$get_experiment_properties('fitted_conditions',flatten=TRUE)
  params_name <- dsf$get_experiment_properties('params_name',flatten=TRUE)
  params_name <- unique(params_name)

  df <- get_sorted_params_table_errors(
    errors_percentage_all,
    fitted_conditions,
    params_name,
    input$sort_table_parameter
  )

  write.csv(df, file,row.names = F,quote = F)
  })

output$download_fitted_cond_table <- downloadHandler(filename = function() { 
  "Fitted_Conditions.csv"},content = function(file) {

  conditions <- dsf$get_experiment_properties('conditions',flatten=TRUE)
    fitted_conditions <- dsf$get_experiment_properties('fitted_conditions',flatten=TRUE)

  write.csv(get_fitted_conditions_table(
      conditions,fitted_conditions),
      file,row.names = F,quote = F)
  })

output$download_signal_plot <- downloadHandler(filename = function() {
  "signal_plot_data.csv"},content = function(file) {
    
    if (input$data_export_format == "row_wise") {

      dfs <- list()
      for (exp in dsf$available_experiments) {
        py_obj <- dsf$experiments[[exp]]
        fluo <- py_obj$fluo
        conds <- py_obj$conditions
        temps <- py_obj$temps
        df <- make_df4plot_row_wise(fluo,conds,temps)
        colnames(df)[colnames(df) == 'temp'] <- 'temperature'
        df$temperature <- df$temperature - 273.15 # Kelvin to Centigrade
        dfs[[length(dfs)+1]] <- df
      }

      df <- colbind_pad(dfs)
      
    } else {
      df <- py_dsf_to_df(dsf)
      colnames(df) <- c("temperature","condition","signal","original_name")
      df$temperature <- df$temperature - 273.15 # Kelvin to Centigrade

    }

    write.csv(df, 
      file,row.names = F,quote = F)
  })


output$download_derivative_plot <- downloadHandler(filename = function() {
  "1st_derivative_plot_data.csv"},content = function(file) {
    
    if (input$data_export_format == "row_wise") {

      dfs <- list()
      for (exp in dsf$available_experiments) {
        py_obj <- dsf$experiments[[exp]]
        fluo <- py_obj$derivative
        conds <- py_obj$conditions
        temps <- py_obj$temps
        df <- make_df4plot_row_wise(fluo,conds,temps)
        colnames(df)[colnames(df) == 'temp'] <- 'temperature'
        df$temperature <- df$temperature - 273.15 # Kelvin to Centigrade
        dfs[[length(dfs)+1]] <- df
      }

      df <- colbind_pad(dfs)
      
    } else {
      df <- py_dsf_to_df(dsf,mode = 'derivative')
      colnames(df) <- c("temperature","condition","derivative","original_name")
      df$temperature <- df$temperature - 273.15 # Kelvin to Centigrade

    }

    write.csv(df, file,row.names = F,quote = F)
    
  })

output$download_2derivative_plot <- downloadHandler(filename = function() { 
  "2nd_derivative_plot_data.csv"},content = function(file) {
    
    if (input$data_export_format == "row_wise") {

      dfs <- list()
      for (exp in dsf$available_experiments) {
        py_obj <- dsf$experiments[[exp]]
        fluo <- py_obj$derivative2
        conds <- py_obj$conditions
        temps <- py_obj$temps
        df <- make_df4plot_row_wise(fluo,conds,temps)
        colnames(df)[colnames(df) == 'temp'] <- 'temperature'
        df$temperature <- df$temperature - 273.15 # Kelvin to Centigrade
        dfs[[length(dfs)+1]] <- df
      }

      df <- colbind_pad(dfs)
      
    } else {
      df <- py_dsf_to_df(dsf,mode = 'derivative2')
      colnames(df) <- c("temperature","condition","derivative","original_name")
      df$temperature <- df$temperature - 273.15 # Kelvin to Centigrade

    }

    write.csv(df, file,row.names = F,quote = F)
  })

output$download_max_derivative_plot <- downloadHandler(filename = function() { 
  "max_derivative_plot_data.csv"},content = function(file) {

    tms_from_deriv <- dsf$get_experiment_properties('tms_from_deriv',flatten=TRUE)
    conditions     <- dsf$get_experiment_properties('conditions',flatten=TRUE)

    df <- generate_max_der_df(tms_from_deriv,conditions)
    colnames(df) <- c("condition","transition_temperature")
    write.csv(df, file,row.names = F,quote = F)
  })

# Download the plots
output$download_fit_plots = downloadHandler(
  filename = 'fitting_plots.zip',
  content = function(file){
    
    # Set temporary working directory
    owd <- setwd( tempdir())
    on.exit( setwd( owd))

    fitted_conditions <- dsf$get_experiment_properties('fitted_conditions',flatten=TRUE)

    withBusyIndicatorServer("download_fit_plots",{
      total_plots <- ceiling(length(fitted_conditions)/global_chunck_n)
      
      iter <- 0
      fns <- c()
      for (selected in 1:total_plots) {
        iter <- iter + 1
        real_data  <- fluo_fit_data()$fluo_fit_real
        model_data <- fluo_fit_data()$fluo_fit_pred
        fname <- paste0("fitting_plot_",iter,".png")
        fns   <- c(fns,fname)
        plot_fluorescence_fit(real_data,model_data,selected,TRUE) %>% 
          ggpubr::ggexport(filename = fname,res=240,width = 1800,height = 1800)
      }
    })
    
    # Zip them up
    zip( file, fns)
  }
)

output$download_unfolded_fraction_plot <- downloadHandler(filename = function() { 
  "fraction_unfolded.csv"},content = function(file) {
    
    selected_indexes <- filter_conditions()

    params_name <- dsf$get_experiment_properties('params_name',flatten=TRUE)
    params_all  <- dsf$get_experiment_properties('params_all',flatten=TRUE)
    params_name <- unique(params_name)
    fitted_conditions <- dsf$get_experiment_properties('fitted_conditions',flatten=TRUE)

    tm_position <- which(params_name == "Tm")
    tms <- sapply(params_all, function(x) x[tm_position])
    
    dh_position <- which(params_name == "DH")
    dhs <- sapply(params_all, function(x) x[dh_position])
    
    dhs        <- dhs[selected_indexes]
    tms        <- tms[selected_indexes]
    conditions <- fitted_conditions[selected_indexes]
    
    dfs <- list()
    
    i <- 0
    for (temp in seq(25,95,1)) {
      
      i <- i+1
      fus <- mapply(get_fraction_unfolded_EquilTwoState, dhs, tms, temp+273.15)
      df_unfolded <- data.frame("fu"=fus,"Condition"=conditions)
      df_unfolded$temp <- temp
      
      df_unfolded$Condition <- sapply(df_unfolded$Condition,trimws,"both")
      df_unfolded$cond_     <- c(paste0(1:nrow(df_unfolded),"_",df_unfolded$Condition))
      
      dfs[[i]] <- df_unfolded
      
    }
    
    df_tog <- do.call(rbind,dfs)
    
    # Avoid duplicate names
    df_tog <- avoid_positions_with_the_same_name_in_df(df_tog)
    
    df_tog <- df_tog[,c(1:4)]
    
    colnames(df_tog) <- c("Condition","Condition_internal_ID","Unfolded_fraction","Temperature")
    
    df_tog <- df_tog %>% select(-Condition_internal_ID)
    
    write.csv(df_tog, file,row.names = F,quote = F)
  })



# The downloadHandler function has an intrinsic timeout limit. We need to generate first the report
observeEvent(input$downloadReport,{
  req(reactives$data_was_fitted)
  reactives$report_was_created             <- FALSE

  withBusyIndicatorServer("downloadReport",{ 
    
    reactives$reportDir <- paste0(tempfile(),'/')
    dir.create(reactives$reportDir)

    filename <- paste0(reactives$reportDir,input$filename_report,".pdf")

    file.copy(paste0(base_dir,"report_template/report.Rmd"), reactives$reportDir, overwrite = TRUE)
    file.copy(paste0(base_dir,"report_template/header.tex"), reactives$reportDir, overwrite = TRUE)
    
    rmarkdown::render(input = paste0(reactives$reportDir,'report.Rmd'), 
                      output_format = 'pdf_document',output_file=filename)
    
    reactives$report_was_created             <- TRUE
  })
})

output$downloadReportHidden <- downloadHandler(
  filename = function() {
    paste0(input$filename_report,".pdf")
  },
  
  content <- function(file) {
    file.copy(paste0(reactives$reportDir,input$filename_report,".pdf"), file)
  }
)




output$download_params_table <- downloadHandler(filename = function() {
  "Fitted_Parameters.csv"},content = function(file) {
    write.csv(get_sorted_params_table(
      dsf$params_all,dsf$fitted_conditions,global_chunck_n,
      dsf$params_name,input$sort_table_parameter), 
      file,row.names = F,quote = F)
  })

output$download_params_errors_table <- downloadHandler(filename = function() { 
    "Fitted_Parameters_Relative_Errors.csv"},content = function(file) {
      write.csv(get_sorted_params_table_errors(
        dsf$errors_percentage_all,dsf$fitted_conditions,
        dsf$params_name,input$sort_table_parameter), 
        file,row.names = F,quote = F)
  })

output$download_fitted_cond_table <- downloadHandler(filename = function() { 
  "Fitted_Conditions.csv"},content = function(file) {
    write.csv(get_fitted_conditions_table(
      dsf$conditions,dsf$fitted_conditions), 
      file,row.names = F,quote = F)
  })

output$download_signal_plot <- downloadHandler(filename = function() {
  "signal_plot_data.csv"},content = function(file) {
    
    if (input$data_export_format == "row_wise") {
      df <- make_df4plot_row_wise(dsf$fluo,dsf$conditions,dsf$temps)
      colnames(df)[colnames(df) == 'temp'] <- 'temperature'
      
    } else {
      df <- make_df4plot(dsf$fluo,dsf$conditions,dsf$temps)
      colnames(df) <- c("temperature","condition","signal","original_name")
    }
    
    df$temperature <- df$temperature - 273.15 # Kelvin to Centigrade
    
    write.csv(df, 
      file,row.names = F,quote = F)
  })


output$download_derivative_plot <- downloadHandler(filename = function() {
  "1st_derivative_plot_data.csv"},content = function(file) {
    
    if (input$data_export_format == "row_wise") {
      df <- make_df4plot_row_wise(dsf$derivative,dsf$conditions,dsf$temps)
      colnames(df)[colnames(df) == 'temp'] <- 'temperature'
      
    } else {
      df <- make_df4plot(dsf$derivative,dsf$conditions,dsf$temps)
      colnames(df) <- c("temperature","condition","derivative","original_name")
    }
    
    df$temperature <- df$temperature - 273.15 # Kelvin to Centigrade
    
    write.csv(df, 
              file,row.names = F,quote = F)
    
  })

output$download_2derivative_plot <- downloadHandler(filename = function() { 
  "2nd_derivative_plot_data.csv"},content = function(file) {
    
    if (input$data_export_format == "row_wise") {
      df <- make_df4plot_row_wise(dsf$derivative2,dsf$conditions,dsf$temps)
      colnames(df)[colnames(df) == 'temp'] <- 'temperature'
      
    } else {
      df <- make_df4plot(dsf$derivative2,dsf$conditions,dsf$temps)
      colnames(df) <- c("temperature","condition","derivative","original_name")
    }
    
    df$temperature <- df$temperature - 273.15 # Kelvin to Centigrade
    
    write.csv(df, 
              file,row.names = F,quote = F)
  })

output$download_max_derivative_plot <- downloadHandler(filename = function() { 
  "max_derivative_plot_data.csv"},content = function(file) {
    
    df <- generate_max_der_df(dsf$tms_from_deriv,dsf$conditions)
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
    
    withBusyIndicatorServer("download_fit_plots",{
      total_plots <- ceiling(length(dsf$fitted_conditions)/global_chunck_n)
      
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
    
    tm_position <- which(dsf$params_name == "Tm")  
    tms <- sapply(dsf$params_all, function(x) x[tm_position])
    
    dh_position <- which(dsf$params_name == "dHm")  
    dhs <- sapply(dsf$params_all, function(x) x[dh_position])
    
    dhs        <- dhs[selected_indexes]
    tms        <- tms[selected_indexes]
    conditions <- dsf$fitted_conditions[selected_indexes]
    
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
  req(fluo_fit_data())
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



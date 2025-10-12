
# Reactive values to track available data
reactives <- reactiveValues(
  tm_fit_done = FALSE
)

# Allow the UI to know the values of these reactives
output$tm_fit_done <- reactive({
  return(reactives$tm_fit_done)
})

outputOptions(output, "tm_fit_done" , suspendWhenHidden = FALSE)

# End of reactive values to track available data

# Render signal plot
output$signal <- renderPlotly({
  
  req(input$table1)
  if (is.null(modify_fluo_temp_cond())) {return(NULL)}
  
  fluo_m <- make_df4plot(dsf$fluo,dsf$concentrations,dsf$temps)
  if (!(is.null(fluo_m))) {
    return( plot_fluo_signal(

      fluo_m=fluo_m,
      which=dsf$signal_type,
      plot_width=input$plot_width,
      plot_height=input$plot_height,
      plot_type=input$plot_type,
      legend_text_size=input$plot_legend_size,
      axis_size=input$plot_axis_size,
      show_x_grid=input$show_x_grid,
      show_y_grid=input$show_y_grid,
      show_axis_lines=input$show_axis_lines,
      tickwidth=input$plot_tickwidth,
      ticklen=input$plot_ticklen,
      line_width=input$plot_line_width,
      color_bar_length=input$color_bar_length,
      color_bar_orientation=input$color_bar_orientation,
      xLegend=input$x_legend_pos,
      yLegend=input$y_legend_pos

      )
    )
    }
  
  return(NULL)
})

# Render signal plot
output$first_der <- renderPlotly({
  
  req(input$table1)
  if (is.null(modify_fluo_temp_cond())) {return(NULL)}
  
  fluo_m <- make_df4plot(dsf$derivative,dsf$concentrations,dsf$temps)

  if (!(is.null(fluo_m))) {
    return( plot_fluo_signal(
      fluo_m=fluo_m,
      which="First derivative",
      plot_width=input$plot_width,
      plot_height=input$plot_height,
      plot_type=input$plot_type,
      legend_text_size=input$plot_legend_size,
      axis_size=input$plot_axis_size,
      show_x_grid=input$show_x_grid,
      show_y_grid=input$show_y_grid,
      show_axis_lines=input$show_axis_lines,
      tickwidth=input$plot_tickwidth,
      ticklen=input$plot_ticklen,
      line_width=input$plot_line_width,
      color_bar_length=input$color_bar_length,
      color_bar_orientation=input$color_bar_orientation,
      xLegend=input$x_legend_pos,
      yLegend=input$y_legend_pos
      )
    )
    }
  
  return(NULL)
})

# Render maximum of derivative plot
output$tm_vs_lig <- renderPlotly({
  
  #req(fluo_signal_loaded())
  req(input$table1)
  if (is.null(modify_fluo_temp_cond())) {return(NULL)}
  
  p <- generate_der_plot(
    dsf$tmsFromDerivative,
    dsf$concentrations,
    plot_width=input$plot_width,
    plot_height=input$plot_height,
    plot_type=input$plot_type,
    axis_size=input$plot_axis_size,
    show_x_grid=input$show_x_grid,
    show_y_grid=input$show_y_grid,
    show_axis_lines=input$show_axis_lines,
    tickwidth=input$plot_tickwidth,
    ticklen=input$plot_ticklen,
    marker_size=input$plot_marker_size,
    log_scale_type=input$logScaleType)

  return(p)
  
}
)

# Render initial signal versus ligand concentration plot
output$initialFluo_vs_lig <- renderPlotly({
  
  req(input$table1)
  if (is.null(modify_fluo_temp_cond())) {return(NULL)}
  
  fluo_m <- make_df4plot(dsf$fluo,dsf$concentrations,dsf$temps)
  if (!(is.null(fluo_m))) {
    return( plot_initialFluo_signal(
      fluo_m = fluo_m,
      which = dsf$signal_type,
      plot_width = input$plot_width,
      plot_height = input$plot_height,
      plot_type = input$plot_type,
      axis_size = input$plot_axis_size,
      show_x_grid= input$show_x_grid,
      show_y_grid= input$show_y_grid,
      show_axis_lines = input$show_axis_lines,
      tickwidth = input$plot_tickwidth,
      ticklen = input$plot_ticklen,
      marker_size = input$plot_marker_size,
      log_scale_type = input$logScaleType
      )
    )
  }
  
  return(NULL)
})

# Fit when the user presses the button
fluo_fit_data <- eventReactive(input$btn_cal, {
  
  req(input$table1)  
  dsf$cp <- input$delta_cp
  
  if ( (max(dsf$temps) - min(dsf$temps)) >=50) {
    shinyalert(paste("Please select a temperature range that covers less than 50 degrees"),
             imageWidth = 180,imageHeight = 180,closeOnClickOutside=TRUE,closeOnEsc=TRUE,
             type="warning")
    
  return(NULL)
  }
  
  # check we have data to fit
  if (ncol(dsf$fluo)>0)   {
    
    withBusyIndicatorServer("btn_cal",{
      
      # ... Fit according to selection ...
      dsf$fit_fluo_local()
      if (input$model_selected != "Local")     { dsf$fit_fluo_global()    }
      if (input$model_selected == "Global_CP") { dsf$fit_fluo_global_cp() }
      
    })
    
    if (input$model_selected == "Global_CP" & dsf$fitting_option == "Global") {
      shinyalert(paste("Number of iterations excedeed for Global_CP. Using Global option."),
                 imageWidth = 180,imageHeight = 180,closeOnClickOutside=TRUE,closeOnEsc=TRUE,
                 type="info")
    } 
    
    # ... export fluorescence fit data ...
    dsf$shiny_export_fit_fluo()
    
    fluo_m      <- make_list_df4plot(dsf$fluo,dsf$concentrations,dsf$temps,global_chunck_n)
    fluo_m_pred <- make_list_df4plot(dsf$fit_fluo_pred,dsf$concentrations,dsf$temps,global_chunck_n)
    
    tms <- data.frame(t(dsf$fit_fluo_params))[,1]
    
    maxT <- ceiling(max(tms))+2
    minT <- floor(min(tms))-2
    
    updateSliderInput(session, "if_range", value = c(minT,maxT),
                      min = minT, max = maxT, step = 1)
    
    #if_range
    
    return(list("fluo_fit_real"=fluo_m,"fluo_fit_pred"=fluo_m_pred,
                "fit_params"=dsf$fit_fluo_params,"fit_errors"=dsf$fit_fluo_errs))
  } else {
    return(NULL)
  }
})

output$data_was_fitted <- reactive({
  # re evaluate expression if the user fitted data again
  req(input$btn_cal)
  req(fluo_fit_data())
  return(TRUE)
  
})

outputOptions(output, "data_was_fitted", suspendWhenHidden = FALSE)

# Change tables and plot accordingly - fluorescence fit
observeEvent(input$btn_cal,{
  
  output$params_table <- renderTable({
    req(fluo_fit_data())
    get_params_df(fluo_fit_data()$fit_params,dsf$concentrations)
  })
  
  output$params_table_errors <- renderTable({
    req(fluo_fit_data())
    get_params_errors_df(fluo_fit_data()$fit_params,dsf$concentrations,fluo_fit_data()$fit_errors)
  })
  
  # Plot the fluorescence fits
  output$fluo_fit_plot <- renderPlot({
    
    req(fluo_fit_data())
    req(input$select_fitting_plot)
    
    selected <- get_selected_from_choice_label(input$select_fitting_plot,global_chunck_n)
    
    real_data  <- fluo_fit_data()$fluo_fit_real
    model_data <- fluo_fit_data()$fluo_fit_pred
    
    return( plot_fluorescence_fit(real_data,model_data,selected) )
  })
  
})


ist_df_real <- reactive({
  
  # Check we have fitted the fluorescence data
  req(fluo_fit_data())
  
  # ... Update isothermals to show ...
  ts <- get_ts(input$if_range[1],input$if_range[2],min(dsf$temps),max(dsf$temps))
  if (length(ts) == 1) {ts <- np_array(ts)}
  dsf$isothermal_ts <- ts
  
  ## Put data in dataframe for plotting. Fraction unfolded versus ligand concentration 
  dsf$pre_fit_isothermal()
  format_ist_data_exp(dsf$isothermal_data,dsf$concentrations,dsf$isothermal_ts)
})

# Fit the isothermals when the user presses the button
ist_fit_data_pred <- eventReactive(input$btn_cal_fit_its, {

  if ( length(unique(dsf$concentrations)) <= 4) {
    shinyalert(paste("Please include information about the ligand concentration in the '1.Load input' tab"),
               imageWidth = 180,imageHeight = 180,closeOnClickOutside=TRUE,closeOnEsc=TRUE,
               type="warning")
    
    return(NULL)
  }
  
  dsf$pconc <- input$protein_conc / (1e6)  # convert protein concentration to molar units
  
  # ... Update isothermals to fit ...
  ts <- get_ts(input$if_range[1],input$if_range[2],min(dsf$temps),max(dsf$temps))
  if (length(ts) == 1) {ts <- np_array(ts)}
  dsf$isothermal_ts <- ts
  
  # ... Perform the fitting and export the data ...
  dsf$fit_isothermal(input$model_sites)
  dsf$shiny_export_isothermal()
  
  ## Put data in dataframe for plotting
  ist_df_data <- format_ist_data_exp_and_pred(dsf$isothermal_data,dsf$concentrations,dsf$isothermal_ts,
                                              dsf$kd_models,dsf$kd_models_lower,dsf$kd_models_upper,
                                              dsf$kd_model_conc,dsf$bind_params,dsf$bind_errors,
                                              dsf$bind_ci95_asymmetric_low, dsf$bind_ci95_asymmetric_up)
  
  return(ist_df_data)
  
})

#Plot the isothermals
output$isf <- renderPlotly({
  
  req(ist_df_real()) # Check we have isothermal data
  
  return(
    plot_isothermal_fitting_exp(
      ist_df_real(),
      plot_width = input$plot_width_isf,
      plot_height = input$plot_height_isf,
      plot_type = input$plot_type_isf,
      axis_size = input$plot_axis_size_isf,
      legend_size = input$plot_legend_size_isf,
      show_x_grid = input$show_x_grid_isf,
      show_y_grid = input$show_y_grid_isf,
      show_axis_lines = input$show_axis_lines_isf,
      tickwidth = input$plot_tickwidth_isf,
      ticklen = input$plot_ticklen_isf,
      marker_size = input$plot_marker_size_isf,
      log_scale_type = input$logScaleType_isf
    )
  )
})

# Change plot accordingly - isothermal fit
observeEvent(input$btn_cal_fit_its,{

  #Adapt the isothermal plot
  output$isf <- renderPlotly({
    
    req(ist_fit_data_pred()) # Check we have isothermal data

    return(
      plot_isothermal_fitting(
        ist_fit_data_pred()$data_exp,
        ist_fit_data_pred()$data_pred,
        plot_width = input$plot_width_isf,
        plot_height = input$plot_height_isf,
        plot_type = input$plot_type_isf,
        axis_size = input$plot_axis_size_isf,
        legend_size = input$plot_legend_size_isf,
        show_x_grid = input$show_x_grid_isf,
        show_y_grid = input$show_y_grid_isf,
        show_axis_lines = input$show_axis_lines_isf,
        tickwidth = input$plot_tickwidth_isf,
        ticklen = input$plot_ticklen_isf,
        marker_size = input$plot_marker_size_isf,
        line_width = input$plot_line_width_isf,
        log_scale_type = input$logScaleType_isf
      )
    )
  })
    
})

tm_fit_data <- eventReactive(input$btn_cal_tm_fit, {
  
  req(input$table1)
  
  df_model <- get_tm_df(dsf$concentrations,dsf$tmsFromDerivative)
  
  if (input$model_selected_tm_shift == "one_site_tm_shift")        (fit_result <- fit_1site_tm_shift(df_model))
  if (input$model_selected_tm_shift == "two_site_one_kd_tm_shift") (fit_result <- fit_2sites_one_kd_tm_shift(df_model))
  if (input$model_selected_tm_shift == "two_site_two_kd_tm_shift") (fit_result <- fit_2sites_two_kd_tm_shift(df_model))
  if (input$model_selected_tm_shift == "three_site_one_kd_tm_shift")      (fit_result <- fit_3sites_one_kd_tm_shift(df_model))
  
  dsf$tms_fit_pred           <- as.list(fit_result$df_pred)  # So it is JSNON serializable
  dsf$tms_fit_info           <- as.list(fit_result$fit_info) # So it is JSNON serializable
  dsf$tms_fit_asymmetricCI95 <- fit_result$asymmetricCI95
  dsf$tms_fit_model          <- input$model_selected_tm_shift

  reactives$tm_fit_done <- TRUE
    
  return(fit_result)
  
})

output$tm_shift_fitting <- renderPlotly({
  
  req(input$table1)
  if (is.null(modify_fluo_temp_cond())) {return(NULL)}
  
  return(
    plot_tm_shift(
      dsf$concentrations,
      dsf$tmsFromDerivative,
      plot_width=input$plot_width_tmShift,
      plot_height=input$plot_height_tmShift,
      plot_type=input$plot_type_tmShift,
      axis_size=input$plot_axis_size_tmShift,
      show_x_grid=input$show_x_grid_tmShift,
      show_y_grid=input$show_y_grid_tmShift,
      show_axis_lines=input$show_axis_lines_tmShift,
      tickwidth=input$plot_tickwidth_tmShift,
      ticklen=input$plot_ticklen_tmShift,
      marker_size=input$plot_marker_size_tmShift,
      log_scale_type=input$logScaleType_tmShift
    )
  )
  
})

observeEvent(input$btn_cal_tm_fit,{

  output$tm_shift_fitting <- renderPlotly({
    
    req(input$table1)
    if (is.null(modify_fluo_temp_cond())) {return(NULL)}
    
    asymmetricCI95 <- NULL
    if (!(input$model_selected_tm_shift == "two_site_two_kd_tm_shift")){
      asymmetricCI95 <- tm_fit_data()$asymmetricCI95
    } 
    
    # Check that the user has pressed the fitting button more than 0 times
    p <- plot_tm_shift_fit(
      dsf$concentrations,
      dsf$tmsFromDerivative,
      tm_fit_data()$df_pred,
      tm_fit_data()$fit_info,
      asymmetricCI95,
      plot_width = input$plot_width_tmShift,
      plot_height = input$plot_height_tmShift,
      plot_type = input$plot_type_tmShift,
      axis_size = input$plot_axis_size_tmShift,
      legend_size = input$plot_legend_size_tmShift,
      show_x_grid = input$show_x_grid_tmShift,
      show_y_grid = input$show_y_grid_tmShift,
      show_axis_lines = input$show_axis_lines_tmShift,
      tickwidth = input$plot_tickwidth_tmShift,
      ticklen = input$plot_ticklen_tmShift,
      marker_size = input$plot_marker_size_tmShift,
      line_width = input$plot_line_width_tmShift,
      log_scale_type = input$logScaleType_tmShift
    )

    return( p )
  })
    
})

output$params_table_tm_shift <- renderTable({
  req(tm_fit_data())
  df <- tm_fit_data()$fit_info[,1:3]
  df <- format_tm_fit_df(df,dsf$tms_fit_model,dsf$tms_fit_asymmetricCI95)
  
  return(df)
})

# Simulations plots

output$kd_vs_temp <- renderPlotly({
  
  if (!(input$runSimulation)) (return(NULL)) # proceed only if the user checks the button 'Launch simulation'
  
  temp_seq <- seq(298,298+75,2)
  if (input$const_kb) {
    kd_seq   <- sapply(temp_seq, function(t){input$kd_simulation_const * 1e-6}) # * 1e-6 to convert micromolar to molar
  } else {
    kd_seq   <- sapply(temp_seq, function(t){get_kd(t,input$dHb,input$dCPb,input$dSb,input$Tb)})
  }
  
  plot_k_seq(data.frame(k=kd_seq,temp=temp_seq),"Kd (M)")
  
})

output$ku_vs_temp <- renderPlotly({
  
  if (!(input$runSimulation)) (return(NULL)) # proceed only if the user checks the button 'Launch simulation'
  
  temp_seq <- seq(300,360,2)
  ku_seq   <- sapply(temp_seq, function(t){get_ku_sim(t,input$dHu,input$dCPu,input$dHu/(input$Tu+273.15),input$Tu)})
  plot_k_seq( data.frame("k"=ku_seq,"temp"=temp_seq) ,"Ku")
  
})

fractions_and_fluo_vs_total_ligand <- reactive({
  
  if (!(input$runSimulation)) (return(NULL)) # proceed only if the user checks the button 'Launch simulation'
  
  p0         <- input$pconc_sim               / 1e6   # From micromolar to molar
  maxLigConc <- input$maxLigandConcSimulation / 1e6   # From micromolar to molar
  
  if (input$const_kb) {
    df <- get_fractions_df_from_thermodynamic_params2_const_kd(
      p0,input$b1,input$b2,input$b3,input$m1,input$m2,input$m3,input$kd_simulation_const / 1e6,
      input$dHu,input$dCPu,input$dHu/(input$Tu+273.15),input$Tu,maxLigConc)
    
  } else {
    
    df <- get_fractions_df_from_thermodynamic_params2(
      p0,input$b1,input$b2,input$b3,input$m1,input$m2,input$m3,input$dHb,input$dCPb,
      input$dSb,input$Tb,input$dHu,input$dCPu,input$dHu/(input$Tu+273.15),input$Tu,maxLigConc)
  }
  
  return(df)
  
})

output$fractions_vs_total_ligand <- renderPlot({

  if (!(input$runSimulation)) (return(NULL)) # proceed only if the user checks the button 'Launch simulation'
  return(plot_fraction_versus_total_ligand(fractions_and_fluo_vs_total_ligand()$melt2))
  
})

output$fluo_vs_ligand <- renderPlot({

  if (!(input$runSimulation)) (return(NULL)) # proceed only if the user checks the button 'Launch simulation'
  return(plot_fluo_vs_ligand(fractions_and_fluo_vs_total_ligand()$tog))
  
})

observeEvent(input$runSimulation,{
  
  if (input$runSimulation) {
    shinyalert(paste("Please take into account that according to this model, 
                     the Tm could be increased as much as desired. 
                     This doesn't happen in practice due to the fact that irreversible effects also take place
                     , and that there exists a solubility limit for the ligand."),
               imageWidth = 400,imageHeight = 400,closeOnClickOutside=TRUE,closeOnEsc=TRUE,
               type="warning")
  }

})

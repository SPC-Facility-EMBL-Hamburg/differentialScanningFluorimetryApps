resetConditionsTable <- function() {

    output$table1 <- NULL
    output$table2 <- NULL
    output$table3 <- NULL
    output$table4 <- NULL

    return(NULL)

}

renderConditionsTable <- function(tables) {

    lapply(seq_along(tables), function(i) {
        output[[paste0("table", i)]] <- tables[[i]]
        return(NULL)
    })

    return(NULL)

}

observeEvent(input$dsf_files,{

    reactives$data_loaded <- FALSE
    reactives$data_was_fitted <- FALSE
    reactives$have_plot_data <- FALSE

    resetConditionsTable()

    dsf$delete_experiments("ALL")

    dsf_data_files <- input$dsf_files$datapath
    names <- input$dsf_files[[1]]

    sorted_indices <- order(names)

    dsf_data_files <- dsf_data_files[sorted_indices]
    names <- names[sorted_indices]

    withBusyIndicatorServer("Go",{

        i <- 0
        # iterate over the files
        for (name in names) {

            i <- i + 1
            dsf_data_file <- dsf_data_files[i]

            read_file_status <- dsf$add_experiment(dsf_data_file,name)

            if (!read_file_status) {

                shinyalert(
                    title = "Error",
                    text = "The selected file could not be read. Please, select a valid DSF data file.",
                    type = "error"
                )

                return(NULL)
            }
        }

        # Update the wavelength range slider for full spectra data
        reactives$full_spectra <- dsf$full_spectrum

        if (dsf$full_spectrum) {
            conditions <- dsf$get_experiment_properties('conditions',flatten=TRUE,mode="full_spectrum")
            reactives$min_wl <- dsf$min_wavelength
            reactives$max_wl <- dsf$max_wavelength

        }

        conditions <- dsf$get_experiment_properties('conditions',flatten=TRUE,mode="all")

        reactives$nconditions <- length(conditions)

        updateSelectInput(session, "which",choices = dsf$all_signals)

        dsf$set_signal(dsf$all_signals[1])

        reactives$global_max_conditions          <- findClosestHigherValue(reactives$nconditions)
        reactives$global_n_rows_conditions_table <- reactives$global_max_conditions / 4

        tables <- get_renderRHandsontable_list(
          conditions,
          reactives$global_n_rows_conditions_table
        )

        renderConditionsTable(tables)

        temps <- dsf$get_experiment_properties('temps',flatten=TRUE,mode="all")

        min_temp <- round(min(temps) - 273.15) # To degree celsius
        max_temp <- round(max(temps) - 273.15) # To degree celsius

        updateSliderInput(
          session,"sg_range",
          NULL,
          min = min_temp,
          max = max_temp,
          value = c(min_temp+1,max_temp-1)
        )

        dsf$estimate_fluo_derivates(input$SG_window2)

        Sys.sleep(0.5)
        reactives$data_loaded <- TRUE

    })

},priority = 10)

observeEvent(
  list(
    input$table1,input$table2,input$table3,input$table4,
    input$which,input$sg_range,input$median_filter,input$SG_window2,
    input$normalization_type,input$selected_cond_series),{

  req(reactives$data_loaded)

  condition_include_list <- get_include_vector(
    input$table1,input$table2,input$table3,input$table4,
    reactives$nconditions,reactives$global_n_rows_conditions_table,
    reactives$global_max_conditions)

  conditions_vector <- as.character(condition_include_list$conditions_vector)
  include_vector <- as.logical(condition_include_list$include_vector)
  series_vector <- condition_include_list$series_vector

  color_vector <- as.character(condition_include_list$color_vector)

  if (input$selected_cond_series != "ALL") {
    include_vector <- include_vector & (series_vector == input$selected_cond_series)
  }

  # Return NULL if no conditions are selected
  if (all(!include_vector))  {
    reactives$have_plot_data <- FALSE
    return(NULL)
  } else {
    reactives$have_plot_data <- TRUE
  }

  # Get signal window range
  left_bound <- input$sg_range[1]
  right_bound <- input$sg_range[2]

  # Give an error if the difference is less than 13 degrees
  if ( (right_bound - left_bound) < 13 ) {
    shinyalert(
      title = "Error",
      text = "The selected temperature range is too small. Please, select a range of at least 13 degrees Celsius.",
      type = "error"
    )
    reactives$have_plot_data <- FALSE
    return(NULL)
  }

  reactives$data_loaded <- FALSE

  dsf$set_signal(input$which)

  # ... Modify in place the python class fluorescence signal according to the selected signal window range ...
  dsf$filter_by_temperature(left_bound,right_bound)

  # ... use only the conditions selected by the user ...
  dsf$set_conditions(conditions_vector)

  dsf$set_colors(color_vector)

  dsf$select_conditions(include_vector)

  median_filter <- get_median_filter(input$median_filter)
  if (median_filter > 0) {dsf$median_filter(median_filter)}

  exps <- dsf$available_experiments

  for (exp in exps) {

    py_obj <- dsf$experiments[[exp]]
    py_obj$fluo <- normalize_fluo_matrix_by_option(input$normalization_type,py_obj$fluo,py_obj$temps)

  }

  dsf$estimate_fluo_derivates(input$SG_window2)

  current_series <- c(input$selected_cond_series)

  for (series in series_vector) {
    if (!(series %in%  current_series)) {
      # Add new series found in the table
      current_series <- c(current_series,series)
    }
  }

  for (series in current_series) {
    if (!(series %in%  series_vector)) {
      # Delete series no longer in the table, except for ALL
      if (series == 'ALL') {next}
      current_series <- current_series[current_series != series]
    }
  }

  # Add the "ALL" option to the end, if not present
    if (!('ALL' %in% current_series)) {
        current_series <- c(current_series,'ALL')
    }


  updateSelectInput(session, "selected_cond_series",choices  = c(current_series))

  reactives$data_loaded <- TRUE

},ignoreInit = TRUE)

observeEvent(input$layout_file$datapath,{
  
  req(input$table1)
  if (!(is.null(input$layout_file))) {
    
    new_conditions <- load_layout(input$layout_file$datapath)[1:reactives$nconditions]

    dsf$set_conditions(new_conditions)
    dsf$set_conditions(new_conditions,original=TRUE)

    tables <- get_renderRHandsontable_list(
      new_conditions,
      reactives$global_n_rows_conditions_table
    )
    
    renderConditionsTable(tables)

  }
  
})

observeEvent(input$sort_conditions,{
  req(input$table1)
  
  resetConditionsTable()

  dsf$sort_by_conditions_name(input$sort_conditions)

  conditions_original <- dsf$get_experiment_properties('conditions_original',flatten=TRUE,mode="all")

  tables <- get_renderRHandsontable_list(
    conditions_original,
    reactives$global_n_rows_conditions_table
  )
  
  renderConditionsTable(tables)
  
})

observeEvent(input$show_colors_column,{

    req(input$table1)

    show_colors_column <- input$show_colors_column

    condition_include_list <- get_include_vector(
        input$table1,input$table2,input$table3,input$table4,
        reactives$nconditions,reactives$global_n_rows_conditions_table,
        reactives$global_max_conditions)

    conditions_vector      <- as.character(condition_include_list$conditions_vector)
    include_vector         <- as.logical(condition_include_list$include_vector)
    series_vector          <- as.character(condition_include_list$series_vector)

    # If the colors are now set to be hidden, we obtain them from the Tables
    if (!show_colors_column){
        color_vector <- as.character(condition_include_list$color_vector)
    } else {
    # If the colors are shown, we use the colors from the dsf object
        color_vector <- dsf$get_experiment_properties('all_colors',flatten=TRUE,mode="all")
    }

    tables <- get_renderRHandsontable_list(
        conditions_vector,
        reactives$global_n_rows_conditions_table,
        colors=color_vector,
        series=series_vector,
        include=include_vector,
        hide_color_column=!show_colors_column)

    renderConditionsTable(tables)

})

# Render signal plot
output$signal <- renderPlotly({
  
  req(input$table1)
  req(reactives$data_loaded)

  if (!reactives$have_plot_data) {
    return(NULL)
  }

  fluo_m <- py_dsf_to_df(dsf)

  colors <- dsf$get_experiment_properties('colors',flatten=TRUE)

  p <- plot_fluo_signal(fluo_m,colors,input$which,
                        input$plot_width, input$plot_height,
                        input$plot_type,input$plot_font_size,input$plot_axis_size,
                        show_x_grid=input$show_x_grid,
                        show_y_grid=input$show_y_grid,
                        show_axis_lines=input$show_axis_lines,
                        tickwidth=input$plot_tickwidth,
                        ticklen=input$plot_ticklen,
                        line_width=input$plot_line_width)
  return(p)


}
)


# Render 1st derivative plot
output$signal_der1 <- renderPlotly({
  
  req(input$table1)
  req(reactives$data_loaded)

  if (!reactives$have_plot_data) {
    return(NULL)
  }

  max_cond <- reactives$global_max_conditions
  n_rows   <- reactives$global_n_rows_conditions_table

  if (reactives$nconditions > max_cond - n_rows*1) {
    req(input$table4)
  }

  fluo_m <- py_dsf_to_df(dsf,mode = 'derivative')
  colors <- dsf$get_experiment_properties('colors',flatten=TRUE)

  p <- plot_fluo_signal(fluo_m,colors,"First derivative",
                        input$plot_width, input$plot_height,
                        input$plot_type,input$plot_font_size,input$plot_axis_size,
                        show_x_grid=input$show_x_grid,
                        show_y_grid=input$show_y_grid,
                        show_axis_lines=input$show_axis_lines,
                        tickwidth=input$plot_tickwidth,
                        ticklen=input$plot_ticklen,
                        line_width=input$plot_line_width)
  return(p)


}
)

# Render 2nd derivative plot
output$signal_der2 <- renderPlotly({
  
  req(input$table1)
  req(reactives$data_loaded)

  if (!reactives$have_plot_data) {
    return(NULL)
  }

  fluo_m <- py_dsf_to_df(dsf,mode = 'derivative2')
  colors <- dsf$get_experiment_properties('colors',flatten=TRUE)
  
  p <- plot_fluo_signal(fluo_m,colors,"Second derivative",
                        input$plot_width, input$plot_height,
                        input$plot_type,input$plot_font_size,input$plot_axis_size,
                        show_x_grid=input$show_x_grid,
                        show_y_grid=input$show_y_grid,
                        show_axis_lines=input$show_axis_lines,
                        tickwidth=input$plot_tickwidth,
                        ticklen=input$plot_ticklen,
                        line_width=input$plot_line_width)
  return(p)


}
)

# Render maximum of derivative plot
output$tm_derivative <- renderPlotly({
  
  req(input$table1)
  req(reactives$data_loaded)

  # Require that at least one derivative is available
  derivatives <- dsf$get_experiment_properties('derivative')
  if (all(sapply(derivatives, is.null))) {
    return(NULL)
  }

  conditions <- dsf$get_experiment_properties('conditions',flatten=TRUE)
  tms_from_deriv <- dsf$get_experiment_properties('tms_from_deriv',flatten=TRUE)

  p <- generate_max_der_plot(
    tms_from_deriv,
    conditions,
    input$plot_width,
    input$plot_height,
    input$plot_type,
    input$plot_font_size,
    input$plot_axis_size
  )
  
  return(p)
  
})

# Fit when the user presses the button
observeEvent(input$btn_cal, {

  model_selected <- input$model_selected

  orders <- list(
    "constant"  = 0,
    "linear"    = 1,
    "quadratic" = 2
  )

  poly_order_native <- orders[[input$native_dependence]]
  poly_order_unfolded <- orders[[input$unfolded_dependence]]

  dsf$set_baseline_types(poly_order_native,poly_order_unfolded)

  req(input$table1)
  req(reactives$have_plot_data)

  reactives$data_was_fitted <- FALSE
  reactives$fluo_fit_real <- NULL
  reactives$fluo_fit_pred <- NULL

  dsf$estimate_baselines_parameters(input$temp_range_baseline_estimation)

  reactives$model_is_two_state <- grepl('TwoState',input$model_selected)
  reactives$model_name <- input$model_selected

  withBusyIndicatorServer("hiddenBtnFit",{

    # ... Fit according to selection ...
    if (model_selected == "EquilibriumTwoState"  )   {

      dsf$equilibrium_two_state(input$delta_cp)

    }

    if (model_selected == "EquilibriumThreeState")   {

      dsf$equilibrium_three_state(
        input$t1min,input$t1max,
        input$t2min,input$t2max
      )
    }

    if (model_selected == "EmpiricalTwoState")  dsf$empirical_two_state()

    if (model_selected == "EmpiricalThreeState")   {

      dsf$empirical_three_state(
        input$t1min,input$t1max,
        input$t2min,input$t2max
      )

    }

    if (model_selected == "IrreversibleTwoState" )   {
      dsf$irreversible_two_state(input$scan_rate)
    }

  })

  fitted_conditions <- dsf$get_experiment_properties('fitted_conditions',flatten=TRUE)

  if (length(fitted_conditions) == 0 ) {
    shinyalert(
      title = "Error",
      text = "No conditions were fitted. Please, check that the selected model is appropriate for your data, and try changing the fitting parameters.",
      type = "error"
    )
    return(NULL)
  }

  fluo_m <- make_list_df4plot(dsf,global_chunck_n,mode="experimental")

  reactives$data_was_fitted <- TRUE

  std_error_estimate_all <- dsf$get_experiment_properties('std_error_estimate_all',flatten=TRUE)
  max_std_err <- max(std_error_estimate_all)

  updateNumericInput(
      session,
      "fitting_std_threshold",
      value = signif(max_std_err*1.01,3),
      min   = 0,
      max   = max_std_err*1.02,
      step  = signif(max_std_err/100,3)
  )

  fluo_m_pred <- make_list_df4plot(dsf,global_chunck_n,mode='fitted')

  reactives$fluo_fit_real <- fluo_m
  reactives$fluo_fit_pred <- fluo_m_pred

  
})

output$three_state_model_selected <- reactive({
  return( grepl("Three",input$model_selected) ) 
})

observeEvent(input$model_selected,{
  
  req(input$table1)
  
  if (grepl("Three",input$model_selected)) {

    temps <- dsf$get_experiment_properties('temps',flatten=TRUE) - 273.15
    tempMax <- round(max(temps))
    tempMin <- round(min(temps))
    tempMiddle <- round( (max(temps) + min(temps)) * 0.5 )
    
    updateNumericInput(session, "t1max", value = tempMiddle+6,  min = tempMin, max = tempMax)
    updateNumericInput(session, "t1min", value = tempMin+5,   min = tempMin, max = tempMax)
    updateNumericInput(session, "t2max", value = tempMax-5,   min = tempMin, max = tempMax)
    updateNumericInput(session, "t2min", value = tempMiddle-6,  min = tempMin, max = tempMax)
    
  }
  
  return(NULL)
})


outputOptions(output, "three_state_model_selected", suspendWhenHidden = FALSE)

output$params_table <- renderTable({

  req(reactives$data_was_fitted)

  params_all <- dsf$get_experiment_properties('params_all',flatten=TRUE)
  params_name <- dsf$get_experiment_properties('params_name',flatten=TRUE)
  params_name <- unique(params_name)

  fitted_conditions <- dsf$get_experiment_properties('fitted_conditions',flatten=TRUE)

  get_sorted_params_table(
    params_all,
    fitted_conditions,
    global_chunck_n,
    params_name,
    input$sort_table_parameter
  )

})

output$params_table_errors <- renderTable({

  req(reactives$data_was_fitted)
  errors_percentage_all <- dsf$get_experiment_properties('errors_percentage_all',flatten=TRUE)
  fitted_conditions <- dsf$get_experiment_properties('fitted_conditions',flatten=TRUE)
  params_name <- dsf$get_experiment_properties('params_name',flatten=TRUE)
  params_name <- unique(params_name)

  return(get_sorted_params_table_errors(
    errors_percentage_all,
    fitted_conditions,
    params_name,
    input$sort_table_parameter)
  )

})

observe({
  
  req(reactives$data_was_fitted)
  fitted_conditions <- dsf$get_experiment_properties('fitted_conditions',flatten=TRUE)

  updateSelectInput(
    session, "select_fitting_plot",
    choices  = get_choices_fluo_fits(length(fitted_conditions),global_chunck_n)
  )
  
})

observe({
  
  req(reactives$data_was_fitted)

  params_name <- dsf$get_experiment_properties('params_name',flatten=TRUE)
  params_name <- unique(params_name)

  updateSelectInput(
    session, "sort_table_parameter",
    choices  = params_name
  )
  
})

# Plot the fluorescence fits
output$fluo_fit_plot <- renderPlot({
  
  req(reactives$data_was_fitted)
  req(input$select_fitting_plot)
  
  selected <- get_selected_from_choice_label(input$select_fitting_plot,global_chunck_n)
  
  real_data  <- reactives$fluo_fit_real
  model_data <- reactives$fluo_fit_pred

  p <- plot_fluorescence_fit(real_data,model_data,selected)
  return(p)
}
)

# Plot the fluorescence fits standarized residuals
output$fluo_residuals_plot <- renderPlot({
  
  req(reactives$data_was_fitted)
  req(input$select_fitting_plot)
  
  selected <- get_selected_from_choice_label(input$select_fitting_plot,global_chunck_n)
  
  real_data  <- reactives$fluo_fit_real
  model_data <- reactives$fluo_fit_pred

  std_error_estimate_all <- dsf$get_experiment_properties('std_error_estimate_all',flatten=TRUE)
  fitted_conditions <- dsf$get_experiment_properties('fitted_conditions',flatten=TRUE)

  p <- plot_fluorescence_residuals(
    real_data,
    model_data,
    selected,
    std_error_estimate_all,
    fitted_conditions
  )

  return(p)
}
)

# End of Plot the fluorescence fits

output$fitted_conditions_table <- renderTable({

  req(reactives$data_was_fitted)

  conditions <- dsf$get_experiment_properties('conditions',flatten=TRUE)
  fitted_conditions <- dsf$get_experiment_properties('fitted_conditions',flatten=TRUE)

  return(
    get_fitted_conditions_table(
    conditions,
    fitted_conditions)
    )

})


filter_conditions <- reactive({
  
  req(reactives$data_was_fitted)
  fitted_conditions <- dsf$get_experiment_properties('fitted_conditions',flatten=TRUE)

  if (input$sd_factor_bool) {

    selected_indexes_1 <- dsf$filter_by_relative_error(input$rel_error_threshold)

  } else {

      selected_indexes_1 <- rep(TRUE,length(fitted_conditions))

  }

  selected_indexes_2 <- dsf$filter_by_fitting_std_error(input$fitting_std_threshold)

  selected_indexes <- unlist(selected_indexes_1) & unlist(selected_indexes_2)

  if (reactives$model_name %in% c("EmpiricalTwoState","EquilibriumTwoState")) {

    selected_indexes_to_add <- dsf$filter_by_param_values(
        'Tm',
        input$lower_Tm_threshold+273.15, # to kelvin
        input$upper_Tm_threshold+273.15
    )

    selected_indexes <- selected_indexes & unlist(selected_indexes_to_add)

  }

  if (reactives$model_name %in% c("EquilibriumTwoState")) {

    selected_indexes_to_add <- dsf$filter_by_param_values(
        'DH',
        input$lower_dh_threshold,
        input$upper_dh_threshold
    )

    selected_indexes <- selected_indexes & unlist(selected_indexes_to_add)

  }

  if (reactives$model_name %in% c("EmpiricalThreeState","EquilibriumThreeState")) {

    selected_indexes_to_add <- dsf$filter_by_param_values(
        'T1',
        input$lower_T1_threshold+273.15, # to kelvin
        input$upper_T1_threshold+273.15
    )

    selected_indexes <- selected_indexes & unlist(selected_indexes_to_add)

    selected_indexes_to_add <- dsf$filter_by_param_values(
        'T2',
        input$lower_T2_threshold+273.15, # to kelvin
        input$upper_T2_threshold+273.15
    )

    selected_indexes <- selected_indexes & unlist(selected_indexes_to_add)

  }

  return(selected_indexes)

})

get_score_table <- reactive({
  
  # Check we have fitted the fluorescence data
  req(reactives$data_was_fitted)

  fitted_conditions <- dsf$get_experiment_properties('fitted_conditions',flatten=TRUE)

  selected_indexes <- filter_conditions()

  if (!(any(selected_indexes))){return(NULL)}
  
  score_table <- NULL
  if ((dsf$model_name == "EquilibriumTwoState")) {

    dG_std <- dsf$get_experiment_properties('dG_std',flatten=TRUE)
    dCp_component <- dsf$get_experiment_properties('dCp_component',flatten=TRUE)

    params_all  <- map2(dG_std,dCp_component,function(x,y) c(x,y))
    params_name <- c("dg_std","cp_comp")
    
    score_table <- get_all_params_df(
      params_all[selected_indexes],
      fitted_conditions[selected_indexes],
      global_chunck_n,params_name
    )
    
    score_table$dg_std <- as.numeric(score_table$dg_std) 
    score_table        <- score_table %>% dplyr::arrange(dg_std)
    colnames(score_table) <- c("ΔG_unfolding (25°C) (kcal/mol) ","ΔCp_component(K)","Condition")
  } 
  
  if ((dsf$model_name == "EquilibriumThreeState")) {
    
    params_name <- c("score")
    dG_comb_std <- dsf$get_experiment_properties('dG_comb_std',flatten=TRUE)

    score_table <- get_all_params_df(
      dG_comb_std[selected_indexes],
      fitted_conditions[selected_indexes],
      global_chunck_n,params_name
    )
    
    score_table$score    <-  as.numeric(score_table$score) 
    score_table          <-  score_table %>% dplyr::arrange(score)
    colnames(score_table) <- c("ΔG_comb (25°C) (kcal/mol)","Condition")
    
  }
  
  if ((dsf$model_name == "EmpiricalTwoState")) {
    
    params_name <- c("score")

    score <- dsf$get_experiment_properties('score',flatten=TRUE)

    score_table <- get_all_params_df(
      score[selected_indexes],
      fitted_conditions[selected_indexes],
      global_chunck_n,params_name
    )
    
    score_table$score    <-  as.numeric(score_table$score) 
    score_table          <-  score_table %>% dplyr::arrange(-score)
    colnames(score_table) <- c("sqrt( T_onset**2 + Tm**2 )","Condition")
    
  }
  
  if ((dsf$model_name == "EmpiricalThreeState")) {
    
    params_name <- c("score")
    T_eucl_comb <- dsf$get_experiment_properties('T_eucl_comb',flatten=TRUE)

    score_table <- get_all_params_df(
      T_eucl_comb[selected_indexes],
      fitted_conditions[selected_indexes],
      global_chunck_n,params_name)
    
    score_table$score    <-  as.numeric(score_table$score) 
    score_table          <-  score_table %>% dplyr::arrange(-score)
    colnames(score_table) <- c("sqrt( T_onset_1**2 + T1**2 ) + sqrt( T_onset_1**2 + T2**2 )","Condition")
    
  }
  
  if ((dsf$model_name == "IrreversibleTwoState")) {
    
    params_name <- c("score")
    pkd <- dsf$get_experiment_properties('pkd',flatten=TRUE)

    score_table <- get_all_params_df(
      pkd[selected_indexes],
      fitted_conditions[selected_indexes],
      global_chunck_n,params_name
    )
    
    score_table$score    <-  as.numeric(score_table$score) 
    score_table          <-  score_table %>% dplyr::arrange(-score)
    colnames(score_table) <- c("pkd","Condition")
    
  }
  
  if (is.null(score_table)) {return(NULL)} 
  return(score_table)
  
})

output$score_table <- renderTable({
  
  req(get_score_table())
  return(get_score_table())
  
})

observe({
  
  req(reactives$data_was_fitted)
  updateSelectInput(
    session,
    "select_plot_type",
    choices  = get_choices_result_plot(dsf$model_name)
    )
  
})


get_results_plot <- reactive({
  
  selected_indexes <- filter_conditions()
  
  if (!(any(selected_indexes))){return(NULL)}

  params_name <- dsf$get_experiment_properties('params_name',flatten=TRUE)
  params_name <- unique(params_name)

  params_all <- dsf$get_experiment_properties('params_all',flatten=TRUE)

  fitted_conditions <- dsf$get_experiment_properties('fitted_conditions',flatten=TRUE)

  # Get the value of the parameter Tm
  if (input$select_plot_type %in% c(
    "Unfolded fraction",
    "The 25 highest Tms",
    "The 25 highest Tms versus Tonset")) {
    
    tm_position <- which(params_name == "Tm")
    tms <- sapply(params_all, function(x) x[tm_position])
    
  }
  
  #Get the value of the parameter T_onset 
  if (input$select_plot_type %in% c("Unfolded fraction","The 25 highest Tms versus Tonset")) {
    
    if (input$model_selected == "EquilibriumTwoState") {
      t_onset <- dsf$get_experiment_properties('T_onset',flatten=TRUE)
    }
    if (input$model_selected == "EmpiricalTwoState")   {
      t_onset_position <- which(params_name == "T_onset")
      t_onset          <- sapply(params_all, function(x) x[t_onset_position])
    }
  }
  
  if (input$select_plot_type == "Unfolded fraction") {
    
    dh_position <- which(params_name == "DH")
    dhs <- sapply(params_all, function(x) x[dh_position])
    
    if (input$model_selected == "EquilibriumTwoState") {

      fig <- generate_fractions_plot(
        dhs[selected_indexes],
        tms[selected_indexes],
        fitted_conditions[selected_indexes],
        t_onset[selected_indexes], # useful only to set colors
        input$plot_width_results, input$plot_height_results,
        input$plot_type_results,input$plot_font_size_results,
        input$plot_axis_size_results,input$fractions_plot_style
      )
    }
    
    return(fig)
  }
  
  if (input$select_plot_type == "The 25 highest Tms") {
    
    fig <- generate_tm_plot(
      tms[selected_indexes],
      fitted_conditions[selected_indexes],
      input$plot_width_results, input$plot_height_results,
      input$plot_type_results,input$plot_axis_size_results
    )

    return(fig)
  }
  
  if (input$select_plot_type == "The 25 highest Tms versus Tonset") {
    
    fig <- generate_tm_tonset_plot(
      tms[selected_indexes],
      fitted_conditions[selected_indexes],
      t_onset[selected_indexes],
      input$plot_width_results, input$plot_height_results,
      input$plot_type_results,input$plot_font_size_results,input$plot_axis_size_results
    )

    return(fig)
  }
  
  if (input$select_plot_type == "The 25 highest combinations of T1 and T2") {
    
    t1_position <- which(params_name == "T1")
    t1s <- sapply(params_all, function(x) x[t1_position])
    
    t2_position <- which(params_name == "T2")
    t2s <- sapply(params_all, function(x) x[t2_position])
    
    if (input$model_selected == "EquilibriumThreeState") {
      dG_comb_std <- dsf$get_experiment_properties('dG_comb_std',flatten=TRUE)
      score <- dG_comb_std
    }
    if (input$model_selected == "EmpiricalThreeState")   {
      T_eucl_comb <- dsf$get_experiment_properties('T_eucl_comb',flatten=TRUE)
      score <- T_eucl_comb
    }
    
    fig <- generate_2t_plot(
      t1s[selected_indexes],
      t2s[selected_indexes],
      fitted_conditions[selected_indexes],
      score[selected_indexes],
      input$plot_width_results, input$plot_height_results,
      input$plot_type_results,input$plot_font_size_results,input$plot_axis_size_results
    )

    return(fig)
    
  }

  if (input$select_plot_type == "Score versus condition") {
    
    req(get_score_table())
    return(generate_score_plot(
      get_score_table(),
      input$plot_width_results,input$plot_height_results,
      input$plot_type_results,input$plot_axis_size_results)
    )
    
  }
  
})

output$results_plot <- renderPlotly({
  
  req(reactives$data_was_fitted)
  req(input$select_plot_type)

  plot <- get_results_plot()
  return(plot)
  
})

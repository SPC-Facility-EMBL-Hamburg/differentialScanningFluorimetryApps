## Plot fluorescence signal versus temperature (color by ligand concentration)
## fluo_m is a 3 column dataframe: temp, fluo and conc. 
## which is the signal (350,330 or Ratio)
## Returns the plot

plot_fluo_signal <- function(
  fluo_m,
  which="Ratio",
  plot_width=12,
  plot_height=8,
  plot_type='png',
  legend_text_size=13,
  axis_size=14,
  show_x_grid=TRUE,
  show_y_grid=TRUE,
  show_axis_lines=TRUE,
  tickwidth=2,
  ticklen=8,
  line_width=2,
  color_bar_length=0.5,
  color_bar_orientation="h",
  xLegend=1,
  yLegend=1) {
  
  x <- list(
    title = "Temperature (°C)",
    titlefont = list(size = axis_size),
    tickfont = list(size = axis_size),
    showgrid = show_x_grid,
    showline = show_axis_lines,
    zeroline = FALSE,
    ticks = "outside",
    tickwidth = tickwidth*show_axis_lines,
    ticklen = ticklen*show_axis_lines,
    automargin = TRUE
    )

  y <- list(
    title = list(text=which,standoff=10),
    titlefont = list(size = axis_size),
    tickfont = list(size = axis_size),
    showgrid = show_y_grid,
    showline = show_axis_lines,
    zeroline = FALSE,
    ticks = "outside",
    tickwidth = tickwidth*show_axis_lines,
    ticklen = ticklen*show_axis_lines,
    automargin = TRUE
    )

  fluo_m0 <- fluo_m %>% filter(conc == 0)
  fluo_m  <- fluo_m %>% filter(conc > 0)

  if (nrow(fluo_m) == 0) { 

    fig <- plot_ly()
    fig <- fig %>% add_trace(data = fluo_m0 ,name = ~conc_, 
                             color = I("#389196"),
                             y    = ~fluo   , x   = ~temp)
    
    # if we ONLY have positions with no ligand return this plot
    return(fig) 
  }

  # Arrange fluo_m by concentration and temperature
  # sorts data in ascending order by default
  fluo_m <- fluo_m %>%
    arrange(conc) %>%
    arrange(temp)

  df_ligColorBar <- fluo_m[c(1,nrow(fluo_m)),]
  df_ligColorBar$conc <- log10(df_ligColorBar$conc)


  fig <- plot_ly(df_ligColorBar,
                 x = ~temp,
                 y = ~fluo,
                 color = ~conc,
                 type = "scatter",
                 mode = "markers",
                 showlegend = FALSE,
                 marker=list(size=0.01))

  for (position in unique(fluo_m0$conc_)) {

    fluo_m0_temp <- fluo_m0 %>% filter(conc_== position) 

    # set as black the positions with no ligand
    fig <- fig %>% add_trace(data  = fluo_m0_temp,
                             color = I("black"),
                             type='scatter',
                             mode = 'lines',
                             line = list(width = line_width),
                             showlegend=FALSE,
                             inherit= FALSE,
                             y     = ~fluo, x=~temp,name=position) 
  }
  
  conc_vector <- unique(fluo_m$conc)
  n_colors   <- sapply(conc_vector,get_color_from_conc,min(fluo_m$conc),max(fluo_m$conc))
  
  for (position in unique(fluo_m$conc_)) {
    
    fluo_m_temp <- fluo_m %>% filter(conc_== position) 
    idx        <- which(conc_vector == unique(fluo_m_temp$conc))

    fig <- fig %>% add_trace(data  = fluo_m_temp,
                             color = I(n_colors[idx]),
                             type='scatter',
                             mode = 'lines',
                             line = list(width = line_width),
                             showlegend=FALSE,
                             inherit= FALSE,
                             y = ~fluo, x=~temp,name=position)

  }

  # Now we create a dataframe only to be able to plot a colorbar
  max_val <- log10(max(fluo_m$conc))
  min_val <- log10(min(fluo_m$conc))

  mid_val <- (max_val + min_val) / 2

  # Create a sequence of tick values from max to min
  tickvals <- c(min_val,mid_val,max_val)

  ticktext <- c(
    min(fluo_m$conc),
    10^mid_val,
    max(fluo_m$conc)
  )

  ticktext <- formatC(ticktext, format = "e", digits = 2)

    # Set layout and position the colorbar
  fig <- fig %>% colorbar(
    title = list(text='[Ligand] (M)',font=list(size=legend_text_size)),
    x = xLegend,   # Horizontal position
    y = yLegend,   # Vertical position
    xanchor = "right",  # Anchoring to the right side
    yanchor = "top",
    tickvals = tickvals,  # Ticks from max to min, in logarith scale
    ticktext = ticktext,  # Use the same tick values as labels
    tickfont = list(size = legend_text_size-1),  # Font size of the ticks
    len = color_bar_length,
    orientation = color_bar_orientation,
    outlinewidth = 0)

  fig <- fig %>% layout(xaxis = x, yaxis = y,font="Roboto")
  
  fig <- fig %>%  config(
    toImageButtonOptions = list(
      format = plot_type,
      filename = "myplot",
      width = plot_width * 50,
      height = plot_height * 50
    ), displaylogo = FALSE)
  
  return(fig)
}


# Plot maximum of derivative
generate_der_plot <- function(
  tms,
  concentrations,
  plot_width=12,
  plot_height=12,
  plot_type='png',
  axis_size=14,
  show_x_grid=TRUE,
  show_y_grid=TRUE,
  show_axis_lines=TRUE,
  tickwidth=2,
  ticklen=8,
  marker_size=10,
  log_scale_type='defaultPlotly'){
  
  #log_scale_typecan be defaultPlotly, micromolar, milimolar, molar

  df2 <- generate_der_df(tms,concentrations)
  
  factorList <- list("defaultPlotly"=1,"molar"=1,"milimolar"=1e3,"micromolar"=1e6)
  
  df2$concentration <- df2$concentration * unlist(factorList[log_scale_type])
  
  titleList <- list("defaultPlotly" = "[Ligand] (M)"   , "molar"      = "[Ligand] (M)",
                    "milimolar"     = "[Ligand] (mM)"  , "micromolar" = "[Ligand] (μM)")
  
  xTitle <- unlist(titleList[log_scale_type])
  
  exponentformat <- ifelse(log_scale_type == "defaultPlotly", "B", 'power')
  
  x <- list(
    title = xTitle,
    titlefont = list(size = axis_size),
    tickfont = list(size = axis_size),
    type = "log",
    exponentformat=exponentformat,
    showgrid = show_x_grid,
    showline = show_axis_lines,
    zeroline = FALSE,
    ticks = "outside",
    tickwidth = tickwidth*show_axis_lines,
    ticklen = ticklen*show_axis_lines,
    automargin = TRUE
  )
  
  y <- list(
    title = "Tm (°C) (using the 1st derivative)",
    titlefont = list(size = axis_size),
    tickfont = list(size = axis_size),
    showgrid = show_y_grid,
    showline = show_axis_lines,
    zeroline = FALSE,
    ticks = "outside",
    tickwidth = tickwidth*show_axis_lines,
    ticklen = ticklen*show_axis_lines,
    automargin = TRUE
  )
  
  fig <- plot_ly(
    data = df2,
    y = ~Tm,
    x = ~concentration,
    type = 'scatter',
    mode = "markers",
    marker = list(size = marker_size,
                  color = "rgba(255, 182, 193, .9)")
  )
  
  fig <- fig %>% layout(xaxis = x, yaxis = y,font="Roboto",showlegend = FALSE)
  
  fig <- fig %>%  config(
    toImageButtonOptions = list(
      format = plot_type,
      filename = "myplot",
      width = plot_width * 50,
      height = plot_height * 50
    ), displaylogo = FALSE)
  
  return(fig)
}

## Plot initial fluorescence signal versus ligand concentration
## fluo_m is a 3 column dataframe: temp, fluo and conc. 
## which is the signal (350,330 or Ratio)
## Returns the plot

plot_initialFluo_signal <- function(
  fluo_m,
  which="Ratio",
  plot_width=12,
  plot_height=12,
  plot_type='png',
  axis_size=14,
  show_x_grid=TRUE,
  show_y_grid=TRUE,
  show_axis_lines=TRUE,
  tickwidth=2,
  ticklen=8,
  marker_size=5,
  log_scale_type='defaultPlotly'
) {

  factorList <- list("defaultPlotly"=1,"molar"=1,"milimolar"=1e3,"micromolar"=1e6)

  fluo_m$conc <- fluo_m$conc * unlist(factorList[log_scale_type])

  titleList <- list("defaultPlotly" = "[Ligand] (M)"   , "molar"      = "[Ligand] (M)",
                    "milimolar"     = "[Ligand] (mM)"  , "micromolar" = "[Ligand] (μM)")

  xTitle <- unlist(titleList[log_scale_type])

  exponentformat <- ifelse(log_scale_type == "defaultPlotly", "B", 'power')

  x <- list(
    title = xTitle,
    titlefont = list(size = axis_size),
    tickfont = list(size = axis_size),
    type = "log",
    exponentformat = exponentformat,
    showgrid = show_x_grid,
    showline = show_axis_lines,
    zeroline = FALSE,
    ticks = "outside",
    tickwidth = tickwidth*show_axis_lines,
    ticklen = ticklen*show_axis_lines,
    automargin = TRUE
  )

  y <- list(
    title = which,
    titlefont = list(size = axis_size),
    tickfont = list(size = axis_size),
    showgrid = show_y_grid,
    showline = show_axis_lines,
    zeroline = FALSE,
    ticks = "outside",
    tickwidth = tickwidth*show_axis_lines,
    ticklen = ticklen*show_axis_lines,
    automargin = TRUE
  )
  
  fluo_m  <- fluo_m %>% filter(conc > 0)
  
  fig <- plot_ly(type = 'scatter', mode = 'markers')
  
  if (nrow(fluo_m) == 0) { return(NULL) }
  
  fluo_m <- fluo_m[fluo_m$temp==min(fluo_m$temp),]
  
  # set as black the positions with no ligand
  fig <- fig %>% add_trace(
    data = fluo_m,
    color = I("black"),
    y = ~fluo,
    x = ~conc,
    marker = list(size = marker_size)
  )
  
  fig <- fig %>% layout(xaxis = x, yaxis = y,font="Roboto",showlegend = FALSE)
  
  fig <- fig %>%  config(
    toImageButtonOptions = list(
      format = plot_type,
      filename = "myplot",
      width = plot_width * 50,
      height = plot_height * 50
    ), displaylogo = FALSE)
  
  return(fig)
}

## Plot the fluorescence fit.
## Requires the list of dataframes from both the experimental and fitted fluorescence
## Plots the selected element from that list

plot_fluorescence_fit <- function(fluo_m_real_list,fluo_m_pred_list,selected) {
  
  fluo_m <- fluo_m_real_list
  
  fluo_m <- fluo_m[[selected]]
  
  if (max(fluo_m$conc) > 1e-3) {
    fluo_m$conc <- fluo_m$conc*1e3
    scale_f <- "mM"
  } else {
    fluo_m$conc <- fluo_m$conc*1e6
    scale_f <- "µM"
  }
  
  fluo_m$legend <- paste0("Rep ",fluo_m$rep, ". ",signif(fluo_m$conc,2)," ",scale_f)
  neworder <- unique(fluo_m$legend)
  fluo_m <- arrange(transform(fluo_m,legend=factor(legend,levels=neworder)),legend)
  
  p <- ggplot(fluo_m,aes(x=temp,y=fluo))+
    geom_point(size=0.4)+
    theme_bw(base_size = 14)+
    xlab("Temperature (ºC)")+
    ylab("Fluorescence (AU)")+
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"))
  
  fluo_m <- fluo_m_pred_list[[selected]]
  
  if (max(fluo_m$conc) > 1e-3) {
    fluo_m$conc <- fluo_m$conc*1e3
    scale_f <- "mM"
  } else {
    fluo_m$conc <- fluo_m$conc*1e6
    scale_f <- "µM"
  }
  
  fluo_m$legend <- paste0("Rep ",fluo_m$rep, ". ",signif(fluo_m$conc,2)," ",scale_f)
  neworder <- unique(fluo_m$legend)
  fluo_m <- arrange(transform(fluo_m,legend=factor(legend,levels=neworder)),legend)
  
  p2 <- p +
    geom_line(data=fluo_m,color="red",linewidth=0.3,
              aes(x=temp,y=fluo),inherit.aes = FALSE)
  
  p3 <- p2 + facet_wrap(~ legend,scales = "free",
                        ncol = global_plot_columns,nrow = global_plot_rows,drop = FALSE)
  
  return(p3)
}

## Plot selected fluorescence fit if possible

plot_sel_fluo_fit <- function(real_data,model_data,selected){
  if (length(real_data) < selected) {return(NULL)}
  
  p <- plot_fluorescence_fit(real_data,model_data,selected)
  return(p)
}

## Plot the isothermal exp data
## iso_real_data has two columns: conc and fu (unfolded_fraction, "experimental")
## Returns the plot

plot_isothermal_fitting_exp <- function(
  iso_real_data,
  plot_width=12,
  plot_height=12,
  plot_type='png',
  axis_size=12,
  legend_size=12,
  show_x_grid=TRUE,
  show_y_grid=TRUE,
  show_axis_lines=TRUE,
  tickwidth=2,
  ticklen=8,
  marker_size=2,
  log_scale_type = 'defaultPlotly') {

  factorList <- list("defaultPlotly"=1,"molar"=1,"milimolar"=1e3,"micromolar"=1e6)
  iso_real_data$conc <- iso_real_data$conc * unlist(factorList[log_scale_type])
  titleList <- list("defaultPlotly" = "[Ligand] (M)"   , "molar"      = "[Ligand] (M)",
                    "milimolar"     = "[Ligand] (mM)"  , "micromolar" = "[Ligand] (μM)")
  xTitle <- unlist(titleList[log_scale_type])
  exponentformat <- ifelse(log_scale_type == "defaultPlotly", "B", 'power')

  xaxis <- list(title = xTitle,titlefont = list(size = axis_size),
            tickfont = list(size = axis_size),type = "log",exponentformat=exponentformat,
            showgrid = show_x_grid,
            showline = show_axis_lines,
            zeroline = FALSE,
            ticks = "outside",
            tickwidth = tickwidth*show_axis_lines,
            ticklen = ticklen*show_axis_lines,
            automargin = TRUE)
  yaxis <- list(title = "Fraction unfolded",titlefont = list(size = axis_size),
            tickfont = list(size = axis_size),
            showgrid = show_y_grid,
            showline = show_axis_lines,
            zeroline = FALSE,
            ticks = "outside",
            tickwidth = tickwidth*show_axis_lines,
            ticklen = ticklen*show_axis_lines,
            automargin = TRUE)
  plot_list <- list()
  x_pos_annot <- min(log10(iso_real_data$conc)) + (max(log10(iso_real_data$conc)) - min(log10(iso_real_data$conc))) / 2
  tot_cond <- length(unique(iso_real_data$legend))
  i <- 0
  for (leg in unique(iso_real_data$legend)) {

    i <- i + 1
    df_temp <- iso_real_data %>% filter(legend == leg)
    fig <- plot_ly(df_temp, x = ~conc, y = ~fu, color = I("#00AFBB"), type = "scatter",
                   marker = list(size = marker_size)) %>%
      layout(xaxis = xaxis,yaxis=yaxis,
             annotations = list(x = x_pos_annot , y = 1.12, 
                                text = ifelse(tot_cond > 1,leg,""), showarrow = F, 
                                xref='x', yref='paper',font = list(size = legend_size)),
             showlegend=FALSE)
    plot_list[[i]] <- fig
  }
  if (i > 1) {
    return( subplot(plot_list,nrows = 2,margin = c(0.03,0.03,0.1,0.1),
                    shareY = TRUE, shareX = TRUE) %>%  config(
      toImageButtonOptions = list(
        format = plot_type,
        filename = "myplot",
        width = plot_width * 50,
        height = plot_height * 50
      ), displaylogo = FALSE,
      modeBarButtonsToRemove = list('sendDataToCloud',
                                    'hoverClosestCartesian','hoverCompareCartesian',
                                    'lasso2d','select2d','zoomIn2d','zoomOut2d',
                                    'zoom2d','pan2d','autoScale2d','resetScale2d')))
  } else {
    return( plot_list[[1]] %>%  layout(
      title=list(text=unique(iso_real_data$legend),font = list(size = legend_size))) %>% config(
      toImageButtonOptions = list(
        format = plot_type,
        filename = "myplot",
        width = plot_width * 50,
        height = plot_height * 50
      ), displaylogo = FALSE,
      modeBarButtonsToRemove = list('sendDataToCloud',
                                    'hoverClosestCartesian','hoverCompareCartesian',
                                    'lasso2d','select2d','zoomIn2d','zoomOut2d',
                                    'zoom2d','pan2d','autoScale2d','resetScale2d')))
  }
}

## Plot the isothermal fitting
## iso_real_data has two columns: conc and fu (unfolded_fraction, "experimental")
## iso_model_data has four columns: conc fu (predicted) fu_low fu_upp (lower and upper bound of the error)
## Returns the plot

plot_isothermal_fitting <- function(
  iso_real_data,
  iso_model_data,
  plot_width=12,
  plot_height=12,
  plot_type='png',
  axis_size=12,
  legend_size=12,
  show_x_grid=T,
  show_y_grid=T,
  show_axis_lines=T,
  tickwidth=2,
  ticklen=8,
  line_width=2,
  marker_size=8,
  log_scale_type= 'defaultPlotly') {
  
  factorList <- list("defaultPlotly"=1,"molar"=1,"milimolar"=1e3,"micromolar"=1e6)
  
  iso_real_data$conc  <- iso_real_data$conc * unlist(factorList[log_scale_type])
  iso_model_data$conc <- iso_model_data$conc * unlist(factorList[log_scale_type])
  
  titleList <- list("defaultPlotly" = "[Ligand] (M)"   , "molar"      = "[Ligand] (M)",
                    "milimolar"     = "[Ligand] (mM)"  , "micromolar" = "[Ligand] (μM)")
  
  xTitle <- unlist(titleList[log_scale_type])

  exponentformat <- ifelse(log_scale_type == "defaultPlotly", "B", 'power')
  
  xaxis <- list(
    title = xTitle,
    titlefont = list(size = axis_size),
    tickfont = list(size = axis_size),
    type = "log",
    exponentformat=exponentformat,
    showgrid = show_x_grid,
    showline = show_axis_lines,
    zeroline = FALSE,
    ticks = "outside",
    tickwidth = tickwidth*show_axis_lines,
    ticklen = ticklen*show_axis_lines,
    automargin = TRUE
  )
  
  yaxis <- list(
    title = "Fraction unfolded",
    titlefont = list(size = axis_size),
    tickfont = list(size = axis_size),
    showgrid = show_y_grid,
    showline = show_axis_lines,
    zeroline = FALSE,
    ticks = "outside",
    tickwidth = tickwidth*show_axis_lines,
    ticklen = ticklen*show_axis_lines,
    automargin = TRUE
  )
  
  plot_list <- list()
  
  x_pos_annot <- min(log10(iso_real_data$conc)) + (max(log10(iso_real_data$conc)) - min(log10(iso_real_data$conc))) / 2
  
  tot_cond <- length(unique(iso_real_data$legend))
  
  i <- 0
  for (leg in unique(iso_model_data$legend)) {
    i <- i + 1
    
    df_temp <- iso_model_data %>% filter(legend == leg) %>% arrange(conc)
    
    fig <- plot_ly(type = 'scatter',mode = 'lines+markers')
    
    fig <- fig %>% add_trace(data=df_temp, x = ~conc,y = ~fu, type = 'scatter', mode = 'lines',
                             line = list(color = 'rgba(255,163,0,0.5)', width = line_width))

    df_temp2 <- iso_real_data %>% filter(legend == leg) %>% arrange(conc)
    
    fig <- fig %>% add_trace(data=df_temp2, x = ~conc, y = ~fu, color = I("#00AFBB"),
                             type = "scatter",mode = "markers",
                             marker = list(size = marker_size))

    # Replace Kd with K in italic and d in subscript
    leg <- gsub("Kd","<i>K</i><sub>d</sub>",leg)

    # Convert CI95 into CI<sub>95</sub>
    leg <- gsub("CI95","CI<sub>95</sub>",leg)

    fig <- fig %>% layout(xaxis = xaxis,yaxis=yaxis,
                          annotations = list(x = x_pos_annot , y = 1.13,
                                             text = ifelse(tot_cond > 1,leg,""), showarrow = F,
                                             xref='x', yref='paper',font = list(size = legend_size)),
                          showlegend=FALSE)
    
    plot_list[[i]] <- fig
    
  }
  
  if (i > 1) {
    return( subplot(plot_list,nrows = 2,margin = c(0.03,0.03,0.1,0.1),
                    shareY = TRUE, shareX = TRUE) %>%  config(
      toImageButtonOptions = list(
        format = plot_type,
        filename = "myplot",
        width = plot_width * 50,
        height = plot_height * 50
      ), displaylogo = FALSE,
      modeBarButtonsToRemove = list('sendDataToCloud',
                                    'hoverClosestCartesian','hoverCompareCartesian',
                                    'lasso2d','select2d','zoomIn2d','zoomOut2d',
                                    'zoom2d','pan2d','autoScale2d','resetScale2d')))
  } else {
    m   <- list(l = 10,r = 10,b = 20,t = 40,pad = 4)
    return( plot_list[[1]] %>%  layout(title=list(text=unique(iso_real_data$legend),
                                                  font = list(size = legend_size)),
                                       margin=m) %>%  config(
      toImageButtonOptions = list(
        format = plot_type,
        filename = "myplot",
        width = plot_width * 50,
        height = plot_height * 50
      ), displaylogo = FALSE,
      modeBarButtonsToRemove = list('sendDataToCloud',
                                    'hoverClosestCartesian','hoverCompareCartesian',
                                    'lasso2d','select2d','zoomIn2d','zoomOut2d',
                                    'zoom2d','pan2d','autoScale2d','resetScale2d')))
  }
  
}

plot_tm_shift <- function(
  conc,
  tms,
  plot_width=12,
  plot_height=12,
  plot_type='png',
  axis_size=14,
  show_x_grid=TRUE,
  show_y_grid=TRUE,
  show_axis_lines=TRUE,
  tickwidth=2,
  ticklen=8,
  marker_size=10,
  log_scale_type='defaultPlotly') {
  
  df_model <- get_tm_df(conc,tms)
  
  df_model$tms <- df_model$tms - 273.15 # To centigrade
  
  factorList <- list("defaultPlotly"=1,"molar"=1,"milimolar"=1e3,"micromolar"=1e6)
  df_model$l_conc <- df_model$l_conc * unlist(factorList[log_scale_type])
  
  titleList <- list("defaultPlotly" = "[Ligand] (M)"   , "molar"      = "[Ligand] (M)",
                    "milimolar"     = "[Ligand] (mM)"  , "micromolar" = "[Ligand] (μM)")
  
  xTitle <- unlist(titleList[log_scale_type])

  exponentformat <- ifelse(log_scale_type == "defaultPlotly", "B", 'power')

  xaxis <- list(
    title = xTitle,
    titlefont = list(size = axis_size),
    tickfont = list(size = axis_size),
    type = "log",
    exponentformat=exponentformat,
    showgrid = show_x_grid,
    showline = show_axis_lines,
    zeroline = FALSE,
    ticks = "outside",
    tickwidth = tickwidth*show_axis_lines,
    ticklen = ticklen*show_axis_lines,
    automargin = TRUE
  )
  
  yaxis <- list(
    title = "Observed T<sub>m</sub> (°C)",
    titlefont = list(size = axis_size),
    tickfont = list(size = axis_size),
    showgrid = show_y_grid,
    showline = show_axis_lines,
    zeroline = FALSE,
    ticks = "outside",
    tickwidth = tickwidth*show_axis_lines,
    ticklen = ticklen*show_axis_lines,
    automargin = TRUE
  )
  
  fig <- plot_ly(
    data=df_model,
    x = ~l_conc,
    y = ~tms,
    color = I("#00AFBB"),
    type = "scatter",
    inherit = FALSE,
    mode = "markers",
    marker = list(size = marker_size)
    ) %>%
    layout(xaxis = xaxis,yaxis=yaxis,showlegend=FALSE) %>%  
    config(toImageButtonOptions = list(
      format = plot_type,
      filename = "myplot",
      width = plot_width * 50,
      height = plot_height * 50
    ), displaylogo = FALSE,
    modeBarButtonsToRemove = list('sendDataToCloud',
                                  'hoverClosestCartesian','hoverCompareCartesian',
                                  'lasso2d','select2d','zoomIn2d','zoomOut2d',
                                  'zoom2d','pan2d'))
  
  return(fig)
}

plot_tm_shift_fit <- function(
  conc,
  tms,
  df_pred,
  fit_info,
  asymmetricCI95,
  plot_width=12,
  plot_height=12,
  plot_type='png',
  axis_size=14,
  legend_size=13,
  show_x_grid=TRUE,
  show_y_grid=TRUE,
  show_axis_lines=TRUE,
  tickwidth=2,
  ticklen=8,
  marker_size=10,
  line_width=2,
  log_scale_type='defaultPlotly') {
  
  legend <- paste0("<i>K</i><sub>d,1</sub> (M) = ",
                   signif(fit_info$estimate[2],2)," ± ",
                   round(fit_info$std.error[2]/fit_info$estimate[2]*100,0),"%"
  )
  
  if (!(is.null(asymmetricCI95))) {
    
    legend <- paste0("<i>K</i><sub>d</sub> (M) : ",
                     signif(fit_info$estimate[2],2)," Asymmetric CI<sub>95</sub> : [",
                     signif(asymmetricCI95$kd_min95,2)," ; ",
                     signif(asymmetricCI95$kd_max95,2),"]"
                     
    )
    
  }
  
  if (nrow(fit_info) == 3){
    legend <- paste0(legend,
                     "<i>K</i><sub>d,2</sub> (M) = ",
                     signif(fit_info$estimate[3],2)," ± ",
                     round(fit_info$std.error[3]/fit_info$estimate[3]*100,0),"%"
    )
  } 
  
  df_pred$tms <- df_pred$tms - 273.15 # To centigrade
  
  df_pred <- df_pred %>% arrange(l_conc)
  
  factorList     <- list("defaultPlotly"=1,"molar"=1,"milimolar"=1e3,"micromolar"=1e6)
  df_pred$l_conc <- df_pred$l_conc * unlist(factorList[log_scale_type])
  
  m   <- list(l = 10,r = 10,b = 20,t = 40,pad = 4)
  
  fig <- plot_tm_shift(
    conc,
    tms,
    plot_width=plot_width,
    plot_height=plot_height,
    plot_type=plot_type,
    axis_size=axis_size,
    show_x_grid=show_x_grid,
    show_y_grid=show_y_grid,
    show_axis_lines=show_axis_lines,
    tickwidth=tickwidth,
    ticklen=ticklen,
    marker_size=marker_size,
    log_scale_type=log_scale_type
  ) %>%
    add_trace(
      data=df_pred,
      x = ~l_conc,
      y = ~tms,
      type = 'scatter',
      mode = 'lines',
      inherit = FALSE,
      showlegend=FALSE,
      line = list(
        width=line_width,
        color = 'rgba(255,163,0,0.5)')) %>%
    layout(title=list(text=legend, font=list(size=legend_size)), margin=m)

  return(fig)
}

box(title = "Tm shift params and data", width = 4, solidHeader = T, status = "primary", 
    fluidRow(
      column(8, p(style = "font-size: 120%",HTML(""),
      downloadLink('download_tm_fit_plot', 'Observed vs Fitted Tm shifts')))
    ))
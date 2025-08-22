box(
  title = "Plot options", width = 12, solidHeader = T, status = "primary",

  conditionalPanel(
    condition = "input.show_plot_download_options_input_tmShift",

    fluidRow(
      column(6,
        p(HTML("<b>Width</b>"),
          span(shiny::icon("info-circle"), id = "info_uu1-13_tmShift"),
          numericInput('plot_width_tmShift', NULL, 18, min = 1, max = 100),
          tippy::tippy_this(
            elementId = "info_uu1-13_tmShift",
            tooltip = "Units are pixels * 50", placement = "right"
          )
        )
      ),

      column(6,
        p(HTML("<b>Height</b>"),
          span(shiny::icon("info-circle"), id = "info_uu1-14_tmShift"),
          numericInput('plot_height_tmShift', NULL, 11, min = 1, max = 100),
          tippy::tippy_this(
            elementId = "info_uu1-14_tmShift",
            tooltip = "Units are pixels * 50", placement = "right"
          )
        )
      ),

      column(6,
        p(HTML("<b>File type</b>"),
          selectInput("plot_type_tmShift", NULL,
            c("PNG" = "png", "SVG" = "svg", "JPEG" = "jpeg")
          )
        )
      )
    ),

    fluidRow(
      column(6,
        p(HTML("<b>Axis font</b>"),
          sliderInput('plot_axis_size_tmShift', NULL, 14, min = 4, max = 40, step = 1)
        )
      ),

      column(6,
        p(HTML("<b>Legend font</b>"),
          sliderInput('plot_legend_size_tmShift', NULL, 13, min = 6, max = 40, step = 1)
        )
      )
    ),

    fluidRow(
      column(6,
        p(HTML("<b>Show x-grid</b>"),
          checkboxInput("show_x_grid_tmShift", NULL, TRUE)
        )
      ),

      column(6,
        p(HTML("<b>Show y-grid</b>"),
          checkboxInput("show_y_grid_tmShift", NULL, TRUE)
        )
      )
    ),

    fluidRow(
      column(6,
        p(HTML("<b>Marker size </b>"),
          sliderInput('plot_marker_size_tmShift', NULL, 8, min = 1, max = 16, step = 0.5)
        )
      ),

      column(6,
        p(HTML("<b>Line width</b>"),
          sliderInput('plot_line_width_tmShift', NULL, 2, min = 1, max = 10, step = 0.5)
        )
      )
    ),

    fluidRow(
      column(6,
        p(HTML("<b>Tick width</b>"),
          sliderInput('plot_tickwidth_tmShift', NULL, 1, min = 1, max = 10, step = 0.5)
        )
      ),

      column(6,
        p(HTML("<b>Tick length</b>"),
          sliderInput('plot_ticklen_tmShift', NULL, 8, min = 1, max = 30, step = 0.5)
        )
      )
    ),

    fluidRow(
      column(6,
        p(HTML("<b>Show axis lines</b>"),
          checkboxInput("show_axis_lines_tmShift", NULL, TRUE)
        )
      ),

      column(6,
        p(HTML("<b>Log Axis type</b>"),
          selectInput("logScaleType_tmShift", NULL,
            c(
              "μM" = "micromolar",
              "mM" = "milimolar",
              "M" = "molar",
              "Plotly" = "defaultPlotly"
            )
          )
        )
      )
    )
  ),

  fluidRow(
    column(6,
      p(HTML("<b>Show/Hide</b>"),
        checkboxInput("show_plot_download_options_input_tmShift", "", FALSE)
      )
    )
  )
)

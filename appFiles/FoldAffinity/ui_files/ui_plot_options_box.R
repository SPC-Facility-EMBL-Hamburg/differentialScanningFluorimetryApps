box(
  title = "Plot options", width = 12, solidHeader = T, status = "primary",

  conditionalPanel(
    condition = "input.show_plot_download_options_input",

    fluidRow(
      column(6,
        p(HTML("<b>Width</b>"),
          span(shiny::icon("info-circle"), id = "info_uu1-13"),
          numericInput('plot_width', NULL, 18, min = 1, max = 100),
          tippy::tippy_this(
            elementId = "info_uu1-13",
            tooltip = "Units are pixels * 50", placement = "right"
          )
        )
      ),

      column(6,
        p(HTML("<b>Height</b>"),
          span(shiny::icon("info-circle"), id = "info_uu1-14"),
          numericInput('plot_height', NULL, 11, min = 1, max = 100),
          tippy::tippy_this(
            elementId = "info_uu1-14",
            tooltip = "Units are pixels * 50", placement = "right"
          )
        )
      ),

      column(6,
        p(HTML("<b>File type</b>"),
          selectInput("plot_type", NULL,
            c("PNG" = "png", "SVG" = "svg", "JPEG" = "jpeg")
          )
        )
      )
    ),

    fluidRow(
      column(6,
        p(HTML("<b>Axis font</b>"),
          sliderInput('plot_axis_size', NULL, 14, min = 4, max = 40, step = 1)
        )
      ),

      column(6,
        p(HTML("<b>Legend font</b>"),
          sliderInput('plot_legend_size', NULL, 13, min = 6, max = 40, step = 1)
        )
      )
    ),

    fluidRow(
      column(6,
        p(HTML("<b>Show x-grid</b>"),
          checkboxInput("show_x_grid", NULL, TRUE)
        )
      ),

      column(6,
        p(HTML("<b>Show y-grid</b>"),
          checkboxInput("show_y_grid", NULL, TRUE)
        )
      )
    ),

    fluidRow(
      column(6,
        p(HTML("<b>Marker size </b>"),
          sliderInput('plot_marker_size', NULL, 8, min = 1, max = 16, step = 0.5)
        )
      ),

      column(6,
        p(HTML("<b>Line width</b>"),
          sliderInput('plot_line_width', NULL, 2, min = 1, max = 10, step = 0.5)
        )
      )
    ),

    fluidRow(
      column(6,
        p(HTML("<b>Colorbar orientation</b>"),
          selectInput("color_bar_orientation", NULL,
            c("vertical" = "v", "horizontal" = "h")
          )
        )
      ),

      column(6,
        p(HTML("<b>Color bar length</b>"),
          sliderInput('color_bar_length', NULL, 0.5, min = 0.2, max = 1, step = 0.1)
        )
      )
    ),

    fluidRow(
      column(6,
        p(HTML("<b>Colorbar x-axis pos.</b>"),
          sliderInput('x_legend_pos', NULL, 1.1, min = 0, max = 1.3, step = 0.1)
        )
      ),

      column(6,
        p(HTML("<b>Colorbar y-axis pos.</b>"),
          sliderInput('y_legend_pos', NULL, 1.1, min = 0, max = 1.3, step = 0.1)
        )
      )
    ),

    fluidRow(
      column(6,
        p(HTML("<b>Tick width</b>"),
          sliderInput('plot_tickwidth', NULL, 1, min = 1, max = 10, step = 0.5)
        )
      ),

      column(6,
        p(HTML("<b>Tick length</b>"),
          sliderInput('plot_ticklen', NULL, 8, min = 1, max = 30, step = 0.5)
        )
      )
    ),

    fluidRow(
      column(6,
        p(HTML("<b>Show axis lines</b>"),
          checkboxInput("show_axis_lines", NULL, TRUE)
        )
      ),

      column(6,
        p(HTML("<b>Log Axis type</b>"),
          selectInput("logScaleType", NULL,
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
        checkboxInput("show_plot_download_options_input", "", FALSE)
      )
    )
  )
)

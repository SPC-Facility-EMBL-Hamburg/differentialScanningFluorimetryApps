box(title = "Plot download options", width = 12, solidHeader = T, status = "primary",

    conditionalPanel(condition = "input.show_plot_download_options_input",

        fluidRow(

            column(6, p(HTML("<b>Width</b>"),
            span(shiny::icon("info-circle"), id = "info_uu1-13"),
            numericInput('plot_width',NULL, 18,min = 1, max = 100),
            tippy::tippy_this(elementId = "info_uu1-13",
            tooltip = "Units are pixels * 50",placement = "right"))),

            column(6, p(HTML("<b>Height</b>"),
            span(shiny::icon("info-circle"), id = "info_uu1-14"),
            numericInput('plot_height',NULL, 11,min = 1, max = 100),
            tippy::tippy_this(elementId = "info_uu1-14",
            tooltip = "Units are pixels * 50",placement = "right")))

        ),

        fluidRow(

            column(6, p(HTML("<b>Legend text size</b>"),
            numericInput('plot_font_size',NULL, 13,min = 6, max = 40))),

            column(6, p(HTML("<b>Axis text size</b>"),
            numericInput('plot_axis_size',NULL, 14,min = 4, max = 40)))

        ),

        fluidRow(

            column(6, p(HTML("<b>Show x-grid</b>"),
            checkboxInput("show_x_grid", "", TRUE))),

            column(6, p(HTML("<b>Show y-grid</b>"),
            checkboxInput("show_y_grid", "", TRUE)))

        ),

        fluidRow(

            column(6, p(HTML("<b>Show axis lines</b>"),
            checkboxInput("show_axis_lines", "", TRUE))),

            column(6, p(HTML("<b>File type</b>"),
            selectInput("plot_type", NULL,
            c("PNG" = "png", "SVG" = "svg","JPEG" = "jpeg"))))

        ),

        fluidRow(

            column(6, p(HTML("<b>Tick width</b>"),
            numericInput('plot_tickwidth', NULL, 1, min = 1, max = 10))),

            column(6, p(HTML("<b>Tick length</b>"),
            numericInput('plot_ticklen', NULL, 8, min = 1, max = 30)))

        ),

         fluidRow(

            column(6, p(HTML("<b>Line width</b>"),
            numericInput('plot_line_width', NULL, 2, min = 0, max = 10,step = 0.5))),

        ),

    ),

    fluidRow(

        column(6, p(HTML("<b>Show options</b>"),
        checkboxInput("show_plot_download_options_input", "", FALSE))),

        column(6, p(HTML("<b>Show colors</b>"),
        checkboxInput("show_colors_column", "", TRUE)))

    )

)

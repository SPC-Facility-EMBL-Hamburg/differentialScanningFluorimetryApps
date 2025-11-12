box(
    title = "Filters", width = 5, solidHeader = T, status = "primary",

    fluidRow(

        column(8,

            p(HTML("<br>"),

                checkboxInput(
                    "sd_factor_bool",
                    "All fitted params: (std / value) * 100 < Threshold", FALSE
                )

            )

        ),

        column(4,

            p(HTML("Threshold"),

                numericInput(
                    "rel_error_threshold",
                    NULL,
                    50,
                    min = 1,
                    max = 1e3,
                    step = 10
                )

            )

        )

    ),

    fluidRow(

        column(4,

            p(HTML("<b>Residuals threshold</b>"),

                numericInput(
                    "fitting_std_threshold",
                    NULL,
                    1e4,
                    min = 0,
                    max = 1e6,
                    step = 10
                )

            )

        )

    ),

    fluidRow(

        conditionalPanel(
            "output.model_name == 'EquilibriumTwoState' ||
            output.model_name == 'EmpiricalTwoState'",

            column(3,

                p(HTML("<b>Lower Tm (°C)</b>"),

                    numericInput(
                        "lower_Tm_threshold",
                        NULL,
                        25,
                        min = 5,
                        max = 100,
                        step = 1
                    )

                )

            ),

            column(3,

                p(HTML("<b>Upper Tm (°C)</b>"),

                    numericInput(
                        "upper_Tm_threshold",
                        NULL,
                        85,
                        min = 5,
                        max = 100,
                        step = 1
                    )

                )

            )

        )

    ),

    fluidRow(

        conditionalPanel(
            "output.model_name == 'EquilibriumTwoState'",

            column(4,

                p(HTML("<b>Lower DH (kcal/mol)</b>"),

                    numericInput(
                        "lower_dh_threshold",
                        NULL,
                        20,
                        min = 5,
                        max = 300,
                        step = 5
                    )

                )

            ),

            column(4,

                p(HTML("<b>Upper DH (kcal/mol)</b>"),

                    numericInput(
                        "upper_dh_threshold",
                        NULL,
                        450,
                        min = 100,
                        max = 600,
                        step = 5
                    )

                )

            )

        )
    )

)
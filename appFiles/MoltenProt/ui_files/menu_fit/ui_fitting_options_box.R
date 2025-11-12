box(title = "2.2 Fitting options", width = 3, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(8, p(HTML("<b>Temperature range for baseline estimation</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu2-2"),
                  numericInput('temp_range_baseline_estimation',NULL, 8,min = 0, max = 20),
                  tippy::tippy_this(
                    elementId = "info_uu2-2",
                    tooltip = "Used to obtain initial guesses for the fluorescence dependence
                     on temperature dependence of
                    the native and unfolded states .",
                    placement = "right")))
    ),

    fluidRow(

        column(6,p(
            HTML("<b>Signal dependence (native)</b>"),
            span(shiny::icon("info-circle"), id = "info_uu-native_dependence"),
            selectInput("native_dependence", NULL,choices =
                c(
                "quadratic"= "quadratic",
                "linear"   = "linear",
                "constant" = "constant"
                )),
            tippy::tippy_this(
            elementId = "info_uu-native_dependence",
            tooltip = "Set to linear if there's a linear dependence
            between the signal and the temperature for the
            native state.
            Set to quadratic if there's a quadratic dependence.",placement = "right"))
        ),

        column(6,p(
            HTML("<b>Signal dependence (unfolded)</b>"),
            span(shiny::icon("info-circle"), id = "info_uu-unfolded_dependence"),
            selectInput("unfolded_dependence", NULL,choices =
                c(
                "quadratic"= "quadratic",
                "linear"   = "linear",
                "constant" = "constant"
                )),
            tippy::tippy_this(
            elementId = "info_uu-unfolded_dependence",
            tooltip = "Set to linear if there's a linear dependence
            between the signal and the temperature for the
            unfolded state.
            Set to quadratic if there's a quadratic dependence.",placement = "right"))
        ),

      conditionalPanel(condition = "output.three_state_model_selected",
                       
         column(6, p(
           HTML("<b>T1 min</b>"),
           span(shiny::icon("info-circle"), id = "info_uu2-t1min"),
           numericInput('t1min',NULL, 0,min = -30, max = 30),
           tippy::tippy_this(elementId = "info_uu2-t1min",
           tooltip = "Left bound for the fitting of the parameter T1.
           Temperature at which [N] = [I] (N and I are respectively
           the native and unfolded intermediate states).",placement = "right"))),
         
         column(6, p(
           HTML("<b>T1 max</b>"),
           span(shiny::icon("info-circle"), id = "info_uu2-t1max"),
           numericInput('t1max',NULL, 0,min = -30, max = 30),
           tippy::tippy_this(elementId = "info_uu2-t1max",
           tooltip = "Right bound for the fitting of the parameter T1.
           Temperature at which [N] = [I] (N and I are respectively
           the native and unfolded intermediate states).",placement = "right"))),
         
         column(6, p(
           HTML("<b>T2 min</b>"),
           span(shiny::icon("info-circle"), id = "info_uu2-t2min"),
           numericInput('t2min',NULL, 0,min = -30, max = 30),
           tippy::tippy_this(elementId = "info_uu2-t2min",
           tooltip = "Left bound for the fitting of the parameter T2.
           Temperature at which [I] = [U] (U and I are respectively
           the unfolded and unfolded intermediate states).",placement = "right"))),
         
         column(6, p(
           HTML("<b>T2 max</b>"),
           span(shiny::icon("info-circle"), id = "info_uu2-t2max"),
           numericInput('t2max',NULL, 0,min = -30, max = 30),
           tippy::tippy_this(elementId = "info_uu2-t2max",
           tooltip = "Right bound for the fitting of the parameter T2.
           Temperature at which [I] = [U] (U and I are respectively
           the unfolded and unfolded intermediate states).",placement = "right"))),
         
                       
      ),
      
      conditionalPanel(condition = "input.model_selected == 'EquilibriumTwoState'",
                       column(4, p(HTML("<b>ΔCP (kcal/°C/mol)</b>"),
                                   span(shiny::icon("info-circle"), id = "info_uu2-3"),
                                   numericInput('delta_cp',NULL, 0,min = -30, max = 30),
                                   tippy::tippy_this(elementId = "info_uu2-3",
                                                     tooltip = "Heat capacity. This value is used to estimate the 'Cp contribution' in the Equilibrium Two State model
                                                     and calculate corrected free energy terms. 
                                                     Without prior knowledge, this value should be left at 0.",placement = "right")))),
      
      conditionalPanel(condition = "input.model_selected == 'IrreversibleTwoState'",
                       column(4, p(HTML("<b>Scan Rate (degrees / minute)</b>"),
                                   span(shiny::icon("info-circle"), id = "info_uu2-32"),
                                   numericInput('scan_rate',NULL, 1,min = 0.2, max = 6),
                                   tippy::tippy_this(elementId = "info_uu2-32",
                                                     tooltip = "Used in ...",placement = "right"))))             
      
    ))
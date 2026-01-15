box(title = "2.1 Model selection", width = 3, solidHeader = T, status = "primary", 
    fluidRow(
      column(10, p(HTML("<b></b>"),
                   span(shiny::icon("info-circle"), id = "info_uu2-1"),
                   selectInput(
                     "model_selected", NULL,
                     c("Equilibrium Two-state v1.1"      = "EquilibriumTwoState",
                       "Empirical Two-state v1.1"        = "EmpiricalTwoState",
                       "Equilibrium Three-state v1.1"    = "EquilibriumThreeState",
                       "Empirical Three-state v1.1"      = "EmpiricalThreeState",
                       "Irreversible Two-state v1.1"     = "IrreversibleTwoState"
                     ),
                     selectize = FALSE),
                   tippy::tippy_this(elementId = "info_uu2-1",
                                     tooltip = "More information about the models can be found in the User Guide section.
                                     The result of the fitting procedure is presented in the 'Fitted conditions' Table. 
                                     If the model fails to fit the data, we suggest first changing the temperature window,
                                     and second, changing the 'Temperature range for baseline estimation'.",placement = "right"))),
      
      column(6,   withBusyIndicatorUI(
        actionButton("btn_cal","Run Fitting",class = "btn-primary")))
    ))
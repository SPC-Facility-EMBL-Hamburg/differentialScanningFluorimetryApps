full_spectrum_dialog <- function() {

  showModal(
    modalDialog(
      fluidRow(
        column(6,
          p(HTML('<p style="margin-bottom:0px;"><br></p>'),
            span(shiny::icon("info-circle"), id = "info_uu-full_spectrum_update_plots"),
            actionButton(
              inputId = "updateSpectralPlots",label = "Create whole spectra plots",
              icon("square", lib = "font-awesome", class = "fa-solid"),
              style="color: #fff; background-color: #337ab7;border-color: #2e6da4"),
            tippy::tippy_this(
              elementId = "info_uu-full_spectrum_update_plots",
              tooltip = "Generate new plots where, for each condition, the spectral data at all
              temperatures is shown.",
              placement = "right")
          )
        ),

        column(6,
          p(HTML("<b>Wavelength  range (nm)</b>"),
            span(shiny::icon("info-circle"), id = "info_uu-full_spectrum_wl_range"),
            sliderInput(
              "wl_range", NULL,
              min = dsf$min_wavelength,
              max = dsf$max_wavelength,
              value = c(reactives$min_wl,reactives$max_wl)
            ),
            tippy::tippy_this(
              elementId = "info_uu-full_spectrum_wl_range",
              tooltip = "Select the wavelength range of interest to show in the whole spectra
              plots and to use for the SVD decomposition.",
              placement = "right")
          )
        )
      ),

      fluidRow(
        column(6,
          p(HTML('<p style="margin-bottom:0px;"><br></p>'),
            span(shiny::icon("info-circle"), id = "info_uu-full_spectrum_svd"),
            actionButton(
              inputId = "decomposeSpectra",
              label = "SVD decomposition",
              icon("square", lib = "font-awesome", class = "fa-solid"),
              style="color: #fff; background-color: #337ab7;border-color: #2e6da4"),
            tippy::tippy_this(
              elementId = "info_uu-full_spectrum_svd",
              tooltip = "Apply SVD decomposition and extract the first SVD coefficient for each condition.",
              placement = "right")
          )
        ),

        # Little hack to use the withBusyIndicatorUI function (loading spinner)
        # This button is not visible in the UI
        column(1,
          p(HTML("<b><br></b>"),
            withBusyIndicatorUI(
              shinyjs::hidden(
                actionButton(
                  "spectraDecompositionHiddenButton",NULL,
                  class = "btn-primary"
                )
              )
            )
          )
        )
        # End of Little hack to use the withBusyIndicatorUI function (loading spinner)
      ),

      fluidRow(
        column(6,
          p(HTML('<p style="margin-bottom:0px;"><br></p>'),
            span(shiny::icon("info-circle"), id = "info_uu-full_spectrum_ratio"),
            actionButton(
              inputId = "generateRatioSignal",label = "Generate ratio signal",
              icon("square", lib = "font-awesome", class = "fa-solid"),
              style="color: #fff; background-color: #337ab7;border-color: #2e6da4"),
            tippy::tippy_this(
              elementId = "info_uu-full_spectrum_ratio",
              tooltip = "Generate a new signal that will be the ratio between the signal at two wavelengths:
                Signal_ratio = Signal(WL1) / Signal(WL2).
              We recommend to use this feature for visualization and qualitative analyses
              (Check the MoltenProt User Documentation at the main page - https://spc.embl-hamburg.de).",
              placement = "right")
          )
        ),

        column(3,
          p(HTML("<b>WL1</b>"),
            selectInput(
              "wl_for_ratio_1",
              NULL,
              choices = rev(dsf$signals),
              selectize = FALSE
            )
          )
        ),

        column(3,
          p(HTML("<b>WL2</b>"),
            selectInput(
              "wl_for_ratio_2",
              NULL,
              choices = dsf$signals[
                !grepl('Ratio', dsf$signals) &
                  !grepl('BCM', dsf$signals) &
                    !grepl('SVD', dsf$signals)
              ],
              selectize = FALSE
            )
          )
        )

      ),

      footer=tagList(
        #actionButton('submitConfig', 'Submit'),
        modalButton('Close')
      )

    )
  )

}

observeEvent(input$wl_range,{
  reactives$min_wl <- input$wl_range[1]
  reactives$max_wl <- input$wl_range[2]
})

observeEvent(input$show_full_spectra_menu, {
  full_spectrum_dialog()
})

renderSpectralPlots <- function() {

  # Delete all previous Tabs
  if (!is.null(reactives$spectra_panel_names)) {
    for (tabPanelTargetName in reactives$spectra_panel_names) {
      removeTab(inputId = "tabset1",
                target = tabPanelTargetName)
    }
  }

  nConditions   <- length(dsf$conditions)

  # Return NULL if no conditions are selected
  if (nConditions == 0) return(NULL)

  tabPanelNames <- generate_tab_panels(nConditions)

  # Create new Tabs - do not use a for loop here, because it does not work
  lapply(tabPanelNames, function(sp){

    appendTab(inputId = "tabset1",
              tabPanel(sp,plotOutput(sp))
    )
    return(NULL)
  })

  reactives$spectra_panel_names <- tabPanelNames

  all_signals <- dsf$signal_data_dictionary
  all_temps   <- dsf$temp_data_dictionary

  tog <- join_all_signals(all_signals,all_temps,
                          c(dsf$conditions),reactives$include_vector,
                          input$sg_range[1],input$sg_range[2])

  maxPanels     <- length(tabPanelNames)

  # Use lapply over for loop to avoid only rendering the last plot
  lapply(0:(maxPanels-1), function (i) {

    tabPanelName <- tabPanelNames[i+1]

    tog2      <- tog[tog$name %in% unique(tog$name)[(20*i+1):(20*(i+1))],]
    tog2$name <- factor(tog2$name,levels=unique(tog2$name))

    fig       <- plot_whole_spectra(
        tog2,font_size = input$plot_axis_size,
        min_wl=input[['wl_range']][1],max_wl=input[['wl_range']][2])

    output[[tabPanelName]] <- renderPlot(fig)

    return(NULL)
  })

  return(NULL)
}

observeEvent(input$updateSpectralPlots,{

  renderSpectralPlots()

})

observeEvent(input$generateRatioSignal,{

  wl1 <- input$wl_for_ratio_1
  wl2 <- input$wl_for_ratio_2

  dsf$create_ratio_signal(wl1, wl2)

  updateSelectInput(session, "which",choices  = dsf$signals)

})

observeEvent(input$decomposeSpectra,{

  withBusyIndicatorServer("spectraDecompositionHiddenButton",{

    dsf$decompose_spectra(
      input$wl_range[1],
      input$wl_range[2]
    )

  })

  updateSelectInput(session, "which",choices  = dsf$signals)

})
observeEvent(input$unify_signals, {

    # Open a modal dialog to select signals to unify
    showModal(modalDialog(
        size = "m",
        easyClose = TRUE,

        # selectInput to select signals to unify
        selectInput("which_to_unify", "Select signals to unify:",
                    choices = dsf$all_signals,
                    multiple = TRUE
        ),

        textInput("new_signal_name", "New unified signal name:", "Signal"),

        footer=tagList(
            actionButton('submitUnify', 'Submit'),
            modalButton('Cancel')
        )

    ))

})

observeEvent(input$submitUnify, {

    req(input$which_to_unify)
    req(input$new_signal_name)

    # Close the modal dialog
    removeModal()

    result <- tryCatch(
        {

            # Perform the unification of signals
            dsf$unify_signals(
                signal_list = input$which_to_unify,
                new_signal_name = input$new_signal_name
            )

        }, error = function(e) {
            if (inherits(e, "python.builtin.ValueError")) {
            err <- py_last_error()
            shinyalert(
                title = "Error",
                text = paste0("⚠ Processing error: ", err$value),
                type = "error"
            )
            return('Error')
            } else {
            stop(e) # rethrow non-Python errors
            }
        }
    )


    # Update the selectInput choices for future unifications
    updateSelectInput(
        session, "which",
        choices = dsf$all_signals
    )

})
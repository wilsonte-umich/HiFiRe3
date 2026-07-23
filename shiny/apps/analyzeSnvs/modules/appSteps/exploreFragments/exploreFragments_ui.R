#----------------------------------------------------------------------
# UI components for the exploreFragments appStep module
#----------------------------------------------------------------------

# module ui function
exploreFragmentsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$exploreFragments)

    # return the UI contents
    standardSequentialTabItem(

        # page header text
        options$longLabel,
        options$leaderText,

        # page header links, uncomment as needed
        id = id,
        # documentation = TRUE,
        # terminal = TRUE,
        console = serverEnv$IS_DEVELOPER,
        code = serverEnv$IS_DEVELOPER,
        settings = TRUE,

        # appStep UI elements, populate as needed
        dataSourceTableUI(
            ns("dataSource"), 
            "Data Source", 
            width = 12, 
            collapsible = TRUE
        ),
        fluidRow(
            box(
                title = NULL,
                width = 12,
                solidHeader = FALSE,
                # status = "primary",
                collapsible = FALSE,
                column(
                    width = 2,
                    numericInput(
                        ns("minCoverage"), 
                        "Min Coverage", 
                        value = 10, 
                        min = 0, 
                        max = 50,
                        step = 5
                    )
                ),
                column(
                    width = 2,
                    numericInput(
                        ns("pixelsPerRead"), 
                        "Read Height Pixels", 
                        value = 6, 
                        min = 1, 
                        max = 10,
                        step = 1
                    )
                ),
                column(
                    width = 2,
                    numericInput(
                        ns("pixelsPerBase"), 
                        "Base Width Pixels", 
                        value = 2, 
                        min = 1, 
                        max = 10,
                        step = 1
                    )
                ),
                column(
                    width = 2,
                    numericInput(
                        ns("windowWidthBases"), 
                        "Window Width Bases", 
                        value = 750, 
                        min = 100, 
                        max = 1000,
                        step = 50
                    )
                ),
                column(
                    width = 2,
                    actionButton(ns("newFragment"), "New Fragment", inline = TRUE)
                ),
                column(
                    width = 2,
                    textOutput(ns("fragmentSpan"))
                ),
                NULL
            )
        ),
        fluidRow(
            mdiInteractivePlotBoxUI(
                ns("encodingPlot"), 
                "Fragment Encoding Plot",
                width = 12, 
                collapsible = FALSE,
                solidHeader = FALSE
            )
        ),
        fluidRow(
            bufferedTableUI(
                ns("variantsTable"),
                title = "Fragment Variants Table",
                width = 12,
                collapsible = TRUE,
                collapsed = TRUE
            )
        ),
        NULL
    )
}

#----------------------------------------------------------------------
# UI components for the summarizeDigestions appStep module
#----------------------------------------------------------------------

# module ui function
summarizeDigestionsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$summarizeDigestions)

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
            "Digestion Source", 
            width = 12, 
            collapsible = TRUE
        ),
        fluidRow(
            bufferedTableUI(
                ns("enzymesTable"),
                title = "Blunt Restriction Enzymes",
                width = 12
            )
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("distributionPlot"),
                "Distribution",
                width = 6,
                collapsible = TRUE,
                collapsed = FALSE
            ),
            staticPlotBoxUI(
                ns("nSitesPlot"),
                "Number of Sites",
                width = 6,
                collapsible = TRUE,
                collapsed = FALSE
            )
        ),
        fluidRow(
            staticPlotBoxUI(
                ns("molarFractionPlot"),
                "Molar Fraction",
                width = 6,
                collapsible = TRUE,
                collapsed = FALSE
            ),
            staticPlotBoxUI(
                ns("massFractionPlot"),
                "Mass Fraction",
                width = 6,
                collapsible = TRUE,
                collapsed = FALSE
            )
        ),
        NULL
    )
}

#----------------------------------------------------------------------
# UI components for the exploreJunctions appStep module
#----------------------------------------------------------------------

# module ui function
exploreJunctionsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$exploreJunctions)

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
            mdiInteractivePlotBoxUI(
                ns('sizePlot'), 
                "Junctions By Size",
                data = FALSE,
                width = 8,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE
            ),
            box(
                title = "Selection Summary",
                width = 4,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE,
                fluidRow(
                    style = "margin: 10px; font-size: 1.2em;",
                    actionLink(
                        ns("clearSelections"), 
                        "Clear Selections", 
                        icon = icon("trash-alt")
                    )
                ),
                bufferedTableUI(
                    ns("selectionSummaryTable"),
                    title = "Selection Junction Counts",
                    width = 12,
                    collapsible = TRUE,
                    collapsed = FALSE
                )
            )
        ),
        fluidRow(
            mdiInteractivePlotBoxUI(
                ns('offsetPlotWide'), 
                "Offset (uhom/insertion) Wide",
                data = FALSE,
                width = 6,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE
            ),
            mdiInteractivePlotBoxUI(
                ns('offsetPlotNarrow'), 
                "Offset (uhom/insertion) Narrow",
                data = FALSE,
                width = 6,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE
            )
        ),
        fluidRow(
            bufferedTableUI(
                ns("junctionsTable"),
                title = "Final Junctions, Filtered",
                width = 12,
                collapsible = TRUE
            )
        ),
        fluidRow(
            box(
                title = "Selected Junction",
                width = 12,
                solidHeader = TRUE,
                status = "primary",
                collapsible = TRUE,
                fluidRow(
                    column(
                        width = 12,
                        style = "text-align: center;  margin-top: 5px;",
                        tags$div(
                            actionButton(ns("prevReadI1"), "<", inline = TRUE),
                            tags$div(
                                style = "display: inline-block; margin: 5px; width: 75px;",
                                numericInput(ns("currReadI1"), label = NULL, value = 1)
                            ),
                            actionButton(ns("nextReadI1"), ">", inline = TRUE),
                            textOutput(ns("nExpansionReads"), inline = TRUE)
                        ),
                        tags$p(
                            style = "display: inline-block; text-align: center; margin: 5px;",
                            uiOutput(ns("expandedReadSummary"))
                        )
                    )
                ),
                plotOutput(ns("readQualPlot"), height = 96 * 0.75),
                plotOutput(ns("refAlnPlot"),   height = 96 * 3),
                uiOutput(ns("expandedJunctionSummary")),
                uiOutput(ns("expandedSequences")),
                plotOutput(ns("chrom1Plot"), height = 96 * 1.5),
                plotOutput(ns("chrom2Plot"), height = 96 * 1.5),
            ),
        ),
        NULL
    )
}

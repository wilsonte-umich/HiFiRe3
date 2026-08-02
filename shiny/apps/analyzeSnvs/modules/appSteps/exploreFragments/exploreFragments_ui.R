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
                collapsible = TRUE,
                collapsed = FALSE,
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
                    selectInput(
                        ns("colorPalette"), 
                        "Show Variants", 
                        choices = c(
                            "show_all",
                            "hide_masked",
                            "hide_indels",
                            "hide_masked_and_indels"
                        ),
                        selected  = "show_all"
                    )
                ),
                column(
                    width = 2,
                    checkboxInput(
                        ns("showClonal"), 
                        "Show Clonal", 
                        value = TRUE
                    )
                ),
                NULL
            )
        ),
        fluidRow(
            box(
                title = NULL,
                width = 12,
                solidHeader = FALSE,
                # status = "primary",
                collapsible = TRUE,
                collapsed = FALSE,
                column(
                    width = 2,
                    actionButton(ns("anyFragment"), "Any Fragment", inline = TRUE, width = "100%")
                ),
                column(
                    width = 2,
                    actionButton(ns("singletonSnv"), "Singleton SNV", inline = TRUE, width = "100%")
                ),
                column(
                    width = 2,
                    actionButton(ns("trueSubclonal"), "True Subclonal", inline = TRUE, width = "100%")
                ),
                column(
                    width = 2,
                    actionButton(ns("snv5Read"), "SNV5+ Read", inline = TRUE, width = "100%")
                ),
                column(
                    width = 2,
                    textInput(ns("jumpToFragment"), "Jump to Fragment", width = "100%")
                ),
                column(
                    width = 4,
                    textOutput(ns("fragmentSpan"))
                ),
                NULL        
            )
        ),
        fluidRow(
            mdiInteractivePlotBoxUI(
                ns("readsOnRefPlot"), 
                "Reads On Reference",
                width = 12,
                solidHeader = FALSE, 
                collapsible = TRUE, 
                collapsed = FALSE
            )
        ),
        fluidRow(
            mdiInteractivePlotBoxUI(
                ns("readsOnHap1Plot"), 
                "Reads On Haplotype 1",
                width = 12,
                solidHeader = FALSE, 
                collapsible = TRUE, 
                collapsed = FALSE
            )
        ),
        fluidRow(
            mdiInteractivePlotBoxUI(
                ns("readsOnHap2Plot"), 
                "Reads On Haplotype 2",
                width = 12,
                solidHeader = FALSE, 
                collapsible = TRUE, 
                collapsed = FALSE
            )
        ),
        fluidRow(
            box(
                title = NULL,
                width = 12,
                solidHeader = FALSE,
                # status = "primary",
                collapsible = TRUE,
                collapsed = FALSE,
                column(
                    width = 2,
                    checkboxInput(
                        ns("tableShowIndels"), 
                        "Show Indels", 
                        value = TRUE
                    )
                ),
                column(
                    width = 2,
                    checkboxInput(
                        ns("tableShowClonal"), 
                        "Show Clonal", 
                        value = TRUE
                    )
                ),
                column(
                    width = 2,
                    checkboxInput(
                        ns("tableShowQnames"), 
                        "Show QNames", 
                        value = TRUE
                    )
                ),
                NULL
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

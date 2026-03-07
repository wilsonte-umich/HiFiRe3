#----------------------------------------------------------------------
# UI components for the exploreVariants appStep module
#----------------------------------------------------------------------

# module ui function
exploreVariantsUI <- function(id, options) {

    # initialize namespace
    ns <- NS(id)
    
    # override missing options to module defaults
    options <- setDefaultOptions(options, stepModuleInfo$exploreVariants)

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
                status = "primary",
                collapsible = FALSE,
                column(
                    width = 10,
                    fluidRow(
                        column(
                            style = "text-align: right;",
                            width = 5,
                            hf3_flagFilterModeUI(ns("alnFlagPassMode"), "pass", "Aln Flag Mode")
                        ),
                        column(
                            width = 7,
                            hf3_alnFailureBitsUI(ns("alnFlagPassBits"), TRUE)
                        )
                    ),
                )
            )
        ),
        NULL
    )
}

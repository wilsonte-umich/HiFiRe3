#----------------------------------------------------------------------
# server components for the exploreVariants appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
exploreVariantsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'exploreVariants'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    settings = file.path(app$sources$suiteGlobalDir, "settings", "jxn_filters.yml"), #id, # for step-level settings
    size = "m"
)
dpi <- 96

#----------------------------------------------------------------------
# load SV junctions
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("dataSource", selection = "single") 
# junctions_all <- reactive({ # all junction from data source, prior to filtering
#     sourceId <- sourceId()
#     req(sourceId)
#     startSpinner(session, message = "loading junctions")
#     hf3_getJunctions_all_source(sourceId)[
#         jxn_type == hf3_junctions$typeToBits["translocation"] | 
#         sv_size > 0 # suppress 0bp inversions
#     ]
# })
# junctions_filtered <- reactive({ # all junctions, filtered by top-level settings
#     junctions_all <- junctions_all()
#     req(junctions_all)
#     startSpinner(session, message = "loading filtered junctions")
#     junctions_all %>% hf3_applyJunctionFilters(settings, input)
# })

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
    # # updateSelectInput(session, "sampleSet-sampleSet", selected = bm$input[['sampleSet-sampleSet']])
    if(!is.null(bm$outcomes)) {
    #     # outcomes <<- listToReactiveValues(bm$outcomes)
        sizePlot$settings$replace(bm$outcomes$sizePlotSettings)
        # offsetPlotWide$settings$replace(bm$outcomes$offsetPlotWideSettings)
        # offsetPlotNarrow$settings$replace(bm$outcomes$offsetPlotNarrowSettings)
    }
    bookmarkObserver$destroy()
})

#----------------------------------------------------------------------
# set return values as reactives that will be assigned to app$data[[stepName]]
#----------------------------------------------------------------------
list(
    input = input,
    settings = settings$all_,
    outcomes = reactive({ list(
        # sizePlotSettings = sizePlot$settings$all_(),
        # offsetPlotWideSettings = offsetPlotWide$settings$all_(),
        # offsetPlotNarrowSettings = offsetPlotNarrow$settings$all_()
    ) }),
    settingsObject = settings,
    # junctions_filtered = junctions_filtered,
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------

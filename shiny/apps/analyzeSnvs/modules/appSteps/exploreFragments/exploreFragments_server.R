#----------------------------------------------------------------------
# server components for the exploreFragments appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
exploreFragmentsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'exploreFragments'
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
encodingBaseColors <- list( # generally follow IGV base color conventions
    M = rgb(0.75, 0.75, 0.75),   # any base match = light grey
    "=" = rgb(0.75, 0.75, 0.75),   # any base match = light grey

    A = rgb(0,    0.8,    0),    # green  
    C = rgb(0,    0,    1),    # blue
    G = rgb(0.82, 0.43, 0),      # orange
    T = rgb(1,    0,    0),    # red
    N = rgb(0.9, 0.9, 0.9),      # N, treated as M
    
    a = rgb(0.5,    0.8,    0.5),    # green  
    c = rgb(0.6,    0,    0.6),    # blue
    g = rgb(0.82, 0.63, 0.4),      # orange
    t = rgb(1,    0.6,    0.6),    # red
    n = rgb(0.9, 0.9, 0.9),      # N, treated as M

    D = rgb(0.1, 0.1, 0.1),      # deleted/missing = black
    I = rgb(0.75,   0,    0.75),  # insertion = purple
    "-" = rgb(0.1, 0.1, 0.1),      # deleted/missing = black
    "+" = rgb(0.75,   0,    0.75),  # insertion = purple
    "?" = rgb(0.9, 0.9, 0.9)
)

#----------------------------------------------------------------------
# load data
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("dataSource", selection = "single") 
encodings <- reactive({
    sourceId <- req(sourceId())
    startSpinner(session, message = "loading encodings")
    hf3_getEncodings_cached(sourceId)
})
variants <- reactive({
    sourceId <- req(sourceId())
    startSpinner(session, message = "loading variants")
    hf3_getVariants_cached(sourceId)
})
fragment <- reactive({
    x <- input$newFragment
    encodings <- req(encodings())
    variants <- req(variants())
    startSpinner(session, message = "selecting fragment")

    encoding <- encodings[n_reads >= input$minCoverage][sample(.N, 1)]

    # snps <- variants[
    #     coverage >= input$minCoverage & 
    #     is_snp == TRUE & 
    #     any_allowed == TRUE & 
    #     simple_repeat == 0 &
    #     count == 1
    # ]
    # snps <- snps[sample(.N, 1)]
    # encoding <- encodings[
    #     snps$chrom_index1 == chrom_index1 &
    #     snps$start0 >= start0 &
    #     snps$start0 <= end1
    # ]

    variants <- variants[
        chrom_index1 == encoding$chrom_index1 &
        start0 >= encoding$start0 &
        start0 <= encoding$end1
    ]
    stopSpinner(session)
    list(
        encoding = encoding,
        variants = variants
    )
})
output$fragmentSpan = renderText({
    sourceId <- req(sourceId())
    fragment <- req(fragment())
    dstr(fragment)
    chrom <- hf3_getChromNames(sourceId, fragment$encoding$chrom_index1)
    dstr(chrom)
    paste(
        chrom,
        fragment$encoding$start0,
        fragment$encoding$end1
    )
})

#----------------------------------------------------------------------
# fragment encoding plot
#----------------------------------------------------------------------
encodingPlot <- mdiInteractivePlotBoxServer(
    "encodingPlot",
    # click = TRUE,
    # brush = TRUE,
    points  = TRUE, # set to TRUE to expose relevant plot options
    lines   = TRUE,
    settings = NULL, # an additional settings template as a list()
    defaults = NULL, # list of default settings values use to inialize settings
    create = function(...) {
        fragment <- req(fragment())
        d <- fragment$encoding
        plotOrder <- if (sum(fragment$variants$any_allowed == 1) >= 1){
            heterozygous <- fragment$variants[any_allowed == 1][which.min(abs(vaf - 0.5))]
            het_qnames <- strsplit(heterozygous$qnames, ",")[[1]]
            qnames <- strsplit(d$qnames, ",")[[1]]
            is_het <- qnames %in% het_qnames
            c(which(is_het), which(!is_het))
        } else 1:d$n_reads

        dpi <- 96
        pointsize <- 7
        px_per_read <- 5
        px_per_base <- 2

        nStackRows <- floor(d$n_bases / input$windowWidthBases) + 1
        nStrackTracks <- d$n_reads * 2 + 1
        ymax <- nStackRows * nStrackTracks

        width_pixels  <- input$pixelsPerBase * input$windowWidthBases
        height_pixels <- input$pixelsPerRead * nStrackTracks * nStackRows
        layout <- list(
            width     = width_pixels,
            height    = height_pixels,
            pointsize = pointsize,
            dpi       = dpi
        )
        
        png(file = encodingPlot$pngFile, width = width_pixels, height = height_pixels, units = "px", 
            pointsize = pointsize, res = dpi, type = "cairo")
        par(mar = c(0, 0, 0, 0))

        encodingPlot$initializeFrame(
            layout,
            xlim = c(0, input$windowWidthBases),
            ylim = c(0, ymax),
            xaxs = "i",
            yaxs = "i",
            xaxt = "n",
            yaxt = "n"
        )

        readStart0s <- as.integer(strsplit(d$read_start0s, ",")[[1]])
        encodings   <-            strsplit(d$encodings, ",")[[1]]
        insertions  <-            strsplit(d$insertions, ",")[[1]]
        seq         <-            strsplit(d$seq, "")[[1]]

        rect(
            xleft   = 0, 
            xright  = input$windowWidthBases, 
            ybottom = 0, 
            ytop    = ymax, 
            col     = encodingBaseColors$M, 
            border  = NA
        ) 

        for (b1 in 1:d$n_bases) {
            b0 <- b1 - 1
            trackRow1 <- floor(b0 / input$windowWidthBases) + 1
            plotTrackRow0 <- nStackRows - trackRow1
            x0 <- b0 %% input$windowWidthBases
            y0 <- plotTrackRow0 * nStrackTracks + d$n_reads
            rect(
                xleft   = x0, 
                xright  = x0 + 1, 
                ybottom = y0,
                ytop    = y0 + 1, 
                col     = encodingBaseColors[[seq[b1]]], 
                border  = NA
            ) 
        }

        for (plotR1 in 1:d$n_reads){
            dataR1 <- plotOrder[plotR1]
            b0 <- readStart0s[dataR1] - d$start0 # can be negative on first trackRow
            encoding <- strsplit(encodings[dataR1], "")[[1]]
            e1 <- 1
            while (e1 <= length(encoding)){
                op <- encoding[e1]
                if (op == "="){
                    nMatch <- ""
                    while (e1 <= length(encoding) && grepl("[0-9]", encoding[e1 + 1])){
                        nMatch <- paste0(nMatch, encoding[e1 + 1])
                        e1 <- e1 + 1
                    }
                    b0 <- b0 + as.integer(nMatch)
                } else {
                    trackRow1 <- floor(b0 / input$windowWidthBases) + 1
                    plotTrackRow0 <- nStackRows - trackRow1
                    x0 <- b0 %% input$windowWidthBases
                    y0 <- plotTrackRow0 * nStrackTracks + d$n_reads + plotR1
                    rect(
                        xleft   = x0, 
                        xright  = x0 + 1, 
                        ybottom = y0,
                        ytop    = y0 + 1, 
                        col     = encodingBaseColors[[op]], 
                        border  = NA
                    ) 
                    b0 <- b0 + 1                
                }
                e1 <- e1 + 1
            }

            b0 <- readStart0s[dataR1] - d$start0
            insertion <- strsplit(insertions[dataR1], "")[[1]]
            i1 <- 1
            while (i1 <= length(insertion)){
                op <- insertion[i1]
                if (op == "="){
                    nMatch <- ""
                    while (i1 <= length(insertion) && grepl("[0-9]", insertion[i1 + 1])){
                        nMatch <- paste0(nMatch, insertion[i1 + 1])
                        i1 <- i1 + 1
                    }
                    b0 <- b0 + as.integer(nMatch)
                } else {
                    trackRow1 <- floor(b0 / input$windowWidthBases) + 1
                    plotTrackRow0 <- nStackRows - trackRow1
                    x0 <- b0 %% input$windowWidthBases
                    y0 <- plotTrackRow0 * nStrackTracks + plotR1
                    rect(
                        xleft   = b0, 
                        xright  = b0 + 1, 
                        ybottom = y0,
                        ytop    = y0 + 1, 
                        col     = encodingBaseColors[[op]], 
                        border  = NA
                    ) 
                    b0 <- b0 + 1
                }
                i1 <- i1 + 1
            }
        } 

        abline(h = 0:ymax, lwd = 0.5, col = rgb(0.5,0.5,0.5))  
        abline(h = 0:nStackRows * nStrackTracks, lwd = 0.5, col = "black") 

        encodingPlot$finishPng(layout)
    }
)

#----------------------------------------------------------------------
# fragment variants tables
#----------------------------------------------------------------------
# observeEvent(input$clearSelections, {
#     startSpinner(session, message = "clearing selections")
#     selectedJxnI$sizePlot   <- c()
#     selectedJxnI$offsetPlot <- c()
#     sizePlotSelection(NULL)
#     offsetPlotSelection(NULL)
#     stopSpinner(session)
# })
variantsTableData <- reactive({
    fragment <- req(fragment())
    d <- fragment$variants[
        is_snp == TRUE & 
        any_allowed == TRUE & 
        simple_repeat == FALSE
    ]
    d$qnames <- NULL
    d$vaf <- round(d$vaf, 3)
    d$vaf_bin <- NULL
    d
})
variantsTable <- bufferedTableServer(
    "variantsTable",
    id,
    input,
    variantsTableData,
    selection = 'single',
    selectionFn = function(selectedRows) NULL,
    options = list()
)

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
        # sizePlot$settings$replace(bm$outcomes$sizePlotSettings)
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

#----------------------------------------------------------------------
# server components for the summarizeLibraries appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
summarizeLibrariesServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'summarizeLibraries'
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
    A = rgb(0,    0.8,    0),    # green  
    C = rgb(0,    0,    1),    # blue
    G = rgb(0.82, 0.43, 0),      # orange
    T = rgb(1,    0,    0),    # red
    N = rgb(0.5, 0.5, 0.5),      # N, treated as M
    D = rgb(0.1, 0.1, 0.1),      # deleted/missing = black
    I = rgb(0.75,   0,    0.75)  # insertion = purple
)

#----------------------------------------------------------------------
# load data
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("dataSource", selection = "single") 
variants <- reactive({
    sourceId <- req(sourceId())
    startSpinner(session, message = "loading variants")
    hf3_getVariants_cached(sourceId)
})
clonal_encodings <- reactive({
    sourceId <- req(sourceId())
    startSpinner(session, message = "loading clonal encodings")
    hf3_getEncodings_cached(sourceId, "clonal")
})

#----------------------------------------------------------------------
# fragment length distribution plot
#----------------------------------------------------------------------
fragLengthPlot <- mdiInteractivePlotBoxServer(
    "fragLength",
    # click = TRUE,
    # brush = TRUE,
    points  = TRUE, # set to TRUE to expose relevant plot options
    lines   = TRUE,
    settings = NULL, # an additional settings template as a list()
    defaults = NULL, # list of default settings values use to inialize settings
    create = function(...) {
        sourceId <- req(sourceId())
        sample_bits <- hf3_sample_bits(sourceId)
        d <- req(clonal_encodings())
        startSpinner(session, message = "rendering length distribution")
        ymax <- 0
        d <- lapply(sample_bits, function(sample_bit){
            dd <- d[bitwAnd(sample_bits, sample_bit) > 0]
            dd <- dd[, 
                .(n_reads = sum(as.integer(unlist(strsplit(sample_bitss, ","))) == sample_bit)), 
                keyby = .(n_bases_bin)
            ]
            dd[, freq := n_reads / sum(n_reads)]
            ymax <<- max(ymax, dd$freq)
            dd
        })
        layout <- fragLengthPlot$initializePng(mar = c(4.1, 4.1, 0.9, 0.9)) %>% 
                  fragLengthPlot$initializeFrame(
            xlim = c(1000, 9000),
            ylim = c(0, ymax * 1.05),
            xlab = "RE Fragment Length",
            ylab = "Frequency (by sample)",
            xaxs = "i",
            yaxs = "i"
        )
        bar_sep <- 250 / length(d)
        offset <- length(d) / 2 * bar_sep
        lapply(1:length(d), function(i){
            dd <- d[[i]]
            fragLengthPlot$addPoints(
                x = dd$n_bases_bin - offset + (i - 1) * bar_sep,
                y = dd$freq,
                typ = "h",
                lwd = 1.5,
                col = i
            )
        })
        stopSpinner(session)
        fragLengthPlot$finishPng(layout)
    }
)

#----------------------------------------------------------------------
# fragment coverage distribution plot
#----------------------------------------------------------------------
fragCoveragePlot <- mdiInteractivePlotBoxServer(
    "fragCoverage",
    # click = TRUE,
    # brush = TRUE,
    points  = TRUE, # set to TRUE to expose relevant plot options
    lines   = TRUE,
    settings = NULL, # an additional settings template as a list()
    defaults = NULL, # list of default settings values use to inialize settings
    create = function(...) {
        d <- req(clonal_encodings())
        startSpinner(session, message = "rendering coverage distribution")
        d <- d[, .(count = .N), keyby = .(n_reads)]
        d[, freq := count / sum(count)]
        layout <- fragCoveragePlot$initializePng(mar = c(4.1, 4.1, 0.9, 0.9)) %>% 
                  fragCoveragePlot$initializeFrame(
            xlim = c(0, d[, max(n_reads)]),
            ylim = c(0, d[, max(freq) * 1.05]),
            xlab = "RE Fragment Coverage",
            ylab = "Frequency",
            xaxs = "i",
            yaxs = "i"
        )
        abline(v = input$minCoverage, col = "blue", lwd = 1.5)
        fragCoveragePlot$addPoints(
            x = d$n_reads,
            y = d$freq,
            typ = "h"
        )
        stopSpinner(session)
        fragCoveragePlot$finishPng(layout)
    }
)

#----------------------------------------------------------------------
# fragment coverage distribution plot
#----------------------------------------------------------------------
lengthVsCoveragePlot <- mdiInteractivePlotBoxServer(
    "lengthVsCoverage",
    # click = TRUE,
    # brush = TRUE,
    points  = TRUE, # set to TRUE to expose relevant plot options
    lines   = TRUE,
    settings = NULL, # an additional settings template as a list()
    defaults = NULL, # list of default settings values use to inialize settings
    create = function(...) {
        d <- req(clonal_encodings())
        startSpinner(session, message = "rendering correlation")
        layout <- lengthVsCoveragePlot$initializePng(mar = c(4.1, 4.1, 0.9, 0.9)) %>% 
                  lengthVsCoveragePlot$initializeFrame(
            xlim = c(1000, 9000),
            ylim = c(0, quantile(d$n_reads, 0.99) * 1.05),
            xlab = "RE Fragment Length",
            ylab = "Read Coverage (all samples)",
            xaxs = "i",
            yaxs = "i"
        )
        d <- d[sample(.N, 5000, replace = FALSE)]
        lengthVsCoveragePlot$addPoints(
            x = jitter(d$n_bases),
            y = jitter(d$n_reads),
            pch = "."
        )
        abline(h = input$minCoverage, col = "blue", lwd = 1.5)
        stopSpinner(session)
        lengthVsCoveragePlot$finishPng(layout)
    }
)

#----------------------------------------------------------------------
# zygosity distribution plot
#----------------------------------------------------------------------
zygosityPlot <- mdiInteractivePlotBoxServer(
    "zygosity",
    # click = TRUE,
    # brush = TRUE,
    points  = TRUE, # set to TRUE to expose relevant plot options
    lines   = TRUE,
    settings = NULL, # an additional settings template as a list()
    defaults = NULL, # list of default settings values use to inialize settings
    create = function(...) {
        d <- req(variants())
        startSpinner(session, message = "rendering zygosity")
        d <- d[
            is_snp == TRUE &
            coverage >= input$minCoverage
        ]
        d <- d[, .(count = .N), keyby = .(vaf_bin)]
        d[, freq := count / sum(count)]
        layout <- zygosityPlot$initializePng(mar = c(4.1, 4.1, 0.9, 0.9)) %>% 
                  zygosityPlot$initializeFrame(
            xlim = c(0, 1),
            ylim = c(0, d[vaf_bin < 1, max(freq) * 1.05]),
            xlab = "Variant Allele Frequency",
            ylab = "Frequency",
            # xaxs = "i",
            yaxs = "i"
        )
        zygosityPlot$addLines(
            x = d$vaf_bin,
            y = d$freq,
            typ = "h",
            lwd = 1
        )
        stopSpinner(session)
        zygosityPlot$finishPng(layout)
    }
)

#----------------------------------------------------------------------
# fragment variants tables
#----------------------------------------------------------------------
variantSummaryTableData <- reactive({
    # fragment <- req(fragment())
    # d <- fragment$variants
    # d$qnames <- NULL
    # d
    sourceId <- req(sourceId())
    smp_bits <- hf3_sample_bits(sourceId)

    clonal_encodings <- req(clonal_encodings())
    variants <- req(variants())

    startSpinner(session, message = "assembling summary table")

    clonal_encodings <- clonal_encodings[n_reads >= input$minCoverage]
    snps <- variants[
        coverage >= input$minCoverage & 
        is_snp == TRUE
    ]
    snps_clonal <- snps[vaf >= 0.2]
    snps_singleton <- snps[n_matching_reads == 1]

    enc_exp <- clonal_encodings[, 
        .(
            n_bases = n_bases,
            sample_bit = as.integer(strsplit(sample_bitss, ",")[[1]])
        ),
        keyby = .(chrom_index1, start0, end1)
    ]

    smp_all <- "all"
    d <- data.table(
        metric = c("n_reads", "n_bases"),
        sample = c(smp_all, smp_all),
        value = clonal_encodings[, c(sum(n_reads), sum(n_reads * n_bases))]
    )
    for (i in 1:length(smp_bits)){
        smp <- paste0("smp", i)
        d <- rbind(d, data.table(
            metric = c("n_reads", "n_bases"),
            sample = c(smp, smp),
            value = enc_exp[, {
                is_smp <- sample_bit == smp_bits[i]
                c(sum(is_smp), sum(n_bases * is_smp))
            }]
        ))
    }

    d <- rbind(d, data.table(
        metric = c("n_snps_clonal", "n_snps_singleton"),
        sample = c(smp_all, smp_all),
        value = c(nrow(snps_clonal), nrow(snps_singleton))
    ))
    for (i in 1:length(smp_bits)){
        smp <- paste0("smp", i)
        d <- rbind(d, data.table(
            metric = c("n_snps_clonal", "n_snps_singleton"),
            sample = c(smp, smp),
            value = c(
                snps_clonal[, sum(bitwAnd(sample_bits, smp_bits[i]) > 0)],
                snps_singleton[, sum(bitwAnd(sample_bits, smp_bits[i]) > 0)]
            )
        ))
    }

    stopSpinner(session)
    dcast(d, metric ~ sample)
})
variantSummaryTable <- bufferedTableServer(
    "variantSummaryTable",
    id,
    input,
    variantSummaryTableData,
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
        fragLengthPlot$settings$replace(bm$outcomes$fragLengthPlotSettings)
        fragCoveragePlot$settings$replace(bm$outcomes$fragCoveragePlotSettings)
        lengthVsCoveragePlot$settings$replace(bm$outcomes$lengthVsCoveragePlotSettings)
        zygosityPlot$settings$replace(bm$outcomes$zygosityPlotSettings)
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
        fragLengthPlotSettings = fragLengthPlot$settings$all_(),
        fragCoveragePlotSettings = fragCoveragePlot$settings$all_(),
        lengthVsCoveragePlotSettings = lengthVsCoveragePlot$settings$all_(),
        zygosityPlotSettings = zygosityPlot$settings$all_()
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

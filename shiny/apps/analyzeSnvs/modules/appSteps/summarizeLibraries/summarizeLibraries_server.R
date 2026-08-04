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

#----------------------------------------------------------------------
# load data
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("dataSource", selection = "single") 
variants <- reactive({
    sourceId <- req(sourceId())
    startSpinner(session, message = "loading variants")
    hf3_getVariants_cached(sourceId)
})
reads_on_reference <- reactive({
    sourceId <- req(sourceId())
    startSpinner(session, message = "loading reads on reference")
    hf3_getFragments_cached(sourceId, "reference")
})
reads_on_haplotype <- reactive({
    sourceId <- req(sourceId())
    startSpinner(session, message = "loading reads on haplotype")
    hf3_getFragments_cached(sourceId, "haplotype")
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
        d <- req(reads_on_reference())
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
        d <- req(reads_on_reference())
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
        d <- req(reads_on_reference())
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
# VAF distribution plot
#----------------------------------------------------------------------
vafPlot <- mdiInteractivePlotBoxServer(
    "vafPlot",
    # click = TRUE,
    # brush = TRUE,
    points  = TRUE, # set to TRUE to expose relevant plot options
    lines   = TRUE,
    settings = NULL, # an additional settings template as a list()
    defaults = NULL, # list of default settings values use to inialize settings
    create = function(...) {
        d <- req(variants())
        startSpinner(session, message = "rendering VAF")
        d <- d[
            is_snv == TRUE &
            n_reads >= input$minCoverage &
            clonal == 1
        ]
        d <- d[, .(count = .N), keyby = .(vaf_bin)]
        d[, freq := count / sum(count)]
        layout <- vafPlot$initializePng(mar = c(4.1, 4.1, 0.9, 0.9)) %>% 
                  vafPlot$initializeFrame(
            xlim = c(0, 1),
            ylim = c(0, d[vaf_bin < 1, max(freq) * 1.05]),
            xlab = "Variant Allele Frequency",
            ylab = "Frequency",
            # xaxs = "i",
            yaxs = "i"
        )
        vafPlot$addLines(
            x = d$vaf_bin,
            y = d$freq,
            typ = "h",
            lwd = 1
        )
        stopSpinner(session)
        vafPlot$finishPng(layout)
    }
)

#----------------------------------------------------------------------
# fragment variants tables
#----------------------------------------------------------------------
variantSummaryTableData <- reactive({
    sourceId <- req(sourceId())
    smp_bits <- hf3_sample_bits(sourceId)
    samples  <- hf3_getSampleNames(sourceId, smp_bits, as_string = FALSE)

    startSpinner(session, message = "loading valid subclonal SNVs")
    valid_subclonal_snvs <- hf3_getValidSubclonalSnvs(sourceId)

    startSpinner(session, message = "loading valid haplotypes")
    valid_haplotypes <- hf3_getValidHaplotypes(sourceId)

    startSpinner(session, message = "loading trustworthy haplotypes")
    trustworthy_haplotypes <- hf3_getTrustworthyHaplotypes(
        sourceId, valid_subclonal_snvs, valid_haplotypes
    )

    startSpinner(session, message = "loading trustworthy subclonal SNVs")
    trustworthy_subclonal_snvs <- hf3_getTrustworthySubclonalSnvs(
        sourceId, valid_subclonal_snvs, trustworthy_haplotypes
    )

    startSpinner(session, message = "loading trustworthy haplotype reads")
    trustworthy_haplotype_reads <- hf3_getHaplotypeReads(
        sourceId, trustworthy_haplotypes
    )

    # assemble a table of metric values
    startSpinner(session, message = "building metrics table")
    smp_all <- "all"
    metrics <- c("n_reads", "n_unmasked_bases", "n_subclonal_snvs")
    n_col <- length(metrics)
    d <- data.table(
        metric     = metrics,
        sample_bit = rep(NA,      n_col),
        sample     = rep(smp_all, n_col),
        value = c(
            trustworthy_haplotype_reads[, c(
                .N, 
                sum(n_unmasked_bases) # includes invariant reads in tally
            )],
            nrow(trustworthy_subclonal_snvs)
        )
    )
    for (i in 1:length(smp_bits)){
        d <- rbind(d, data.table(
            metric     = metrics,
            sample_bit = rep(smp_bits[i], n_col),
            sample     = rep(samples[i],  n_col),
            value = c(
                trustworthy_haplotype_reads[sample_bit == smp_bits[i], c(
                    .N, 
                    sum(n_unmasked_bases)
                )],
                nrow(trustworthy_subclonal_snvs[sample_bits == smp_bits[i]]) # since n_samples == 1
            )
        ))
    }

    # cast to one row per sample, metrics in columns
    d <- dcast(d, sample_bit + sample ~ metric)

    # calculate subclonal SNV rates and return the result
    d[, ":="(
        subclonal_snv_rate = sprintf("%.2e", n_subclonal_snvs / n_unmasked_bases)
    )]
    stopSpinner(session)
    d
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
        # fragLengthPlot$settings$replace(bm$outcomes$fragLengthPlotSettings)
        # fragCoveragePlot$settings$replace(bm$outcomes$fragCoveragePlotSettings)
        # lengthVsCoveragePlot$settings$replace(bm$outcomes$lengthVsCoveragePlotSettings)
        # vafPlot$settings$replace(bm$outcomes$vafPlotSettings)
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
        # fragLengthPlotSettings = fragLengthPlot$settings$all_(),
        # fragCoveragePlotSettings = fragCoveragePlot$settings$all_(),
        # lengthVsCoveragePlotSettings = lengthVsCoveragePlot$settings$all_(),
        # vafPlotSettings = vafPlot$settings$all_()
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

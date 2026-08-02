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
comp <- c("A" = "T", "C" = "G", "G" = "C", "T" = "A")
add_canonical_mutation <- function(dt) {
  dt[, mutation := fcase(
    tgt_bases %in% c("C", "T"), paste0(     tgt_bases,  ">",      alt_bases),
    tgt_bases %in% c("A", "G"), paste0(comp[tgt_bases], ">", comp[alt_bases])
  )]
}
variantSummaryTableData <- reactive({
    sourceId <- req(sourceId())
    smp_bits <- hf3_sample_bits(sourceId)
    reads_on_hap <- req(reads_on_haplotype())
    variants <- req(variants())

    startSpinner(session, message = "assembling summary table")

    # haplotype FILTER 1: reject heterozygous haplotypes with too few reads to establish a consensus
    startSpinner(session, message = "applying read filter 1")
    message(paste(nrow(reads_on_hap), " = number of fragment haplotypes"))
    print(reads_on_hap[, .(
        .N, 
        n_bases          = sprintf("%2e", sum(n_bases)), 
        n_unmasked_bases = sprintf("%2e", sum(n_unmasked_bases))
    ), keyby = .(haplotype)])
    reads_on_hap <- reads_on_hap[n_reads >= 3]
    message(paste(nrow(reads_on_hap), " = number of allowed haplotypes, n_reads >= 3"))
    print(reads_on_hap[, .(
        .N, 
        n_bases          = sprintf("%2e", sum(n_bases)), 
        n_unmasked_bases = sprintf("%2e", sum(n_unmasked_bases))
    ), keyby = .(haplotype)])
    print(reads_on_hap[n_reads <= 65, .N, keyby = .(haplotype, n_reads)] %>% dcast(n_reads ~ haplotype))

    # TODO: filter away problematic fragment haplotypes?

    # expand haplotypes to one read per row with metadata for further filtering and grouping
    startSpinner(session, message = "expanding reads")
    reads_on_hap_expanded <- reads_on_hap[, 
        .(
            sample_bit      = as.integer(strsplit(sample_bitss, ",")[[1]]),
            qname           =            strsplit(qnames, ",")[[1]],
            n_read_variants = as.integer(strsplit(n_variantss, ",")[[1]])
        ),
        keyby = .(
            chrom_index1, 
            start0, 
            end1, 
            haplotype, 
            n_unmasked_bases
        )
    ]
    message(paste(nrow(reads_on_hap_expanded), " = number of reads in allowed haplotypes"))
    print(reads_on_hap_expanded[n_read_variants <= 10, .N, keyby = .(n_read_variants)])

    # haplotype FILTER 2: filter to reads suitable for SNV calling
    startSpinner(session, message = "applying read filter 2")
    reads_on_hap_expanded <- reads_on_hap_expanded[
        haplotype != 3 |     # all heterozygous haplotype read are informative
        n_read_variants <= 1 # cannot trust multi-variant homozogyous reads
    ]
    message(paste(nrow(reads_on_hap_expanded), " = number of allowed reads in allowed haplotypes"))
    print(reads_on_hap_expanded[n_read_variants <= 10, .N, keyby = .(n_read_variants)])

    # variant FILTER 1: restrict to subclonal SNVs with the same coverage threshold as above
    startSpinner(session, message = "applying variant filter 1")
    snvs_subclonal <- variants[
        is_snv == TRUE & 
        clonal == 0 &
        n_samples == 1 & 
        n_haplotype_reads >= 3 &
        tgt_bases != "N" & 
        alt_bases != "N"
    ]    
    message(paste(nrow(snvs_subclonal), " = number of single-sample subclonal SNVs, n_reads >= 3"))
    print(snvs_subclonal[, .N, keyby = .(n_matching_reads)])

    # variant FILTER 2: filter to reads suitable for SNV calling
    startSpinner(session, message = "applying variant filter 2")
    snvs_subclonal[, n_passing_reads := fcase(
        haplotype == 3,
        n_matching_reads - n_multivariant_reads, # matches read FILTER 2 above
        default = n_matching_reads
    )]
    print(snvs_subclonal[, .N, keyby = .(n_passing_reads)])
    snvs_subclonal <- snvs_subclonal[n_passing_reads > 0]
    snvs_singleton <- snvs_subclonal[n_passing_reads == 1]
    message(paste(nrow(snvs_subclonal), " = number of allowed single-sample subclonal SNVs"))
    print(snvs_subclonal[, .N, keyby = .(tgt_bases, alt_bases)] %>% dcast(tgt_bases ~ alt_bases))
    add_canonical_mutation(snvs_subclonal)
    print(snvs_subclonal[, .N, keyby = .(mutation, sample_bits)] %>% dcast(mutation ~ sample_bits))

    smp_all <- "all"
    d <- data.table(
        metric = c("n_reads", "n_bases"),
        sample_bit = c(NA, NA),
        sample = c(smp_all, smp_all),
        value = reads_on_hap_expanded[, c(.N, sum(n_unmasked_bases))]
    )
    for (i in 1:length(smp_bits)){
        d <- rbind(d, data.table(
            metric = c("n_reads", "n_bases"),
            sample_bit = c(smp_bits[i], smp_bits[i]),
            sample = hf3_getSampleNames(sourceId, c(smp_bits[i], smp_bits[i]), as_string = FALSE),
            value = reads_on_hap_expanded[
                sample_bit == smp_bits[i], 
                c(.N, sum(n_unmasked_bases))
            ]
        ))
    }

    d <- rbind(d, data.table(
        metric = c("n_snvs_subclonal", "n_snvs_singleton"),
        sample_bit = c(NA, NA),
        sample = c(smp_all, smp_all),
        value = c(nrow(snvs_subclonal), nrow(snvs_singleton))
    ))
    for (i in 1:length(smp_bits)){
        d <- rbind(d, data.table(
            metric = c("n_snvs_subclonal", "n_snvs_singleton"),
            sample_bit = c(smp_bits[i], smp_bits[i]),
            sample = hf3_getSampleNames(sourceId, c(smp_bits[i], smp_bits[i]), as_string = FALSE),
            value = c(
                snvs_subclonal[, sum(bitwAnd(sample_bits, smp_bits[i]) > 0)],
                snvs_singleton[, sum(bitwAnd(sample_bits, smp_bits[i]) > 0)]
            )
        ))
    }

    stopSpinner(session)
    d <- dcast(d, sample_bit + sample ~ metric)
    d[, ":="(
        subclonal_rate = sprintf("%.2e", n_snvs_subclonal / n_bases),
        singleton_rate = sprintf("%.2e", n_snvs_singleton / n_bases)
    )]
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

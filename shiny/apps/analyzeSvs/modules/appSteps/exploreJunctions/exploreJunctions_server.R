#----------------------------------------------------------------------
# server components for the exploreJunctions appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
exploreJunctionsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'exploreJunctions'
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
jxnFileType <- "svJunctions1Bgz"
dpi <- 96

#----------------------------------------------------------------------
# load SV junctions
#----------------------------------------------------------------------
sourceIds <- dataSourceTableServer("dataSource", selection = "multiple") 
junctions_all <- reactive({ # all junction from data source, prior to filtering
    sourceIds <- sourceIds()
    req(sourceIds)
    startSpinner(session, message = "loading junctions")
    af_getJunctions_all_sources(sourceIds)[
        jxnType == af_junctions$typeToIndex$translocation | 
        svSize > 0 # suppress 0bp inversions
    ]
})
junctions_filtered <- reactive({ # all junctions, filtered by top-level settings
    junctions_all <- junctions_all()
    req(junctions_all)
    startSpinner(session, message = "loading filtered junctions")
    junctions_all %>% af_applyJunctionFilters(settings)
})
# readPaths <- reactive({
#     sourceIds <- sourceIds()
#     req(sourceIds)
#     x <- af_getReadPaths_sources(sourceIds)
#     stopSpinner(session)
#     x
# })

#----------------------------------------------------------------------
# count SV junctions
#----------------------------------------------------------------------
count_junctions_by_type <- function(jxns){
    agg <- jxns[, .N, keyby = .(jxnType, isIntergenome)]
    agg[, jxnTypeLabel := af_junctions$getJxnTypeLabels(jxnType)]
    agg[isIntergenome == TRUE, ":="(jxnTypeLabel = "int-gen", jxnType = 5L)]
    agg[, color := ifelse(
        isIntergenome,
        CONSTANTS$plotlyColors$purple,
        af_junctions$getIndexColor(jxnType)
    )]
    agg[order(jxnType)]
}
count_junctions_by_class <- function(jxns){
    agg <- jxns[, .N, keyby = .(jxnClass)]
    agg[, jxnClassName := af_junctions$getJxnClassLabels(jxnClass)]
    agg[, color := af_junctionClasses$getClassColor(jxnClass)]
    agg
}
count_junctions_by_both <- reactive({
    jxns <- junctions_selected_both()
    req(nrow(jxns) > 0)
    agg <- jxns[, .N, keyby = .(jxnType, isIntergenome, jxnClass)]
    agg[, jxnTypeLabel := af_junctions$getJxnTypeLabels(jxnType)]
    agg[isIntergenome == TRUE, ":="(jxnTypeLabel = "int-gen", jxnType = 5L)]
    agg[, jxnClassName := af_junctionClasses$getJxnClassLabels(jxnClass)]
    dcast(agg, jxnClassName ~ jxnTypeLabel, value.var = "N", fun.aggregate = sum, fill = 0)
})

#----------------------------------------------------------------------
# SV junction selections
#   sizePlot is filtered by offset plot selections only
#   offsetPlots are filtered by size plot selections only
#   downstream tables and plots are filtered by both, requiring at least one to be selected
#----------------------------------------------------------------------
sizePlotSelection   <- reactiveVal(NULL)
offsetPlotSelection <- reactiveVal(NULL) # only one offset selection is applied to both plots, from either plot
selectedJxnI  <- reactiveValues(
    sizePlot   = c(),
    offsetPlot = c()
)
junctions_selected_both <- reactive({
    req(length(selectedJxnI$sizePlot) > 0 || length(selectedJxnI$offsetPlot) > 0)
    junctions_filtered()[
        jxnI %in% if(length(selectedJxnI$sizePlot) > 0 && length(selectedJxnI$offsetPlot) > 0){
            intersect(selectedJxnI$sizePlot, selectedJxnI$offsetPlot)
        } else {
            c(selectedJxnI$sizePlot, selectedJxnI$offsetPlot)
        }
    ]
})

#----------------------------------------------------------------------
# SV size distributions
#----------------------------------------------------------------------
sizePlotJunctions <- reactive({ # all junctions, prior to filtering for offset selection
    jxns <- junctions_filtered()
    req(nrow(jxns) > 0)
    startSpinner(session, message = "grouping SVs")
    samples <- unique(jxns$sample)
    jxns[, stratumLabel := af_junctionStrata$getJxnStratumLabels(jxnStratum)]
    stratumLabels <- af_junctionStrata$indexToStratumLabel
    strata <- 1:length(stratumLabels)
    names(strata) <- stratumLabels
    jitterAmount <- 0.4
    typeJitterAmount <- jitterAmount / 3
    jxns <- jxns[, {
        stratumI <- strata[stratumLabel]
        jxnTypeI <- unlist(af_junctions$typeToIndex[jxnType + 1])
        jxnTypeSubRowOffset <- ifelse(
            isIntergenome,
            1/3,
            af_junctions$indexToRowOffset[jxnTypeI]
        )
        workingLog10Size <- ifelse(
            jxnType == af_junctions$typeToIndex$translocation, 
            jitter(rep(9.25, .N), amount = 0.75), 
            svSize %>% log10
        )
        pointColor <- ifelse(
            isIntergenome,
            CONSTANTS$plotlyColors$purple,
            af_junctions$getIndexColor(jxnType)
        )
        .(
            stratumI = stratumI,
            jxnI = jxnI, # for interactive selection and mapping to source data.table
            x = workingLog10Size,
            y_ = stratumI + jxnTypeSubRowOffset, 
            color = pointColor,
            jxnType = jxnType,
            isIntergenome = isIntergenome,
            hasExcluded = isExcluded_1 == 1 | isExcluded_2 == 1,
            hasAltAlignment = hasAltAlignment == 1
        )
    }, by = .(stratumLabel)]
    jxns[, y := jitter(y_, amount = typeJitterAmount)] # spread the subrow points through a y-axis jitter
    list(
        samples       = samples,
        jxns          = jxns,
        stratumLabels = stratumLabels,
        strata        = strata,
        nStrata       = length(strata),
        jitterAmount  = 1/3 + typeJitterAmount,
        labelHeight   = 2 * (1/3 + typeJitterAmount) / 5
    )
})
createSizePlot <- function(settings, plot) {
    d <- sizePlotJunctions()
    if(length(selectedJxnI$offsetPlot) > 0) {
        d$jxns <- d$jxns[jxnI %in% selectedJxnI$offsetPlot]
    }
    if(nrow(d$jxns) == 0) {
        stopSpinner(session)
        req(FALSE)
    }
    startSpinner(session, message = "plotting sizes")

    # initialize plot
    lwd  <- settings$get("Points_and_Lines","Line_Width")
    xlim <- c(-0.25, 10.25) # 1e9 to 1e10 used for plotting translocations
    ylim <- c(0.45, d$nStrata + 0.55)
    layout <- plot$initializePng() %>% plot$initializeFrame(
        xlim = xlim,
        ylim = ylim,
        title = d$samples,
        xlab = "log10 SV Size",
        ylab = "",
        yaxt = "n",
        xaxs = "i",
        yaxs = "i"
    )

    # size guide rules
    Alu1_size <- 300
    L1_size   <- 6000
    abline(v = log10(c(1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9)), col = "grey80")
    abline(v = log10(c(Alu1_size, L1_size)), col = CONSTANTS$plotlyColors$purple, lwd = lwd * 1.25)
    abline(h = 1:d$nStrata,       col = "grey60")
    abline(h = 1:d$nStrata + 1/3, col = "grey60")
    abline(h = 1:d$nStrata - 1/3, col = "grey60")

    # add individual junction points
    colors <- ifelse(
        d$jxns$hasExcluded | d$jxns$hasAltAlignment,
        CONSTANTS$plotlyColors$black,
        d$jxns$color
    )
    plot$addPoints(
        x   = d$jxns$x,
        y   = d$jxns$y,
        col = colors %>% addAlphaToColors(settings$get("Junctions_Plot","Point_Alpha"))
    )

    # # overplot violin distributions
    # sizeSpan <- diff(xlim)
    # binSize <- sizeSpan / settings$get("Junctions_Plot","N_Violin_Bins")
    # minObsX <- d$jxns[, min(x, na.rm = TRUE)]
    # d$jxns[, bin := round((x - minObsX) / binSize, 0) * binSize + minObsX]
    # agg <- d$jxns[ , .N, keyby = .(stratumI, bin)] %>%
    #        dcast(bin ~ stratumI, fun.aggregate = sum, fill = 0, value.var = "N")
    # agg <- merge(
    #     data.table(bin = seq(min(d$jxns$bin), max(d$jxns$bin), binSize)),
    #     agg,
    #     by = "bin",
    #     all.x = TRUE
    # )
    # agg[is.na(agg)] <- 0
    # for (i in 1:d$nStrata) {
    #     col <- as.character(i)
    #     N <- sum(agg[[col]])
    #     if(N == 0) next
    #     agg[[col]] <- agg[[col]] / N
    #     agg[[col]] <- agg[[col]] / max(agg[[col]], na.rm = TRUE)
    #     lines(
    #         y = i + agg[[col]] * d$jitterAmount,
    #         x = agg$bin,
    #         col = CONSTANTS$plotlyColors$brown,
    #         lwd = lwd
    #     )
    #     lines(
    #         y = i - agg[[col]] * d$jitterAmount,
    #         x = agg$bin,
    #         col = CONSTANTS$plotlyColors$brown,
    #         lwd = lwd
    #     )
    # }

    # overplot selected area rectangle
    sizePlotSelection <- sizePlotSelection()
    if(!is.null(sizePlotSelection)) {
        rect(
            xleft   = sizePlotSelection$x1,
            xright  = sizePlotSelection$x2,
            ybottom = sizePlotSelection$y1,
            ytop    = sizePlotSelection$y2,
            col     = NA,
            border  = "black",
            lwd     = 1.25,
            lty     = 2
        )
    }

    # stratum axis label
    mtext(
        text = names(d$strata),
        side = 2,
        at   = d$strata,
        line = 0.5,
        las  = 1,
        adj  = 1,
        cex  = 0.9
    )

    # event counts by class and type
    for (i in 1:d$nStrata) {
        agg <- count_junctions_by_type(d$jxns[stratumI == i])
        for (jxnType_ in agg$jxnType) {
            x <- agg[jxnType == jxnType_]
            mtext(
                text = x[, paste(jxnTypeLabel, sprintf("%5d", x$N), sep = " ")],
                side = 4,
                at   = i + d$jitterAmount - d$labelHeight * (x$jxnType - 1),
                line = 0.5,
                las  = 1,
                adj  = 0,
                col  = x$color,
                cex  = 0.8,
                family = "monospace"
            )
        }
    }

    # finish plot and return layout
    plot$finishPng(layout)
}
sizePlot <- mdiInteractivePlotBoxServer(
    'sizePlot',
    click    = FALSE,
    brush    = TRUE,
    create   = function(...) createSizePlot(..., sizePlot), # a function or reactive that creates the plot as a png file using settings and helpers
    points   = TRUE, # set to TRUE to expose relevant plot options
    lines    = TRUE,
    legend   = FALSE,
    margins  = TRUE,
    title    = TRUE,
    data     = FALSE,
    size     = "m",
    settings = list(
        Junctions_Plot = list(
            N_Violin_Bins = list(
                type = "numericInput",
                value = 100,
                min = 25,
                max = 500,
                step = 10
            ),
            Point_Alpha = list(
                type = "numericInput",
                value = 0.6,
                min = 0.1,
                max = 1,
                step = 0.1
            )
        )
    ),
    defaults = list(
        Plot_Frame = list(
            Width_Inches  = 10,
            Height_Inches = 5,
            Font_Size     = 12,
            Bottom_Margin = 4.1,
            Left_Margin   = 6.6,
            Top_Margin    = 1.9,
            Right_Margin  = 6.1
        ),
        Points_and_Lines = list(
            Point_Size = 0.5,
            Point_Type = 19,
            Line_Width = 1.75
        )
    )
)
observeEvent(sizePlot$plot$brush(), {
    coord <- sizePlot$plot$brush()$coord
    minX <- min(coord$x1, coord$x2)
    maxX <- max(coord$x1, coord$x2)
    minY <- min(coord$y1, coord$y2)
    maxY <- max(coord$y1, coord$y2)
    jxnI <- sizePlotJunctions()$jxns[
        x  %between% c(minX, maxX) & 
        y_ %between% c(minY, maxY),
        jxnI
    ]
    sizePlotSelection(list(x1 = minX, x2 = maxX, y1 = minY, y2 = maxY))
    selectedJxnI$sizePlot <- jxnI
})

#----------------------------------------------------------------------
# alignment offset histogram
#----------------------------------------------------------------------
randomOffsetProbs <- c( # as calculate from Roth/Wilson equation
    0.5625,
    0.28125,
    0.10546875,
    0.03515625,
    0.010986328,
    0.003295898,
    0.000961304,
    0.000274658
)
offsetPlotSettings <- list(
    Offset_Plot = list(
        Min_Offset = list(
            type = "numericInput",
            value = -40,
            min = -100,
            max = 0,
            step = 10
        ),
        Max_Offset = list(
            type = "numericInput",
            value = 40,
            min = 0,
            max = 100,
            step = 10
        )
    )
)
createOffsetPlot <- function(settings, plot, v) {

    # get the data
    d <- junctions_filtered()[!is.na(alnOffset)] # may be missing when SEQ was dropped for known SVs
    if(length(selectedJxnI$sizePlot) > 0) {
        d <- d[jxnI %in% selectedJxnI$sizePlot]
    }
    if(nrow(d) == 0) {
        stopSpinner(session)
        req(FALSE)
    }
    startSpinner(session, message = "plotting offsets")

    # aggregate by offset
    d <- d[, .(y = .N), keyby = .(x = alnOffset, jxnTypeI = jxnType + isIntergenome)]
    dd <- d[, .(y = sum(y)), keyby = .(x)]
    xlim <- c(
        settings$get("Offset_Plot","Min_Offset"),
        settings$get("Offset_Plot","Max_Offset")
    )

    # initialize plot
    ylim <- c(0, max(dd$y, na.rm = TRUE) * 1.05)
    layout <- plot$initializePng() %>% plot$initializeFrame(
        xlim = xlim,
        ylim = ylim,
        xlab = "Alignment Offset",
        ylab = "# Junctions",
        xaxs = "i",
        yaxs = "i",
        title = sizePlotJunctions()$samples
    )
    abline(v = v, col = "grey80")
    abline(v = 0, col = "black")
    offsetPlotSelection <- offsetPlotSelection()

    # highlight selected area
    if(!is.null(offsetPlotSelection)) {
        rect(
            xleft   = offsetPlotSelection$x1 - 0.5,
            xright  = offsetPlotSelection$x2 + 0.5,
            ybottom = ylim[1],
            ytop    = ylim[2] / 1.025,
            col     = NA,
            border  = "black",
            lwd     = 1.25,
            lty     = 2
        )
    }

    # add barplot
    d <- d[x %between% xlim]
    if(nrow(d) > 0) {
        d <- dcast(d, x ~ jxnTypeI, value.var = "y", fun.aggregate = sum, fill = 0)
        jxnTypeIs_ <- names(d)[2:ncol(d)]
        jxnTypeIs  <- as.integer(jxnTypeIs_)
        for(i in 1:nrow(d)) {
            ybottom <- 0
            for(jnxTypeI in jxnTypeIs){
                ytop <- ybottom + d[i][[as.character(jnxTypeI)]]
                if(ytop == ybottom) next
                lines(# ensure at least a line is visible in case rectangle disappears in wide plot
                    x = c(d$x[i], d$x[i]),
                    y = c(ybottom, ytop),
                    col = af_junctions$getIndexColor(jnxTypeI),
                    lwd = 1.5
                )
                rect(
                    xleft   = d$x[i] - 0.5,
                    xright  = d$x[i] + 0.5,
                    ybottom = ybottom,
                    ytop    = ytop,
                    col     = af_junctions$getIndexColor(jnxTypeI),
                    border  = NA
                )
                ybottom <- ytop
            }
        }

        # random offset probabilities
        if(xlim[1] >= -21){
            maxRandomMH <- length(randomOffsetProbs) - 1
            nJxnsInRandomWindow <- sum(dd$y[dd$x %between% c(-maxRandomMH, 0)])
            lines(
                x = 0:-maxRandomMH,
                y = randomOffsetProbs * nJxnsInRandomWindow,
                col = CONSTANTS$plotlyColors$black,
                lwd = 2
            )            
        }

        # event counts by type
        d <- colSums(d[, .SD, .SDcols = jxnTypeIs_])
        inc <- ylim[2] / 10
        for(j in 1:length(jxnTypeIs_)) {
            mtext(
                text = paste(
                    af_junctions$getJxnTypeLabels(jxnTypeIs[j]),
                    sprintf("%5d", d[jxnTypeIs_[j]]),
                    sep = " "
                ),
                side = 4,
                at   = ylim[2] - inc * j,
                line = 0.5,
                las  = 1,
                adj  = 0,
                col  = af_junctions$getIndexColor(jxnTypeIs[j]),
                cex  = 0.8,
                family = "monospace"
            )
        }
    }
    plot$finishPng(layout)
}
offsetPlot <- function(id, xlim, lwd, v){
    settings <- offsetPlotSettings
    settings$Offset_Plot$Min_Offset$value <- xlim[[1]]
    settings$Offset_Plot$Max_Offset$value <- xlim[[2]]
    plot <- mdiInteractivePlotBoxServer(
        id,
        click    = FALSE,
        brush    = TRUE,
        create   = function(...) createOffsetPlot(..., plot, v), # a function or reactive that creates the plot as a png file using settings and helpers
        points   = FALSE, # set to TRUE to expose relevant plot options
        lines    = TRUE,
        legend   = FALSE,
        margins  = TRUE,
        title    = TRUE,
        data     = FALSE,
        size     = "m",
        settings = settings,
        defaults = list(
            Plot_Frame = list(
                Width_Inches  = 5,
                Height_Inches = 3,
                Font_Size     = 12,
                Bottom_Margin = 4.1,
                Left_Margin   = 4.1,
                Top_Margin    = 1.9,
                Right_Margin  = 6.25
            ),
            Points_and_Lines = list(
                Line_Width = lwd
            )
        )
    )
    observeEvent(plot$plot$brush(), {
        coord <- plot$plot$brush()$coord
        minX <- min(coord$x1, coord$x2)
        maxX <- max(coord$x1, coord$x2)
        jxnI <- if(plot$plot$brush()$keys$shift) junctions_filtered()[
            !(alnOffset %between% c(minX, maxX)),
            jxnI
        ] else junctions_filtered()[
            alnOffset %between% c(minX, maxX),
            jxnI
        ]
        offsetPlotSelection(list(x1 = minX, x2 = maxX))
        selectedJxnI$offsetPlot <- jxnI
    })
    plot
}
offsetPlotWide <- offsetPlot(
    id   = "offsetPlotWide",
    xlim = c(-210, 210),
    lwd  = 1.5, 
    v    = c(seq(-1000, 1000, 100), -10, 10)
)
offsetPlotNarrow <- offsetPlot(
    id   = "offsetPlotNarrow",
    xlim = c(-21, 21),
    lwd  = 3, 
    v    = c(seq(-100, 100, 10), c(-5, -2, -1))
)

#----------------------------------------------------------------------
# junction selection summary
#----------------------------------------------------------------------
observeEvent(input$clearSelections, {
    startSpinner(session, message = "clearing selections")
    selectedJxnI$sizePlot   <- c()
    selectedJxnI$offsetPlot <- c()
    sizePlotSelection(NULL)
    offsetPlotSelection(NULL)
    stopSpinner(session)
})
junctionsTable <- bufferedTableServer(
    "selectionSummaryTable",
    id,
    input,
    count_junctions_by_both,
    selection = 'single',
    selectionFn = function(selectedRows) NULL,
    options = list()
)

#----------------------------------------------------------------------
# junctionsTable table, filtered by both sizePlot and offsetPlot selections
#----------------------------------------------------------------------
junctionsTableData <- reactive({
    jxns <- junctions_selected_both()
    x <- jxns[, .SD, .SDcols = names(af_bgzColumns_display[[jxnFileType]])]
    setnames(x, af_bgzColumns_display[[jxnFileType]])
    x[, ":="(
        chrom1 = af_getChromNames(jxns$sourceId, chrom1),
        chrom2 = af_getChromNames(jxns$sourceId, chrom2),
        strand1 = ifelse(strand1 == 0, "+", "-"),
        strand2 = ifelse(strand2 == 0, "+", "-"),
        jxnType = af_junctions$getJxnTypeLabels(jxnType),
        expected = ifelse(expected == 0, TRUE, FALSE)
    )]
    stopSpinner(session)
    x
})
junctionsTable <- bufferedTableServer(
    "junctionsTable",
    id,
    input,
    junctionsTableData,
    selection = 'single',
    selectionFn = function(selectedRows) NULL,
    options = list()
)
selectedJunction <- reactive({
    i <- junctionsTable$rows_selected()
    jxn <- if(isTruthy(i)) junctions_selected_both()[i] else NULL
    updateNumericInput(session, "currReadI1", value = 1)
    fjxn <- af_getFullJunction(jxn) # get full junction with SEQ, QUAL and CIGAR
    jxn[, ":="(
        quals  = fjxn$quals,
        seqs   = fjxn$seqs,
        cigars = fjxn$cigars
    )]
})

#----------------------------------------------------------------------
# single selected junction (from table) expansion
#----------------------------------------------------------------------
ALNI0_1 <- 0
ALNI0_2 <- 1
ALNI1_1 <- 1
ALNI1_2 <- 2
TOP_STRAND0 <- 0
BOTTOM_STRAND0 <- 1
INT_TRUE <- 1
commify <- function(x) format(x, big.mark = ",", scientific = FALSE)
junctionExpandData <- reactive({
    jxn <- selectedJunction()
    req(jxn)
    summary <- jxn[, .(
        sample_ = paste0("SV in ", sample),
        jxnType_ = paste0(
            " is a ", commify(svSize), " bp ",
            if(isValidated) "validated " else "",
            if(isExpected) "expected " else "unexpected ",
            if(isIntergenome) " inter-genomic " else "",
            af_junctions$getJxnTypeNames(jxnType)
        ),
        nodes_ = paste0(
            " from ",
            paste0(af_getChromNames(sourceId, chromIndex1_1), ":", commify(refPos1_1), if(strandIndex0_1 == TOP_STRAND0) "+" else "-" ),
            " to ",
            paste0(af_getChromNames(sourceId, chromIndex1_2), ":", commify(refPos1_2), if(strandIndex0_2 == TOP_STRAND0) "+" else "-" )
        ),
        observed_ = paste0(
            " observed ", nObserved, " times out of coverage ",
            bkptCoverage_1, "/", bkptCoverage_2
        ), 
        jxn_ = paste0(
            " with a ", 
            if(alnOffset < 0) {
                paste0(-alnOffset, " bp junction microhomology = ", jxnBases)
            } else if (alnOffset > 0) {
                paste0(alnOffset, " bp junction insertion = ", jxnBases)
            } else {
                "blunt junction"
            }
        ),
        target1_ = paste0(
            " where ",
            " node 1 is ", commify(targetDist1), " bp from center of target ", target1, " in gene ", genes1
        ),
        target2_ = paste0(
            " and ", 
            " node 2 is ", commify(targetDist2), " bp from center of target ", target2, " in gene ", genes2
        ),
        quality_ = paste0(
            " with min flanking MAPQ = ",  mapQ, 
            " and max base divergence = ", deTag
        ),
        qName = strsplit(qNames, ",")[[1]], # one or more source reads per output row
        seq   = strsplit(seqs,   ",")[[1]], # seq and qual match read aln1, i.e., the first of the two conjoined cigar strings
        qual  = strsplit(quals,  ",")[[1]],
        cigars = strsplit(cigars, ",")[[1]], # alignments were not flipped upstream, always come to us in original read order
        orientation = as.integer(strsplit(orientations, ",")[[1]]), # whether the alignments were/should be flipped for this read to match the junction canonical node order
        strands = .(c(strandIndex0_1, strandIndex0_2)), # junction strands _after_ junction was flipped to canonical order
        alnOffset = alnOffset,
        jxnColor = ifelse(isIntergenome, CONSTANTS$plotlyColors$purple, af_junctions$getIndexColor(jxnType))
    )]
    stopSpinner(session)
    summary
})
selectedReadPath <- reactive({
    jxn <- selectedJunction()
    # rps <- readPaths()
    qName_ <- junctionExpandData()$qName[currReadI1()]
    # rps[sourceId == jxn$sourceId & qName == qName_]
    af_getReadPath(jxn$sourceId, qName_)
})
readPathExpandData <- reactive({
    rp <- selectedReadPath()
    x <- rp[, .(
        qName       = qName,
        qLen        = qLen,
        chroms      = .(strsplit(chroms, ",")[[1]]),
        pos1s       = .(as.integer(strsplit(pos1s, ",")[[1]])),
        strand0s    = .(as.integer(strsplit(strand0s, ",")[[1]])),
        blockNs     = .(as.integer(strsplit(blockNs, ",")[[1]])),
        nRefBases   = .(as.integer(strsplit(nRefBases, ",")[[1]])),
        qryStart0s  = .(as.integer(strsplit(qryStart0s, ",")[[1]])),
        qryEnd1s    = .(as.integer(strsplit(qryEnd1s, ",")[[1]])),
        cigars      = .(strsplit(cigars, ",")[[1]]),
        qual        = .(sapply(strsplit(qual, "")[[1]], function(x) as.integer(charToRaw(x)) - 33L))
    )]
    stopSpinner(session)
    req(nrow(x) > 0)
    x
})
output$nExpansionReads <- renderText({
    paste("out of ", nrow(junctionExpandData()), "reads")
})
currReadI1 <- reactive({
    i1 <- input$currReadI1
    req(i1 >= 1 && i1 <= nrow(junctionExpandData()))
    i1
})  
observeEvent(input$prevReadI1, {
    i1 <- currReadI1()
    req(i1 > 1)
    updateNumericInput(session, "currReadI1", value = i1 - 1)
})
observeEvent(input$nextReadI1, {
    i1 <- currReadI1()
    req(i1 < nrow(junctionExpandData()))
    updateNumericInput(session, "currReadI1", value = i1 + 1)
})
# parseExtendeQName <- function(qName){
#     qName <- strsplit(qName, ":")[[1]]
#     n <- length(qName) - N_QNAME_EXTENSIONS
#     extensions <- qName[(n+1):length(qName)]
#     qName <- qName[1:n]
#     list(
#         qName = qName,
#         extensions = extensions
#     )
#     #AV241004:Aviti2_Run_073_12_19_24:2414452353:2:21201:3910:1360 - 0:0:0:0:0:0:0:1
# }
output$expandedReadSummary <- renderText({
    junctionExpandData()$qName[currReadI1()]
    # qName <- junctionExpandData()$qName[currReadI1()]
    # x <- parseExtendeQName(qName)
    # paste(paste(x$qName, collapse = ":"), " - ", paste(x$extensions, collapse = ":"))
})
parseCigar <- function(cigar, strand0_jxn, wasFlipped, alnI0 = NULL) {
    ops <- data.table(
        n  = as.numeric(unlist(regmatches(cigar, gregexpr('\\d+', cigar)))),
        op =            unlist(regmatches(cigar, gregexpr('\\D',  cigar)))
    )
    n_read_bases_total <- ops[op %in% c("S", "M", "I"), sum(n)]
    n_read_bases_aln   <- ops[op %in% c("M", "I"),      sum(n)]
    n_ref_bases_aln    <- ops[op %in% c("M", "D"),      sum(n)]
    strand0_read <- if(wasFlipped == INT_TRUE) 1 - strand0_jxn else strand0_jxn
    if(strand0_read == BOTTOM_STRAND0) ops <- ops[.N:1]
    list(
        n_read_bases_total = n_read_bases_total,
        n_read_bases_aln   = n_read_bases_aln,
        n_ref_bases_aln    = n_ref_bases_aln,
        ops = ops,
        strand0_jxn  = strand0_jxn,
        strand0_read = strand0_read,
        breakpointV = if(is.null(alnI0)) NULL else if(alnI0 == ALNI0_1){
            n_read_bases_aln + (if(ops[1, op == "S"]) ops[1, n] else 0) + 0.5
        } else {
            n_read_bases_total - n_read_bases_aln - (if(ops[nrow(ops), op == "S"]) ops[nrow(ops), n] else 0) + 1 - 0.5
        }
    )
}
cigars_jxn <- reactive({
    jxn <- junctionExpandData()
    readI1 <- currReadI1()
    cigars <- strsplit(jxn[readI1, cigars], "::")[[1]] # aln1 and aln2 in read order
    wasFlipped <- jxn[readI1, orientation]
    cigars <- lapply(ALNI0_1:ALNI0_2, function(alnI0) {
        strandI0 <- if(wasFlipped == INT_TRUE) 1 - alnI0 else alnI0
        parseCigar( 
            cigars[alnI0 + 1], 
            jxn[readI1, strands[[1]][strandI0 + 1]],
            wasFlipped,
            alnI0
        )
    })
    cigars
})
cigars_rp <- reactive({
    rp <- readPathExpandData()
    lapply(1:length(rp$cigars[[1]]), function(alnI1) {
        parseCigar(
            rp$cigars[[1]][alnI1], 
            rp$strand0s[[1]][alnI1],
            FALSE
        )
    })
})

output$readQualPlot <- renderPlot({
    rp <- readPathExpandData()
    par(mar = c(0.1, 4.1, 0.1, 0.1))
    plot(
        x = 1:rp$qLen,
        y = pmin(rp$qual[[1]], 50), 
        xlim = c(0.5, rp$qLen + 0.5),
        ylim = c(0, 50),
        type = "p", 
        pch = 19,
        cex = 0.25, 
        col = "blue", 
        ylab = "Qual",
        xaxs = "i",
        xaxt = "n"
    )
    abline(h = seq(0, 50, 10), col = "grey80")
    abline(v = c(rp$qryStart0s[[1]] + 1 - 0.5, rp$qryEnd1s[[1]] + 0.5), col = "black")
}, height = dpi * 0.75, res = dpi)

# plot the read to reference alignment in read coordinates
output$refAlnPlot <- renderPlot({
    rp <- readPathExpandData()
    jxn <- junctionExpandData()[currReadI1()]
    cigars_rp <- cigars_rp()
    cigars_jxn <- cigars_jxn()
    ymax <- max(sapply(1:length(cigars_rp), function(alnI1) cigars_rp[[alnI1]]$n_ref_bases_aln))
    par(mar = c(4.1, 4.1, 0.1, 0.1))
    plot(
        NA, NA,
        xlim = c(0.5, rp$qLen + 0.5),
        ylim = c(0.5, ymax + 0.5), 
        ylab = "Ref. Position",
        xlab = "Read Position",
        xaxs = "i",
        yaxs = "i"
    )
    abline(v = c(rp$qryStart0s[[1]] + 1 - 0.5, rp$qryEnd1s[[1]] + 0.5), col = "black")
    for(alnI1 in ALNI1_1:ALNI1_2){
        abline(v = cigars_jxn[[alnI1]]$breakpointV, col = jxn$jxnColor)
    }
    for (alnI1 in 1:length(cigars_rp)){
        cg <- cigars_rp[[alnI1]]
        readOffset0 <- 0
        if(cg$strand0_jxn == TOP_STRAND0) {
            refOffset0 <- 0
            refMultiplier <- 1
        } else {
            refOffset0 <- cg$n_ref_bases_aln 
            refMultiplier <- -1
        }
        for (rowI1 in 1:nrow(cg$ops)) {
            cgg <- cg$ops[rowI1]
            if(cgg$op == "M") lines( # thus, lines plot from start1 - 0.5 to end1 + 0.5
                x = c(readOffset0 + 1 - 0.5, readOffset0 + cgg$n + 0.5),
                y = c(refOffset0  + 1 - 0.5, refOffset0  + cgg$n * refMultiplier + 0.5),
                col = CONSTANTS$plotlyColors$black,
                lwd = 2
            ) else if(cgg$op == "I") lines(
                x = c(readOffset0 + 1 - 0.5, readOffset0 + cgg$n + 0.5),
                y = c(refOffset0  + 1 - 0.5, refOffset0 + 1 - 0.5),
                col = CONSTANTS$plotlyColors$red,
                lwd = 2
            ) else if(cgg$op == "D") lines(
                x = c(readOffset0 + 1 - 0.5, readOffset0 + 1 - 0.5),
                y = c(refOffset0  + 1 - 0.5, refOffset0  + cgg$n * refMultiplier + 0.5),
                col = CONSTANTS$plotlyColors$blue,
                lwd = 2
            )
            if(cgg$op %in% c("M", "I", "S")) readOffset0 <- readOffset0 + cgg$n
            if(cgg$op %in% c("M", "D"))      refOffset0  <- refOffset0  + cgg$n * refMultiplier
        }
    }
}, height = dpi * 3, res = dpi)

# provide a more human readable summary of the junction
output$expandedJunctionSummary <- renderUI({
    summary <- junctionExpandData()[1] # relevant values are the same in all summary rows
    req(summary)
    fluidRow(
        style = "margin: 10px; font-size: 1.2em;",
        lapply(c("sample_", "jxnType_", "nodes_", "observed_", 
                    "jxn_", "target1_", "target2_","quality_"), function(x) {
            tags$p(summary[[x]], style = "margin: 0px; padding: 0px;")
        })
    )
})

# print the parts of the selected read sequence for copy paste
rc_bases <- list(A = "T", C = "G", G = "C", T = "A", N = "N")
rc <- function(SEQ){
    paste(rev(unlist(rc_bases[ strsplit(SEQ, '')[[1]] ])), collapse = "")
}
output$expandedSequences <- renderUI({
    summary <- junctionExpandData()[currReadI1()] 
    cigars <- cigars_jxn()
    seq <- strsplit(
        if(cigars[[ALNI1_1]]$strand0_read == BOTTOM_STRAND0) rc(summary$seq) else summary$seq,
        ""
    )[[1]]
    nBases <- length(seq)
    subSeq <- function(subSeq, label){
        label <- paste0(label, " (", commify(length(subSeq)), " bp)")
        tagList(
            tags$p(label, style = "margin: 0px; padding: 0px; margin-top: 10px; font-weight: bold;"),
            tags$p(paste0(subSeq, collapse = ""), style = "margin: 0px; padding: 0px;")
        )
    }
    nClip5_1 <- if(cigars[[ALNI1_1]]$ops[ 1, op == "S"]) cigars[[ALNI1_1]]$ops[ 1, n] else 0
    nClip5_2 <- if(cigars[[ALNI1_2]]$ops[ 1, op == "S"]) cigars[[ALNI1_2]]$ops[ 1, n] else 0
    nClip3 <- if(cigars[[ALNI1_2]]$ops[.N, op == "S"]) cigars[[ALNI1_2]]$ops[.N, n] else 0
    nBases_1 <- cigars[[ALNI1_1]]$n_read_bases_aln
    fluidRow(
        style = "margin: 10px; font-size: 1.2em;",
        tagList(
            if(nClip5_1 > 0) subSeq(seq[1:nClip5_1], "5' clip") else NULL,
            subSeq(seq[nClip5_1 + (1:nBases_1)], "alignment 1"),
            if(summary$alnOffset < 0) {
                subSeq(seq[nClip5_1 + ((nBases_1 + summary$alnOffset + 1):nBases_1)], "junction microhomology")
            } else if(summary$alnOffset > 0) {
                subSeq(seq[nClip5_1 + nBases_1 + 1:summary$alnOffset], "junction insertion")
            } else NULL,
            subSeq(seq[nClip5_2 + (1:cigars[[ALNI1_2]]$n_read_bases_aln)], "alignment 2"),
            if(nClip3 > 0) subSeq(seq[(nBases - nClip3 + 1):nBases], "3' clip") else NULL,
            subSeq(seq, "entire read")
        )
    )
})

# plot the location of all breakpoints on the two chromosomes for this junction
chromExpandData <- reactive({
    jxn <- selectedJunction()
    req(jxn)
    summary <- jxn[, .(
        sourceId      = sourceId,
        chromIndex1_1 = chromIndex1_1,
        refPos1_1     = refPos1_1,
        chromIndex1_2 = chromIndex1_2,
        refPos1_2     = refPos1_2,
        chrom_1       = af_getChromNames(sourceId, chromIndex1_1),
        chrom_2       = af_getChromNames(sourceId, chromIndex1_2),
        chromSize_1   = af_getChromSize(sourceId, chromIndex1_1),
        chromSize_2   = af_getChromSize(sourceId, chromIndex1_2)
    )]
    stopSpinner(session)
    summary
})
createChromPlot <- function(jxn, chromIndex1, refPos1, chrom_, chromSize){

    # parse the plot width and bin size
    maxChromSize <- jxn[, max(chromSize_1, chromSize_2)]
    nChromBins <- 1000
    chromBinSize <- maxChromSize / nChromBins

    # collect all filtered junction breakpoints on the chromosome
    jxns <- junctions_filtered()
    # jxns <- junctions_selected_both()
    bkpts <- data.table(bin = floor(c(
        jxns[chromIndex1_1 == chromIndex1, refPos1_1],
        jxns[chromIndex1_2 == chromIndex1, refPos1_2]
    ) / chromBinSize) * chromBinSize + chromBinSize / 2)
    bkpts <- bkpts[, .N, keyby = .(bin)]

    # plot the number of junctions in each bin
    par(mar = c(4.1, 4.1, 0.1, 0.1))
    plot(
        NA, 
        xlim = c(0.5, maxChromSize + 0.5),
        ylim = c(0, max(bkpts$N) * 1.05),
        ylab = "Count",
        xaxs = "i",
        xlab = paste(chrom_, "position (bp)"),
    )

    # denote exclusion regions as grey areas
    excls <- af_getExcludedRegions(jxn$sourceId)[chrom == chrom_]
    rect(
        xleft   = excls$start0 + 1,
        xright  = excls$end1,
        ybottom = 0,
        ytop    = max(bkpts$N) * 1.025,
        col     = "grey80",
        border  = NA
    )

    # highlight this junction breakpoint
    jxn <- junctionExpandData()[currReadI1()]
    abline(v = refPos1, col = jxn$jxnColor, lwd = 1.5)

    # draw the chrom bin points
    points(
        x = bkpts$bin,
        y = bkpts$N, 
        pch = 19,
        cex = 0.25, 
        col = "blue"
    )
}
output$chrom1Plot <- renderPlot({
    jxn <- chromExpandData()[currReadI1()] 
    createChromPlot(jxn, jxn$chromIndex1_1, jxn$refPos1_1, jxn$chrom_1, jxn$chromSize_1)
}, height = dpi * 1.5, res = dpi)
output$chrom2Plot <- renderPlot({
    jxn <- chromExpandData()[currReadI1()] 
    createChromPlot(jxn, jxn$chromIndex1_2, jxn$refPos1_2, jxn$chrom_2, jxn$chromSize_2)
}, height = dpi * 1.5, res = dpi)

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
        offsetPlotWide$settings$replace(bm$outcomes$offsetPlotWideSettings)
        offsetPlotNarrow$settings$replace(bm$outcomes$offsetPlotNarrowSettings)
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
        sizePlotSettings = sizePlot$settings$all_(),
        offsetPlotWideSettings = offsetPlotWide$settings$all_(),
        offsetPlotNarrowSettings = offsetPlotNarrow$settings$all_()
    ) }),
    settingsObject = settings,
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------

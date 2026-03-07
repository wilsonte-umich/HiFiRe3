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
sourceId <- dataSourceTableServer("dataSource", selection = "single") 
junctions_all <- reactive({ # all junction from data source, prior to filtering
    sourceId <- sourceId()
    req(sourceId)
    startSpinner(session, message = "loading junctions")
    hf3_getJunctions_all_source(sourceId)[
        jxn_type == hf3_junctions$typeToBits["translocation"] | 
        sv_size > 0 # suppress 0bp inversions
    ]
})
junctions_filtered <- reactive({ # all junctions, filtered by top-level settings
    junctions_all <- junctions_all()
    req(junctions_all)
    startSpinner(session, message = "loading filtered junctions")
    junctions_all %>% hf3_applyJunctionFilters(settings, input)
})

#----------------------------------------------------------------------
# count SV junctions
#----------------------------------------------------------------------
count_junctions_by_type <- function(jxns){
    agg <- jxns[, .N, keyby = .(jxn_type, is_intergenomic)]
    agg[,  ":="(
        jxnTypeLabel = hf3_junctions$getTypeLabelsFromBits(jxn_type, is_intergenomic),
        color        = hf3_junctions$getColorsFromBits(    jxn_type, is_intergenomic),
        labelOffset  = hf3_junctions$bitsToIndex[jxn_type] + is_intergenomic
    )]
    agg[order(jxn_type)]
}
count_junctions_by_class <- function(jxns){
    agg <- jxns[, .N, keyby = .(jxnClass)]
    agg[, jxnClassName := hf3_junctions$getJxnClassLabels(jxnClass)]
    agg[, color := hf3_junctionClasses$getClassColor(jxnClass)]
    agg
}
count_junctions_by_both <- reactive({
    jxns <- junctions_selected_both()
    req(nrow(jxns) > 0)
    agg <- jxns[, .N, keyby = .(jxn_type, is_intergenomic, jxnClass)]
    agg[, jxnTypeLabel := hf3_junctions$getTypeLabelsFromBits(jxn_type, is_intergenomic)]
    agg[, jxnClassName := hf3_junctionClasses$getJxnClassLabels(jxnClass)]
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
    # samples <- unique(jxns$sample_names)
    jxns[, stratumLabel := hf3_junctionStrata$getJxnStratumLabels(jxnStratum)]
    stratumLabels <- hf3_junctionStrata$indexToStratumLabel
    strata <- 1:length(stratumLabels)
    names(strata) <- stratumLabels
    jitterAmount <- 0.4
    typeJitterAmount <- jitterAmount / 3
    jxns <- jxns[, {
        stratumI <- strata[stratumLabel]
        jxnTypeSubRowOffset <- hf3_junctions$getOffsetsFromBits(jxn_type, is_intergenomic)
        workingLog10Size <- ifelse(
            jxn_type == hf3_junctions$typeToBits["translocation"], 
            jitter(rep(9.25, .N), amount = 0.75), 
            sv_size %>% log10
        )
        pointColor <- hf3_junctions$getColorsFromBits(jxn_type, is_intergenomic)
        .(
            stratumI = stratumI,
            jxnI = jxnI, # for interactive selection and mapping to source data.table
            x = workingLog10Size,
            y_ = stratumI + jxnTypeSubRowOffset, 
            color = pointColor,
            jxn_type = jxn_type,
            is_intergenomic = is_intergenomic,
            hasExcluded = is_excluded_1 == 1 | is_excluded_2 == 1
        )
    }, by = .(stratumLabel)]
    jxns[, y := jitter(y_, amount = typeJitterAmount)] # spread the subrow points through a y-axis jitter
    list(
        # samples       = samples,
        jxns          = jxns,
        stratumLabels = stratumLabels,
        strata        = strata,
        nStrata       = length(strata),
        jitterAmount  = 1/3 + typeJitterAmount,
        labelHeight   = 2 * (1/3 + typeJitterAmount) / 5
    )
})
createSizePlot <- function(settings, plot) {
    req(!input$suspendPlotting)
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
        # title = d$samples,
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
    abline(h = 1:(d$nStrata - 1) + 0.5, col = "black")

    # add individual junction points, with downsampling for speed
    plot_jxns <- if(input$maxPlottedPoints > 0 && nrow(d$jxns) > input$maxPlottedPoints){
        d$jxns[sample(1:nrow(d$jxns), input$maxPlottedPoints)]
    } else {
        d$jxns
    }
    colors <- ifelse(
        plot_jxns$hasExcluded,
        CONSTANTS$plotlyColors$black,
        plot_jxns$color
    )
    plot$addPoints(
        x   = plot_jxns$x,
        y   = plot_jxns$y,
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
        for (jxn_type_ in agg$jxn_type) {
            x <- agg[jxn_type == jxn_type_]
            mtext(
                text = x[, paste(jxnTypeLabel, sprintf("%5d", x$N), sep = " ")],
                side = 4,
                at   = i + d$jitterAmount - d$labelHeight * x$labelOffset,
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
            Height_Inches = 4, # was 5
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
    req(!input$suspendPlotting)

    # get the data
    d <- junctions_filtered()[!is.na(offset)] # may be missing when SEQ was dropped for known SVs
    if(length(selectedJxnI$sizePlot) > 0) {
        d <- d[jxnI %in% selectedJxnI$sizePlot]
    }
    if(nrow(d) == 0) {
        stopSpinner(session)
        req(FALSE)
    }
    startSpinner(session, message = "plotting offsets")

    # aggregate by offset
    d <- d[, .(y = .N), keyby = .(
        x = offset, 
        jxn_type = jxn_type, 
        is_intergenomic = is_intergenomic,
        jxnTypeI = hf3_junctions$bitsToIndex[jxn_type] + is_intergenomic)
    ]
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
        yaxs = "i"
        # ,
        # title = sizePlotJunctions()$samples
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
                color <- hf3_junctions$typeToColor[[hf3_junctions$indexToType[jnxTypeI]]]
                lines(# ensure at least a line is visible in case rectangle disappears in wide plot
                    x = c(d$x[i], d$x[i]),
                    y = c(ybottom, ytop),
                    col = color,
                    lwd = 1.5
                )
                rect(
                    xleft   = d$x[i] - 0.5,
                    xright  = d$x[i] + 0.5,
                    ybottom = ybottom,
                    ytop    = ytop,
                    col     = color,
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
            color <- hf3_junctions$typeToColor[[hf3_junctions$indexToType[jxnTypeIs[j]]]]
            mtext(
                text = paste(
                    hf3_junctions$indexToTypeLabel[jxnTypeIs[j]],
                    sprintf("%5d", d[jxnTypeIs_[j]]),
                    sep = " "
                ),
                side = 4,
                at   = ylim[2] - inc * j,
                line = 0.5,
                las  = 1,
                adj  = 0,
                col  = color,
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
            !(offset %between% c(minX, maxX)),
            jxnI
        ] else junctions_filtered()[
            offset %between% c(minX, maxX),
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
summaryTableData <- reactive({
    req(!input$suspendPlotting)
    count_junctions_by_both()
})
junctionsTable <- bufferedTableServer(
    "selectionSummaryTable",
    id,
    input,
    summaryTableData,
    selection = 'single',
    selectionFn = function(selectedRows) NULL,
    options = list()
)

#----------------------------------------------------------------------
# junctionsTable table, filtered by both sizePlot and offsetPlot selections
#----------------------------------------------------------------------
junctionsTableData <- reactive({
    req(!input$suspendPlotting)
    jxns <- junctions_selected_both()
    if (input$maxTableRows > 0 && nrow(jxns) > input$maxTableRows) {
        jxns <- jxns[sample(1:nrow(jxns), input$maxTableRows)]
    }
    x <- jxns[, .SD, .SDcols = names(hf3_bgzColumns_display[[jxnFileType]])]
    setnames(x, hf3_bgzColumns_display[[jxnFileType]])
    x[, ":="(
        chrom1   = hf3_getChromNames(jxns$sourceId, chrom1),
        chrom2   = hf3_getChromNames(jxns$sourceId, chrom2),
        strand1  = ifelse(strand1 == 0, "+", "-"),
        strand2  = ifelse(strand2 == 0, "+", "-"),
        jxnType  = hf3_junctions$getTypeLabelsFromBits(jxnType, interGen)
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
    jxn
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
        sample_ = paste0("SV in ", sample_names),
        jxnType_ = paste0(
            " is a ", commify(sv_size), " bp ",
            if(is_validated) "validated " else "",
            if(is_intergenomic) " inter-genomic " else "",
            hf3_junctions$getTypesFromBits(jxn_type, is_intergenomic)
        ),
        nodes_ = paste0(
            " from ",
            paste0(hf3_getChromNames(sourceId, chrom_index1_1), ":", commify(ref_pos1_1), if(strand_index0_1 == TOP_STRAND0) "+" else "-" ),
            " to ",
            paste0(hf3_getChromNames(sourceId, chrom_index1_2), ":", commify(ref_pos1_2), if(strand_index0_2 == TOP_STRAND0) "+" else "-" )
        ),
        observed_ = paste0(
            " observed ", n_instances_dedup, " times out of coverage ",
            bkpt_coverage_1, "/", bkpt_coverage_2
        ), 
        jxn_ = paste0(
            " with a ", 
            if(offset < 0) {
                paste0(-offset, " bp junction microhomology")
            } else if (offset > 0) {
                paste0(offset, " bp junction insertion")
            } else {
                "blunt junction"
            }
        ),
        target1_ = paste0(
            " where ",
            " node 1 is ", commify(target_dist_1), " bp from center of target ", target_1, " in gene ", genes_1
        ),
        target2_ = paste0(
            " and ", 
            " node 2 is ", commify(target_dist_2), " bp from center of target ", target_2, " in gene ", genes_2
        ),
        quality_ = paste0(
            " with min flanking MAPQ = ",  max_min_mapq, 
            " and max base divergence = ", min_max_divergence
        ),
        qName = strsplit(sub(",", "", q_names), ",")[[1]], # one or more source reads per output row
        isDuplex = strsplit(sub(",", "", is_duplexes), ",")[[1]],
        aln5I0 = as.integer(strsplit(sub(",", "", aln5_is), ",")[[1]]), # alignment index (0-based) of 5' side of junction within read
        orientation = as.integer(strsplit(sub(",", "", jxn_orientations), ",")[[1]]), # whether the alignments were/should be flipped for this read to match the junction canonical node order
        strands = .(c(strand_index0_1, strand_index0_2)), # junction strands _after_ junction was flipped to canonical order
        offset = offset,
        jxnColor = hf3_junctions$getColorsFromBits(jxn_type, is_intergenomic),
        aln_failure_flag,
        jxn_failure_flag
    )]
    stopSpinner(session)
    summary
})
selectedReadPath <- reactive({
    jxn <- selectedJunction()
    qName <- junctionExpandData()$qName[currReadI1()]
    hf3_getReadPath(jxn$sourceId, qName)
})
readPathExpandData <- reactive({
    rp <- selectedReadPath()
    if(nrow(rp) == 0){
        stopSpinner(session)
        req(FALSE)
    }
    rp[, .(
        qName      = qname,
        qLen       = read_len,
        chroms     = .(strsplit(chroms, ",")[[1]]),
        pos1s      = .(as.integer(strsplit(pos1s, ",")[[1]])),
        strand0s   = .(as.integer(strsplit(strand0s, ",")[[1]])),
        nRefBases  = .(as.integer(strsplit(n_ref_bases, ",")[[1]])),
        qryStart0s = .(as.integer(strsplit(qry_start0s, ",")[[1]])),
        qryEnd1s   = .(as.integer(strsplit(qry_end1s, ",")[[1]])),
        blockNs    = .(as.integer(strsplit(block_ns, ",")[[1]])),
        mapqs      = .(sapply(strsplit(mapqs, ",")[[1]], function(x) as.integer(x))),
        divergences= .(sapply(strsplit(divergences, ",")[[1]], function(x) as.numeric(x))),
        alnFailureFlags = .(as.integer(strsplit(aln_failure_flags, ",")[[1]])),
        jxnFailureFlags = .(as.integer(strsplit(jxn_failure_flags, ",")[[1]])),
        cigars     = .(strsplit(cigars, ",")[[1]]),
        seqStrand0 = seq_strand0,
        seq        = seq,
        qual       = .(sapply(strsplit(qual, "")[[1]], function(x) as.integer(charToRaw(x)) - 33L))
    )]
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
output$expandedReadQname <- renderText({
    i <- currReadI1()
    req(i)
    jed <- junctionExpandData()
    duplexStatus <- if(jed$isDuplex[i] == "1") " (duplex)" else ""
    paste(jed$qName[i], duplexStatus)
})
parseCigar <- function(cigar, strand0, alnI0 = NULL) {
    ops <- data.table(
        n  = as.numeric(unlist(regmatches(cigar, gregexpr('\\d+', cigar)))),
        op =            unlist(regmatches(cigar, gregexpr('\\D',  cigar)))
    )
    n_read_bases_total <- ops[op %in% c("S", "M", "I"), sum(n)]
    n_read_bases_aln   <- ops[op %in% c("M", "I"),      sum(n)]
    n_ref_bases_aln    <- ops[op %in% c("M", "D"),      sum(n)]
    if(strand0 == BOTTOM_STRAND0) ops <- ops[.N:1]
    list(
        n_read_bases_total = n_read_bases_total,
        n_read_bases_aln   = n_read_bases_aln,
        n_ref_bases_aln    = n_ref_bases_aln,
        ops = ops,
        strand0  = strand0,
        breakpointV = if(is.null(alnI0)) NULL else if(alnI0 == ALNI0_1){
            n_read_bases_aln + (if(ops[1, op == "S"]) ops[1, n] else 0) + 0.5
        } else {
            n_read_bases_total - n_read_bases_aln - (if(ops[nrow(ops), op == "S"]) ops[nrow(ops), n] else 0) + 1 - 0.5
        }
    )
}
cigars_jxn <- reactive({
    jxn <- junctionExpandData()[currReadI1()]
    rp <- readPathExpandData()
    lapply(jxn$aln5I0 + 0:1, function(alnI0) {
        parseCigar(
            rp$cigars[[1]][alnI0 + 1], 
            rp$strand0s[[1]][alnI0 + 1],
            alnI0 - jxn$aln5I0
        )
    })
})
cigars_rp <- reactive({
    rp <- readPathExpandData()
    lapply(1:length(rp$cigars[[1]]), function(alnI1) {
        parseCigar(
            rp$cigars[[1]][alnI1], 
            rp$strand0s[[1]][alnI1]
        )
    })
})

output$readQualPlot <- renderPlot({
    rp <- readPathExpandData()
    qual <- rp$qual[[1]]
    if (rp$seqStrand0 == BOTTOM_STRAND0) qual <- rev(qual)
    par(mar = c(0.1, 4.1, 0.1, 0.1))
    plot(
        x = 1:rp$qLen,
        y = jitter(qual, a = 2), 
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
        abline(v = cigars_jxn[[alnI1]]$breakpointV, col = jxn$jxnColor, lwd = 2)
    }
    for (alnI1 in 1:length(cigars_rp)){
        cg <- cigars_rp[[alnI1]]
        readOffset0 <- 0
        if(cg$strand0 == TOP_STRAND0) {
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
    tags$div(
        lapply(c("sample_", "jxnType_", "nodes_", "observed_", 
                    "jxn_", "target1_", "target2_","quality_"), function(x) {
            tags$p(summary[[x]], style = "margin: 0px; padding: 0px;")
        })
    )
})
output$expandedJunctionFlagBits <- renderUI({
    summary <- junctionExpandData()[1] # relevant values are the same in all summary rows
    req(summary)
    fluidRow(
        column(
            width = 6,
            tags$p(tags$strong("Alignment Failure Flag")),
            lapply(names(hf3_alnFailureBits), function(flagName) {
                failed <- bitwAnd(summary$aln_failure_flag, hf3_alnFailureBits[[flagName]]) != 0
                color <- if(failed) CONSTANTS$plotlyColors$red else CONSTANTS$plotlyColors$black
                tags$p(style = paste0("margin: 0px; padding: 0px; color: ", color, ";"), paste(
                    flagName, ":", 
                    if(failed) "failed" else "passed/NA"
                ))
            })
        ),
        column(
            width = 6,
            tags$p(tags$strong("Junction Failure Flag")),
            lapply(names(hf3_jxnFailureBits), function(flagName) {
                failed <- bitwAnd(summary$jxn_failure_flag, hf3_jxnFailureBits[[flagName]]) != 0
                color <- if(failed) CONSTANTS$plotlyColors$red else CONSTANTS$plotlyColors$black
                tags$p(style = paste0("margin: 0px; padding: 0px; color: ", color, ";"), paste(
                    flagName, ":", 
                    if(failed) "failed" else "passed/NA"
                ))
            })
        )
    )
})

# print the parts of the selected read sequence for copy paste
rc_bases <- list(A = "T", C = "G", G = "C", T = "A", N = "N")
rc <- function(SEQ){
    paste(rev(unlist(rc_bases[ strsplit(SEQ, '')[[1]] ])), collapse = "")
}
output$expandedSequences <- renderUI({
    summary <- junctionExpandData()[currReadI1()] 
    rp <- readPathExpandData()
    cigars <- cigars_jxn()
    seq <- strsplit(
        if(rp$seqStrand0 == BOTTOM_STRAND0) rc(rp$seq) else rp$seq,
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
    # q1 <- rp$qual[[1]][nClip5_1 + (1:nBases_1)]
    # q2 <- rp$qual[[1]][nClip5_2 + (1:cigars[[ALNI1_2]]$n_read_bases_aln)]
    fluidRow(
        style = "margin: 10px; font-size: 1.2em;",
        tagList(
            if(nClip5_1 > 0) subSeq(seq[1:nClip5_1], "5' clip") else NULL,
            subSeq(seq[nClip5_1 + (1:nBases_1)], "alignment 1"),
            if(summary$offset < 0) {
                subSeq(seq[nClip5_1 + ((nBases_1 + summary$offset + 1):nBases_1)], "junction microhomology")
            } else if(summary$offset > 0) {
                subSeq(seq[nClip5_1 + nBases_1 + 1:summary$offset], "junction insertion")
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
        chromIndex1_1 = chrom_index1_1,
        refPos1_1     = ref_pos1_1,
        chromIndex1_2 = chrom_index1_2,
        refPos1_2     = ref_pos1_2,
        chrom_1       = hf3_getChromNames(sourceId, chrom_index1_1),
        chrom_2       = hf3_getChromNames(sourceId, chrom_index1_2),
        chromSize_1   = hf3_getChromSize( sourceId, chrom_index1_1),
        chromSize_2   = hf3_getChromSize( sourceId, chrom_index1_2)
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
    bkpts <- data.table(bin = floor(c(
        jxns[chrom_index1_1 == chromIndex1, ref_pos1_1],
        jxns[chrom_index1_2 == chromIndex1, ref_pos1_2]
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
    excls <- hf3_getExcludedRegions(jxn$sourceId)[chrom == chrom_]
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
    abline(v = refPos1, col = jxn$jxnColor, lwd = 2)

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
    junctions_filtered = junctions_filtered,
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------

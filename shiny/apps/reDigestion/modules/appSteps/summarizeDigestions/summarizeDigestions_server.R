#----------------------------------------------------------------------
# server components for the summarizeDigestions appStep module
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# BEGIN MODULE SERVER
#----------------------------------------------------------------------
summarizeDigestionsServer <- function(id, options, bookmark, locks) { 
    moduleServer(id, function(input, output, session) {    
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# initialize module
#----------------------------------------------------------------------
module <- 'summarizeDigestions'
appStepDir <- getAppStepDir(module)
options <- setDefaultOptions(options, stepModuleInfo[[module]])
settings <- activateMdiHeaderLinks( # uncomment as needed
    session,
    # url = getDocumentationUrl("path/to/docs/README", domain = "xxx"), # for documentation
    # dir = appStepDir, # for terminal emulator
    envir = environment(), # for R console
    baseDirs = appStepDir, # for code viewer/editor
    settings = id, # for step-level settings
    # immediate = TRUE # plus any other arguments passed to settingsServer()
    size = "m"
)

#----------------------------------------------------------------------
# load the data source and it requested objects
#----------------------------------------------------------------------
sourceId <- dataSourceTableServer("dataSource", selection = "single") 
summaries <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    summariesFile <- getSourceFilePath(sourceId, "summariesFile")
    readRDS(summariesFile)
})
distributions <- reactive({
    sourceId <- sourceId()
    req(sourceId)
    distributionsFile <- getSourceFilePath(sourceId, "distributionsFile")
    readRDS(distributionsFile)
})
enzymesData <- reactive({
    summaries()$siteSummaries[
        CpG_priority >= settings$get("RE_Settings","Min_Priority") & 
        CpG_priority <= settings$get("RE_Settings","Max_Priority")
    ]
})
siteSummariesData <- reactive({
    summaries()$siteSummaries
})

#----------------------------------------------------------------------
# sweet spot size range
#----------------------------------------------------------------------
sweetSpot <- reactive({
    sweetSpotMin_kb <- settings$get("RE_Settings", "Sweet_Spot_Min_kb")
    sweetSpotMax_kb <- sweetSpotMin_kb * 2
    sweetSpotMin <- sweetSpotMin_kb * 1000
    sweetSpotMax <- sweetSpotMin * 2
    list(
        min_kb = sweetSpotMin_kb, 
        max_kb = sweetSpotMax_kb,
        suffix = paste(sweetSpotMin, sweetSpotMax, sep = "_"),
        min = sweetSpotMin, 
        max = sweetSpotMax
    )
})

#----------------------------------------------------------------------
# enzymes table
#----------------------------------------------------------------------
enzymesTableData <- reactive({
    sweetSpot <- req(sweetSpot())
    enzymes <- req(enzymesData())
    nSweet <- unlist(enzymes[, .SD, .SDcols = paste("N", sweetSpot$suffix, sep = "_")])
    fSweet <- unlist(enzymes[, .SD, .SDcols = paste("massFrac", sweetSpot$suffix, sep = "_")])
    enzymes[, .(
        RE = enzyme,
        site = cut_site,
        N50 = N50,
        fSweet = round(fSweet, 3),
        nSweet = nSweet,
        priority = CpG_priority,
        high_fid = high_fidelity,
        weight = effective_length,
        CpG_sens = CpG_sensitive,
        CpG_level = CpG_level,
        cents_unit = cents_unit,
        temp = temperature,
        buffer = buffer,
        star = star_activity
    )]
})
enzymesTable <- bufferedTableServer(
    "enzymesTable",
    id,
    input,
    enzymesTableData,
    selection = 'single',
    selectionFn = function(selectedRows) NULL,
    options = list()
)
selectedEnzyme <- reactive({
    i <- enzymesTable$rows_selected()
    req(i)
    enzymesData()[i]
})
selectedEnzymeDistribution <- reactive({
    selectedEnzyme <- selectedEnzyme()
    req(selectedEnzyme)
    distributions()[[selectedEnzyme$reKey]]
})

#----------------------------------------------------------------------
# short and long read windowed size distribution plots
#----------------------------------------------------------------------
plotDistributions <- function(plot, xScalar, xlabUnit){
    sweetSpot <- req(sweetSpot())
    enzymes <- siteSummariesData()[
        CpG_priority >= settings$get("RE_Settings","Min_Priority") & 
        CpG_priority <= settings$get("RE_Settings","Max_Priority")
    ][order(N50)] 
    nEnzymes <- nrow(enzymes)
    maxY <- nrow(enzymes) + 1
    par(mar = c(4.1,6.1,0,0))
    plot$initializeFrame(
        xlim = c(0, sweetSpot$max * 2) / xScalar,
        ylim = c(0, maxY),
        xlab = paste0("Fragment Size (", xlabUnit, ")"),
        ylab = "",
        yaxt = "n",
        xaxs = "i",
        yaxs = "i"
    )
    axis(2, at = 1:nEnzymes, labels = paste0(enzymes$enzyme, " (", enzymes$effective_length, ")"), tick = TRUE, las = 1, cex = 0.75)
    rect(sweetSpot$min / xScalar, 0, sweetSpot$max / xScalar, maxY, col = "grey90", border = "NA")
    abline(v = c(sweetSpot$min, sweetSpot$max) / xScalar, col = CONSTANTS$plotlyColors$red, lwd = 1)
    for(i in 1:nEnzymes){
        enz <- enzymes[i]
        col <- if(enz$CpG_priority >= 4) "blue3" else "black"
        lines(c(enz$N2.5, enz$N97.5) / xScalar, c(i, i), col = col)
        rect(enz$N5 / xScalar, i - 0.25, enz$N95 / xScalar, i + 0.25, col = "white", border = col)
        lines(c(enz$N10, enz$N10) / xScalar, c(i - 0.25, i + 0.25), col = col)
        lines(c(enz$N90, enz$N90) / xScalar, c(i - 0.25, i + 0.25), col = col)
        lines(c(enz$N50, enz$N50) / xScalar, c(i - 0.4, i + 0.4), col = col)
        points(c(enz$maxLen) / xScalar, c(i), pch = 19, cex = 0.5, col = col)
    }
    rect(0, 0, sweetSpot$max * 2, maxY)
}
distributionPlot <- staticPlotBoxServer(
    "distributionPlot",
    settings = NULL,
    maxHeight = "600px",
    create = function() plotDistributions(distributionPlot, 1000, "kb")
)

#----------------------------------------------------------------------
# number of sites vs. N50
#----------------------------------------------------------------------
nSitesPlot <- staticPlotBoxServer(
    "nSitesPlot",
    settings = NULL,
    maxHeight = "600px",
    create = function() {
        sweetSpot <- req(sweetSpot())
        enzymes <- siteSummariesData()[
            CpG_priority >= settings$get("RE_Settings","Min_Priority") & 
            CpG_priority <= settings$get("RE_Settings","Max_Priority")
        ][order(N50)]
        par(mar = c(4.1,6.1,0,0))
        scale <- 1000
        unit <- if(scale == 1) "" else if(scale == 1e3) "(K)" else "(M)"
        barplot(
            enzymes$nSites / scale,
            log = "x",
            width = 0.75,
            names.arg = enzymes$enzyme,
            horiz = TRUE,
            las = 1,
            xlab = paste("# of sites", unit),
            cex.names = 0.8,
            space = 0.4
        )
    }
)

#----------------------------------------------------------------------
# single-enzyme distributions
#----------------------------------------------------------------------
binEnzymePlotData <- function(dist, maxLength, fractionCol){
    dist <- dist[length <= maxLength]
    nBins <- settings$get("RE_Settings", "Bins_Per_Plot")
    # binSize <- diff(log10(range(dist$length))) / nBins
    # dist[, 
    #     bin := floor(log10(length) / binSize) * binSize
    # ][, 
    #     .(fraction = sum(.SD[[fractionCol]])), 
    #     keyby = .(bin)
    # ]  
    binSize <- diff(range(dist$length)) / nBins
    dist[, 
        bin := floor(length / binSize) * binSize
    ][, 
        .(fraction = sum(.SD[[fractionCol]])), 
        keyby = .(bin)
    ]    
}
binPlotV <- function(selectedEnzymeSummary, prefix){
    levels <- 50 # c(1, 5, 50, 95, 99)
    names <- paste0(prefix, levels)
    v <- unlist(selectedEnzymeSummary[, .SD, .SDcols = names])
    names(v) <- names
    v
}
plotEnzymeDistribution <- function(plot, distType, prefix){
    xScalar <- 1e3
    xlabUnit <- "kb"
    sweetSpot <- req(sweetSpot())
    maxX <- sweetSpot$max * 2
    selectedEnzyme <- selectedEnzyme()
    selectedEnzymeSummary <- siteSummariesData()[reKey == selectedEnzyme$reKey]
    selectedEnzymeDistribution <- selectedEnzymeDistribution()
    fractionCol <- paste(distType, "fraction", sep = "_")
    nCol <- "N"
    bins <- binEnzymePlotData(selectedEnzymeDistribution, maxX, fractionCol)
    maxY <- bins[, max(fraction)]
    plotMaxY <- maxY * 1.05
    par(mar = c(4.1,6.1,2.1,0))
    plot$initializeFrame(
        xlim = c(0, maxX) / xScalar,
        ylim = c(0, plotMaxY),
        xlab = paste0("Fragment Size (", xlabUnit, ")"),
        ylab = stringr::str_to_title(paste(distType, "fraction", sep = " ")),
        xaxs = "i",
        yaxs = "i",
        title = selectedEnzyme$reKey
    )
    rect(sweetSpot$min / xScalar, 0, sweetSpot$max / xScalar, plotMaxY, col = "grey90", border = "NA")
    v <- binPlotV(selectedEnzymeSummary, prefix) / xScalar
    abline(v = v, col = CONSTANTS$plotlyColors$blue, lwd = 1.25)
    abline(v = c(sweetSpot$min, sweetSpot$max) / xScalar, col = CONSTANTS$plotlyColors$red, lwd = 1)
    text(
        x = v,
        y = rep(maxY, length(v)),
        labels = names(v), 
        col = CONSTANTS$plotlyColors$blue,
        pos = 4
    )
    plot$addPoints(
        x = bins$bin / xScalar,
        # x = bins$bin,
        y = bins$fraction,
        col = CONSTANTS$plotlyColors$blue,
        cex = 0.75
    )
    addSpanLabel <- function(x, I) text(
        x = x,
        y = maxY * 0.8,
        labels = paste(
            paste0(round(sum(selectedEnzymeDistribution[I][[fractionCol]]) * 100, 1), "%"),
            paste0(round(sum(selectedEnzymeDistribution[I][[nCol]]) / 1000, 0), "K"),
            sep = "\n"
        ), 
        col = CONSTANTS$plotlyColors$red,
        pos = 1
    )
    addSpanLabel(
        (sweetSpot$min + sweetSpot$max) / 2 / xScalar, 
        selectedEnzymeDistribution[, length >= sweetSpot$min & length < sweetSpot$max]
    )
    addSpanLabel(
        sweetSpot$min / 2 / xScalar, 
        selectedEnzymeDistribution[, length < sweetSpot$min]
    )
    rect(0, 0, maxX / xScalar, maxY * 1.05, border = "black")
}
molarFractionPlot <- staticPlotBoxServer(
    "molarFractionPlot",
    settings = NULL,
    create = function() plotEnzymeDistribution(molarFractionPlot, "molar", "q")
)
massFractionPlot <- staticPlotBoxServer(
    "massFractionPlot",
    settings = NULL,
    create = function() plotEnzymeDistribution(massFractionPlot, "mass", "N")
)

#----------------------------------------------------------------------
# define bookmarking actions
#----------------------------------------------------------------------
bookmarkObserver <- observe({
    bm <- getModuleBookmark(id, module, bookmark, locks)
    req(bm)
    settings$replace(bm$settings)
    # updateSelectInput(session, "sampleSet-sampleSet", selected = bm$input[['sampleSet-sampleSet']])
    if(!is.null(bm$outcomes)) {
        # outcomes <<- listToReactiveValues(bm$outcomes)
        distributionPlot$settings$replace(bm$outcomes$distributionPlotSettings)
        nSitesPlot$settings$replace(bm$outcomes$nSitesPlotSettings)
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
        distributionPlotSettings = distributionPlot$settings$all_(),
        nSitesPlotSettings = nSitesPlot$settings$all_()
    ) }),
    # isReady = reactive({ getStepReadiness(options$source, ...) }),
    NULL
)

#----------------------------------------------------------------------
# END MODULE SERVER
#----------------------------------------------------------------------
})}
#----------------------------------------------------------------------

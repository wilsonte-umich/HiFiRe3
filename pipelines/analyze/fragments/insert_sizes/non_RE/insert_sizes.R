# action:
#   assess and plot the insert size distribution in a non-RE sequenced library
# input:
#     $INSERT_SIZES_FILE
# outputs:
#     distribution plot

#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("initializing")
suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
    library(parallel)
}))
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'INSERT_SIZES_FILE',
        'READ_LENGTH_TYPE',
        'IS_FIXED_READ_LENGTH',
        'PLOT_PREFIX'
    ),
    integer = c(
        'FIXED_READ_LENGTH',
        'MIN_INSERT_SIZE'
    )
))
env$LONG_READ_BIN_SIZE <- 25
env$MIN_RAW_COVERAGE   <- 0.05 # used for peak finding
#-------------------------------------------------------------------------------------
# set variable derived from environment variables
if(env$READ_LENGTH_TYPE == "short"){
    isShortRead <- TRUE
    peakSizeMultiple <- 3
    xlab <- "Insert Size (bp)"
    xScalar <- 1
    genomeXlab <- "Bin Index (M)"
    genomeXScalar <- 1e6
    stepSize <- 1
    maxFitSize <- 1000
    maxX <- as.integer(env$MAX_INSERT_SIZE)
} else {
    isShortRead <- FALSE
    peakSizeMultiple <- 4
    xlab <- "Insert Size (kb)"
    xScalar <- 1e3
    genomeXlab <- "Bin Index (K)"
    genomeXScalar <- 1e3
    stepSize <- env$LONG_READ_BIN_SIZE
    maxFitSize <- 35000
    maxX <- maxFitSize
}
#-------------------------------------------------------------------------------------
# source R scripts
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
sourceScripts(rUtilDir, c('utilities'))
sourceScripts(file.path(rUtilDir, 'genome'), c('chroms', 'stats'))
sourceScripts(file.path(rUtilDir, 'sequence'), c('general'))
setCanonicalChroms()
nonNuclear   <- c("chrM", "chrEBV")
nonAutosomes <- c("chrX", "chrY", nonNuclear)
nuclearChroms <- canonicalChroms[!(canonicalChroms %in% nonNuclear)]
autosomes <- canonicalChroms[!(canonicalChroms %in% nonAutosomes)]
env$MAX_INSERT_SIZE <- getMaxInsertSize()
#-------------------------------------------------------------------------------------
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#=====================================================================================

#=====================================================================================
# feedback functions
#-------------------------------------------------------------------------------------
insertColor_represented <- "blue"
readLenColor <- "purple2"
insertColor_transparent <- addAlphaToColor(insertColor_represented, 0.1)
maxXMulitiplier <- 5
plotIt <- function(plotName, plotFn, width, ...){
    png(filename = paste(env$PLOT_PREFIX, plotName, "png", sep = "."),
        width = width, height = 3, units = "in", pointsize = 8,
        bg = "white",  res = 600, type = "cairo")
    result <- plotFn(plotName, ...)
    dev.off()
    result
}
plotSizeDistributions <- function(plotName, fracPrefix, ylab, peak = NULL){
    fracColumn_inserts   <- paste(fracPrefix, "inserts",   sep = "_")

    # plot the distributions
    maxY <- max(sizeDists[[fracColumn_inserts]], na.rm = TRUE)
    plot(
        x = NA,
        y = NA,
        xlab = xlab,
        ylab = ylab,
        xlim = c(0, maxX) / xScalar,
        ylim = c(0, maxY),
        main = plotName
    )

    # plot guidance
    abline(v = N50_inserts / xScalar, col = insertColor_represented)
    if(env$IS_FIXED_READ_LENGTH == "TRUE") abline(v = env$FIXED_READ_LENGTH / xScalar, col = readLenColor)
    # abline(h = peak$minFrac, col = insertColor_represented)

    # plot distributions
    points(
        sizeDists$binnedSize / xScalar,
        sizeDists[[fracColumn_inserts]],
        col = insertColor_represented,
        pch = 19,
        cex = 0.5
    )

    # add legend
    legend(
        "topright",
        legend = c("obs. inserts"),
        col    = c(insertColor_represented),
        pch = 19
    )
}
#=====================================================================================

#=====================================================================================
# assess the observed insert size distribution
#-------------------------------------------------------------------------------------
message("loading observed insert sizes")
observedSizes <- fread(env$INSERT_SIZES_FILE)
setnames(observedSizes, c("size", "nInserts"))
observedSizes <- observedSizes[size %between% c(env$MIN_INSERT_SIZE, env$MAX_INSERT_SIZE)]

message("aggregating observed insert size distribution")
observedSizes[, binnedSize := if(isShortRead) size else as.integer(size / env$LONG_READ_BIN_SIZE) * env$LONG_READ_BIN_SIZE]
observedSizes <- observedSizes[, .(nInserts = sum(nInserts)), keyby = .(binnedSize)]
observedSizes[, molarFrac_inserts := nInserts / sum(nInserts)]
observedSizes[, massFrac_inserts  := {
    x <- as.numeric(nInserts) * as.numeric(binnedSize) # prevent integer overflow
    x / sum(x)
}]
#-------------------------------------------------------------------------------------
sizeDists <- observedSizes[, .(binnedSize, molarFrac_inserts, massFrac_inserts)]
sizeDists[is.na(molarFrac_inserts),   molarFrac_inserts   := 0]
sizeDists[is.na(massFrac_inserts),    massFrac_inserts    := 0]
#=====================================================================================

#=====================================================================================
# assess and plot the size distributions
#-------------------------------------------------------------------------------------
message("plotting insert size distributions")
N50_inserts <- N50_weighted(sizeDists$binnedSize, sizeDists$molarFrac_inserts)
# maxX <- env$MAX_INSERT_SIZE # max(N50_inserts * maxXMulitiplier, coverageSizeRange[2] * 1.1)
invisible(plotIt("insert_sizes_molar", plotSizeDistributions, width = 3, "molarFrac", "Molar Fraction"))
invisible(plotIt("insert_sizes_mass",  plotSizeDistributions, width = 3, "massFrac",  "Mass Fraction"))
#=====================================================================================

#=====================================================================================
# write output for fragment coverage assessment 
#-------------------------------------------------------------------------------------
# representation <- sizeDists[
#     representation >= env$MIN_REPRESENTATION & # set the fragments sizes allowed to have coverage tracking
#     binnedSize %between% coverageSizeRange, 
#     .(binnedSize, representation) # the correction factor to apply to observed fragment counts
# ]
# eventual output must include 
#   binSize <- 1 || env$LONG_READ_BIN_SIZE
#   representation
#=====================================================================================

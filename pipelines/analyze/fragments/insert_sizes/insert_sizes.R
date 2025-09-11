# action:
#   assess and plot the insert size distribution in the genome and sequenced ligFree library, where:
#       fragment representation in library = observed insert size frequency / expected fragment size frequency
# input:
#     $INSERT_SIZES_FILE
# outputs:
#     $INSERT_SIZE_REPRESENTATION

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
        'SITE_COLLECTION_SITES',
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
env$MIN_REPRESENTATION <- 0.25 # used for restricted usable insert size range
#-------------------------------------------------------------------------------------
# set variable derived from environment variables
if(env$READ_LENGTH_TYPE == "short"){
    isShortRead <- TRUE
    peakSizeMultiple <- 3
    xlab <- "Fragment Size (bp)"
    xScalar <- 1
    genomeXlab <- "Bin Index (M)"
    genomeXScalar <- 1e6
    stepSize <- 1
    maxFitSize <- 1000
} else {
    isShortRead <- FALSE
    peakSizeMultiple <- 4
    xlab <- "Fragment Size (kb)"
    xScalar <- 1e3
    genomeXlab <- "Bin Index (K)"
    genomeXScalar <- 1e3
    stepSize <- env$LONG_READ_BIN_SIZE
    maxFitSize <- 35000
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
fragmentColor <- "grey80"
insertColor_represented <- "blue"
insertColor_not_represented <- "cornflowerblue"
readLenColor <- "purple2"
insertColor_transparent <- addAlphaToColor(insertColor_represented, 0.1)
maxPlottedFragments <- 25000
minFragmentInserts <- 2
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
    fracColumn_fragments <- paste(fracPrefix, "fragments", sep = "_")
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
    abline(v = N50_fragments / xScalar, col = fragmentColor)
    abline(v = N50_inserts   / xScalar, col = insertColor_represented)
    if(env$IS_FIXED_READ_LENGTH == "TRUE") abline(v = env$FIXED_READ_LENGTH / xScalar, col = readLenColor)
    # abline(h = peak$minFrac, col = insertColor_represented)

    # plot distributions
    points(
        sizeDists$binnedSize / xScalar,
        sizeDists[[fracColumn_fragments]],
        col = fragmentColor,
        pch = 19,
        cex = 0.5
    )
    points(
        sizeDists$binnedSize / xScalar,
        sizeDists[[fracColumn_inserts]],
        col = ifelse(sizeDists$representation >= env$MIN_REPRESENTATION, insertColor_represented, insertColor_not_represented),
        pch = 19,
        cex = 0.5
    )

    # add legend
    legend(
        "topright",
        legend = c("obs. inserts", "ref. fragments"),
        col    = c(insertColor_represented, fragmentColor),
        pch = 19
    )
}
plotRepresentation <- function(plotName){
    # plot the distributions
    maxY <- max(sizeDists$representation, na.rm = TRUE)
    plot(
        x = NA,
        y = NA,
        xlab = xlab,
        ylab = "Relative Representation",
        xlim = c(0, maxX) / xScalar,
        ylim = c(0, maxY),
        main = plotName
    )

    # plot guidance
    abline(v = N50_fragments / xScalar, col = fragmentColor)
    abline(h = 1, col = fragmentColor)    
    abline(v = N50_inserts   / xScalar, col = insertColor_represented)
    if(env$IS_FIXED_READ_LENGTH == "TRUE") abline(v = env$FIXED_READ_LENGTH / xScalar, col = readLenColor)
    abline(v = coverageSizeRange / xScalar, col = "red3")    
    abline(h = env$MIN_REPRESENTATION, col = insertColor_represented)

    # plot distributions
    points(
        sizeDists$binnedSize / xScalar,
        sizeDists$representation,
        col = ifelse(sizeDists$representation >= env$MIN_REPRESENTATION, insertColor_represented, insertColor_not_represented),
        pch = 19,
        cex = 0.5
    )
    tryCatch({
        lines(
            sizeDists_fitLow$binnedSize / xScalar,
            predict(fitLow, sizeDists_fitLow),
            col = "darkorange"
        )
        lines(
            sizeDists_fitHigh$binnedSize / xScalar,
            predict(fitHigh, sizeDists_fitHigh),
            col = "darkorange"
        )
    }, error = function(e) NULL)
}
#=====================================================================================

#=====================================================================================
# load and analyze the RE sites on the reference genome (sufficient for assessing the fragment size distribution)
#-------------------------------------------------------------------------------------
message("loading in silico RE site positions")
sites <- fread(env$SITE_COLLECTION_SITES)
sites <- sites[bitwAnd(haplotypes, 1) == 1]
#=====================================================================================

#=====================================================================================
# assess the expected fragment size distribution
#-------------------------------------------------------------------------------------
message("calculating resulting RE fragment lengths")
fragments <- do.call(rbind, mclapply(nuclearChroms[nuclearChroms %in% unique(sites$chrom)], function(chrom_){
    sites[chrom == chrom_, .(
        chrom = chrom_,
        size = diff(sitePos1)
    )]
}, mc.cores = env$N_CPU))
fragments <- fragments[size %between% c(env$MIN_INSERT_SIZE, env$MAX_INSERT_SIZE)]

message("determining size distribution of autosomal RE fragments")
fragments[, binnedSize := if(isShortRead) size else as.integer(size / env$LONG_READ_BIN_SIZE) * env$LONG_READ_BIN_SIZE]
expectedSizes <- fragments[chrom %in% autosomes, .(nFragments = .N), keyby = .(binnedSize)]
expectedSizes[, molarFrac_fragments := nFragments / sum(nFragments)]
expectedSizes[, massFrac_fragments  := (nFragments * binnedSize) / sum(nFragments * binnedSize)]
#=====================================================================================

#=====================================================================================
# assess the expected fragment size distribution
#-------------------------------------------------------------------------------------
message("loading observed insert sizes")
observedSizes <- fread(env$INSERT_SIZES_FILE)
setnames(observedSizes, c("size", "nInserts"))
observedSizes <- observedSizes[size %between% c(env$MIN_INSERT_SIZE, env$MAX_INSERT_SIZE)]

message("aggregating observed insert size distribution")
observedSizes[, binnedSize := if(isShortRead) size else as.integer(size / env$LONG_READ_BIN_SIZE) * env$LONG_READ_BIN_SIZE]
observedSizes <- observedSizes[, .(nInserts = sum(nInserts)), keyby = .(binnedSize)]
observedSizes[, molarFrac_inserts := nInserts / sum(nInserts)]
observedSizes[, massFrac_inserts  := (nInserts * binnedSize) / sum(nInserts * binnedSize)]
#=====================================================================================

#=====================================================================================
# compare the expected fragment and observed insert size distributions
#-------------------------------------------------------------------------------------
message("calculating fragment size representation in library inserts, observed / expected")
sizeDists <- merge(
    expectedSizes[, .(binnedSize, molarFrac_fragments, massFrac_fragments)],
    observedSizes[, .(binnedSize, molarFrac_inserts,   massFrac_inserts)],
    by = "binnedSize",
    all = TRUE
)
sizeDists[is.na(molarFrac_fragments), molarFrac_fragments := 0]
sizeDists[is.na(massFrac_fragments),  massFrac_fragments  := 0]
sizeDists[is.na(molarFrac_inserts),   molarFrac_inserts   := 0]
sizeDists[is.na(massFrac_inserts),    massFrac_inserts    := 0]
sizeDists[molarFrac_fragments > 0, representation := molarFrac_inserts / molarFrac_fragments]

message("setting size boundaries usable for rendering coverage maps")
N50_fragments <- N50_weighted(sizeDists$binnedSize, sizeDists$molarFrac_fragments)
N50_inserts   <- N50_weighted(sizeDists$binnedSize, sizeDists$molarFrac_inserts)
sizeDists_usable <- sizeDists[!is.na(representation) & representation > 0 & binnedSize < maxFitSize]
sizeDists_usable[, ":="(
    binnedSize2 = binnedSize ** 2,
    binnedSize3 = binnedSize ** 3
)]
sizeDists_fitLow <- sizeDists_usable[binnedSize < N50_inserts]
sizeDists_fitHigh <- sizeDists_usable[binnedSize > N50_inserts]

# catch error on low complexity cases, e.g.nCATS
fitSucceeded <- TRUE
coverageSizeRange <- tryCatch({
    fitLow <- lm(representation ~ binnedSize + binnedSize2, data = sizeDists_fitLow)
    yLow <- predict(fitLow, sizeDists_fitLow)
    fitHigh <- lm(representation ~ binnedSize + binnedSize2, data = sizeDists_fitHigh)
    yHigh <- predict(fitHigh, sizeDists_fitHigh)
    c(
        sizeDists_fitLow[ min(which(yLow  > env$MIN_REPRESENTATION)), binnedSize],
        sizeDists_fitHigh[max(which(yHigh > env$MIN_REPRESENTATION)), binnedSize]
    )
}, error = function(e){
    fitSucceeded <- FALSE
    if(sum(sizeDists$binnedSize > env$MIN_REPRESENTATION) > 2){
        c(
            sizeDists[binnedSize > env$MIN_REPRESENTATION, min(binnedSize)],
            sizeDists[binnedSize > env$MIN_REPRESENTATION, max(binnedSize)]
        )
    } else {
        c(
            sizeDists[, min(binnedSize)],
            sizeDists[, max(binnedSize)]
        )
    }
})
#=====================================================================================

#=====================================================================================
# assess and plot the size distributions
#-------------------------------------------------------------------------------------
message("plotting fragment vs. insert size distributions")
maxX <- max(N50_inserts * maxXMulitiplier, coverageSizeRange[2] * 1.1)
invisible(plotIt("insert_sizes_molar", plotSizeDistributions, width = 3, "molarFrac", "Molar Fraction"))
invisible(plotIt("insert_sizes_mass",  plotSizeDistributions, width = 3, "massFrac",  "Mass Fraction"))
invisible(plotIt("representation",     plotRepresentation,    width = 3))
#=====================================================================================

#=====================================================================================
# write output for fragment coverage assessment 
#-------------------------------------------------------------------------------------
representation <- sizeDists[
    representation >= env$MIN_REPRESENTATION & # set the fragments sizes allowed to have coverage tracking
    binnedSize %between% coverageSizeRange, 
    .(binnedSize, representation) # the correction factor to apply to observed fragment counts
]
# eventual must include 
#   binSize <- 1 || env$LONG_READ_BIN_SIZE
#   representation
#=====================================================================================

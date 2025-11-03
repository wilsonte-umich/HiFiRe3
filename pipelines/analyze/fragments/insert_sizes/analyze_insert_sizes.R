# action:
#     create plots of insert size distributions before and after filtering and RE site projection
#     establish automated thresholds for allowed insert sizes working with the hint from --min-selected-size
# output: 
#     
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
        'UNFILTERED_INSERT_SIZES_FILE',
        'FILTERED_INSERT_SIZES_FILE',
        'FILTERED_STEM_LENGTHS_FILE',
        'PLOT_PREFIX',
        'DATA_NAME',
        'TARGETS_BED'
    ),
    integer = c(
        'MIN_SELECTED_SIZE',
        'N_CPU'
    ),
    double = c(
        'SELECTED_SIZE_CV'
    )
))
#-------------------------------------------------------------------------------------
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#-------------------------------------------------------------------------------------
env$ENOUGH_READS <- 50
#=====================================================================================

#=====================================================================================
# derive minAllowedSize and maxAllowedSize from MIN_SELECTED_SIZE
#-------------------------------------------------------------------------------------
minAllowedSize <- env$MIN_SELECTED_SIZE / (1 + env$SELECTED_SIZE_CV)
maxAllowedSize <- minAllowedSize * 2
#=====================================================================================

#=====================================================================================
# load the insert sizes
#-------------------------------------------------------------------------------------
insertSizeCols  <- c("type","sizeBin","nObserved")
insertSizeTypes <- list(
    unfiltered      = "black",   # pre-alignment insert sizes
    nonSV.actual    = "green4",  # observed insert sizes from non-SV reads passing all quality and chimeric filters
    nonSV.projected = "blue2",   # non-SV actual projected as needed to nearest RE site
    SV.chimeric     = "purple3", # insert sizes from all reads rejected as chimeric, e.g., end(s) match RE site, jxn has adapter
    SV.passed       = "red3",    # insert sizes from SV reads passing all quality and chimeric filters
    #----------------------------------
    intra.chimeric  = "blue2",
    intra.passed    = "cyan3",
    inter.chimeric  = "red3",
    inter.passed    = "orange2"
)

message("loading insert sizes")
insertSizeCounts <- rbind(
    fread(env$UNFILTERED_INSERT_SIZES_FILE), # all pre-alignment reads, includes on- and off-target ("unfiltered")
    fread(env$FILTERED_INSERT_SIZES_FILE)    # filtered upstream to on-target only (insertSizeTypes other than "unfiltered")
)
setnames(insertSizeCounts, insertSizeCols)
if(env$TARGETS_BED != "NA") insertSizeCounts <- insertSizeCounts[type != "unfiltered" | sizeBin > 1500] # don't tally unfiltered off-target fragments
insertSizeCounts <- dcast(
    insertSizeCounts,
    sizeBin ~ type,
    value.var = "nObserved",
    fill = 0
)[order(sizeBin)]
availableInsertSizeTypeNames <- names(insertSizeTypes[names(insertSizeTypes) %in% names(insertSizeCounts)])
sumInsertSizeCounts <- sapply(availableInsertSizeTypeNames, function(x) sum(insertSizeCounts[[x]]))
names(sumInsertSizeCounts) <- availableInsertSizeTypeNames
hasEnoughReads <- sumInsertSizeCounts >= env$ENOUGH_READS
availableInsertSizeTypeNames <- availableInsertSizeTypeNames[hasEnoughReads]
insertSizeFreqs_molar <- insertSizeCounts[, lapply(.SD, function(x) x / sum(x)), .SDcols = availableInsertSizeTypeNames]
insertSizeFreqs_mass  <- insertSizeCounts[, lapply(.SD, function(x) {
    x <- as.numeric(x)
    x * sizeBin / sum(x * sizeBin)
}), .SDcols = availableInsertSizeTypeNames]
#=====================================================================================

#=====================================================================================
# insert size plotting functions
#-------------------------------------------------------------------------------------
compressInteger <- Vectorize(function(x){
    if (x < 1e3) {
        as.character(x)
    } else if (x < 1e6) {
        paste0(signif(x / 1e3, 3), "K")
    } else if (x < 1e9) {
        paste0(signif(x / 1e6, 3), "M")
    } else {
        paste0(signif(x / 1e9, 3), "B")
    }
})
plotInsertSizes <- function(level, type, insertSizeFreqs, insertSizeTypeNames, x2xFreqs = NULL){
    pngFile <- paste(env$PLOT_PREFIX, level, type, "png", sep = ".")
    png(filename = pngFile,
        width = 3, height = 3, units = "in", pointsize = 8,
        bg = "white", res = 600, type = "cairo")
    maxX <- minAllowedSize * 4
    insertSizeFreqs <- insertSizeFreqs[, ..insertSizeTypeNames]
    maxY <- max(unlist(insertSizeFreqs))
    plot(
        x = NA,
        y = NA,
        xlab = "Insert Size (kb)",
        ylab = "Frequency",
        xlim = c(0, maxX) / 1000,
        ylim = c(0, maxY) * 1.05,
        main = paste(env$DATA_NAME, level, type)
    )
    abline(h = 0)
    abline(v = 1:100, col = "grey80")
    abline(v = c(minAllowedSize, maxAllowedSize) / 1000, col = "black")
    for(type in insertSizeTypeNames){
        points(
            x = insertSizeCounts$sizeBin / 1000,
            y = insertSizeFreqs[[type]],
            col = insertSizeTypes[[type]],
            pch = 19,
            cex = 0.75
        )
    }
    if(!is.null(x2xFreqs)){
        lines(
            x = x2xSizes$sizeBin / 1000,
            y = x2xFreqs$e2e,
            col = insertSizeTypes[["trans.intra.chimeric"]],
            lwd = 2
        )
        lines(
            x = x2xSizes$sizeBin / 1000,
            y = x2xFreqs$m2m,
            col = insertSizeTypes[["trans.intra.passed"]],
            lwd = 2
        )
    }
    legend(
        "topright", 
        legend = paste0(
            sub("trans.", "", insertSizeTypeNames), 
            " (", compressInteger(sumInsertSizeCounts[insertSizeTypeNames]), ")"
        ), 
        col = unlist(insertSizeTypes[insertSizeTypeNames]), 
        pch = 19, 
        cex = 0.8
    )
    dev.off()
}
#=====================================================================================

#=====================================================================================
# make insert size distribution plot with overlaid types unfiltered..projected
#-------------------------------------------------------------------------------------
message("plotting insert sizes")
insertSizeTypeNames <- names(insertSizeTypes)[1:5]
insertSizeTypeNames <- insertSizeTypeNames[insertSizeTypeNames %in% availableInsertSizeTypeNames]
sink <- plotInsertSizes("insert_sizes", "molar", insertSizeFreqs_molar, insertSizeTypeNames)
sink <- plotInsertSizes("insert_sizes", "mass",  insertSizeFreqs_mass,  insertSizeTypeNames)
#=====================================================================================

#=====================================================================================
# establish expected distributions for post-size-selection artifacts:
#     end-to-end, where two size-selected fragments are ligated together
#     middle-to-middle, two pieces of broken size-selected fragments are ligated together
# working from the observed projected insert size distribution that best defines the size of true fragments
#-------------------------------------------------------------------------------------
message("modeling artifact insert size distributions")
minPieceSize <- 100
sizePlotBinSize <- 250
nSamples <- 5e4
message("  sampling projected insert sizes")
x2xSizes <- insertSizeCounts[sizeBin > minPieceSize, .(sizeBin, nonSV.projected)]
I <- 1:nrow(x2xSizes)
sampleI_1 <- sample(I, nSamples, replace = TRUE, prob = x2xSizes$nonSV.projected)
sampleI_2 <- sample(I, nSamples, replace = TRUE, prob = x2xSizes$nonSV.projected)
message("  modeling end-to-end")
e2eSizes <- x2xSizes[sampleI_1, sizeBin] + x2xSizes[sampleI_2, sizeBin]
message("  modeling middle-to-middle")
m2mSizes <- sapply(1:nSamples, function(j){
    size1 <- x2xSizes[sampleI_1[j], sample(minPieceSize:(sizeBin-minPieceSize), 1)] # the size of a piece taken from a broken sizeBin fragment
    size2 <- x2xSizes[sampleI_2[j], sample(minPieceSize:(sizeBin-minPieceSize), 1)]
    as.integer(floor((size1 + size2) / sizePlotBinSize) * sizePlotBinSize)
})
x2xSizes <- merge(
    x2xSizes,
    data.table(sizeBin = e2eSizes)[, .(e2e = .N), keyby = sizeBin],
    by = "sizeBin",
    all.x = TRUE
)
x2xSizes <- merge(
    x2xSizes,
    data.table(sizeBin = m2mSizes)[, .(m2m = .N), keyby = sizeBin],
    by = "sizeBin",
    all.x = TRUE
)
x2xSizes[is.na(e2e), e2e := 0]
x2xSizes[is.na(m2m), m2m := 0]
x2xFreqs <- x2xSizes[, lapply(.SD, function(x) x / sum(x)), .SDcols = c("e2e", "m2m")]
#=====================================================================================

#=====================================================================================
# make insert size distribution plot for exploring (intergenomic) translocation 1N to 2N logic
#-------------------------------------------------------------------------------------
message("plotting translocation insert sizes")
insertSizeTypeNames <- c("nonSV.actual", names(insertSizeTypes)[6:9])
insertSizeTypeNames <- insertSizeTypeNames[insertSizeTypeNames %in% availableInsertSizeTypeNames]
sink <- plotInsertSizes("trans_insert_sizes", "molar", insertSizeFreqs_molar, insertSizeTypeNames, x2xFreqs)
sink <- plotInsertSizes("trans_insert_sizes", "mass",  insertSizeFreqs_mass,  insertSizeTypeNames)
#=====================================================================================

#=====================================================================================
# make min stem length distribution plot for exploring (intergenomic) translocation <1N logic
#-------------------------------------------------------------------------------------
message("loading stem lengths")
stemLengths <- fread(env$FILTERED_STEM_LENGTHS_FILE)
setnames(stemLengths, insertSizeCols)
stemLengths <- dcast(
    stemLengths,
    sizeBin ~ type,
    value.var = "nObserved",
    fill = 0
)[order(sizeBin)]
availableStemLengthTypeNames <- names(insertSizeTypes[names(insertSizeTypes) %in% names(stemLengths)])
sumStemLengthCounts <- sapply(availableStemLengthTypeNames, function(x) sum(stemLengths[[x]]))
names(sumStemLengthCounts) <- availableStemLengthTypeNames
hasEnoughReads <- sumStemLengthCounts >= env$ENOUGH_READS
availableStemLengthTypeNames <- availableStemLengthTypeNames[hasEnoughReads]
if(length(availableStemLengthTypeNames) > 0){
    stemLengthFreqs <- stemLengths[, lapply(.SD, function(x) x / sum(x)), .SDcols = availableStemLengthTypeNames]

    message("plotting translocation stem lengths")
    plotType <- "trans_stem_lengths"
    pngFile <- paste(env$PLOT_PREFIX, plotType, "png", sep = ".")
    png(filename = pngFile,
        width = 3, height = 3, units = "in", pointsize = 8,
        bg = "white", res = 600, type = "cairo")
    maxX <- minAllowedSize * 4
    maxY <- max(unlist(stemLengthFreqs))
    plot(
        x = NA,
        y = NA,
        xlab = "Min Stem Length (kb)",
        ylab = "Frequency",
        xlim = c(0, maxX) / 1000,
        ylim = c(0, maxY) * 1.05,
        main = paste(env$DATA_NAME, plotType)
    )
    abline(h = 0)
    abline(v = 1:100, col = "grey80")
    abline(v = c(minAllowedSize, maxAllowedSize) / 1000, col = "black")
    for(type in availableStemLengthTypeNames){
        points(
            x = stemLengths$sizeBin / 1000,
            y = stemLengthFreqs[[type]],
            col = insertSizeTypes[[type]],
            pch = 19,
            cex = 0.75
        )
    }
    legend(
        "topright", 
        legend = paste0(
            availableStemLengthTypeNames, 
            " (", compressInteger(sumStemLengthCounts[availableStemLengthTypeNames]), ")"
        ), 
        col = unlist(insertSizeTypes[availableStemLengthTypeNames]), 
        pch = 19, 
        cex = 0.8
    )
    dev.off()
}
#=====================================================================================

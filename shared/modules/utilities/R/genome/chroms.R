#----------------------------------------------------------
# genome manipulations
#----------------------------------------------------------
canonicalChroms <- character()
chromIndex      <- list() # 1-reference chromosome indices, e.g., chr3 => 3 ...
revChromIndex   <- list() # ... and vice versa

# parsers for restricting work to properly ordered canonical chromosomes
setCanonicalChroms <- function(){ 
    GENOME_CHROMS <- strsplit(env$GENOME_CHROMS, "\\s+")[[1]]
    USE_ALL_CHROMS <- !(env$USE_ALL_CHROMS == "" || env$USE_ALL_CHROMS == "0")
    IS_COMPOSITE_GENOME <- !(env$IS_COMPOSITE_GENOME == "" || env$IS_COMPOSITE_GENOME == "0")
    if(!is.null(env$GENOME_CHROMS) &&
       (USE_ALL_CHROMS || IS_COMPOSITE_GENOME)){
        canonicalChroms <<- GENOME_CHROMS
    } else {
        addChrom <- function(chr){
            chrom <- paste0("chr", chr)
            if(chrom %in% GENOME_CHROMS) canonicalChroms <<- c(canonicalChroms, chrom)
        } 
        sapply(c(1:90,'X','Y','M','EBV'), addChrom)
    }
    # if(is.null(env$GENOME_CHROMS) || env$USE_ALL_CHROMS == "" || env$USE_ALL_CHROMS == "0"){
    #     addChrom <- function(chr){
    #         chrom <- paste0("chr", chr)
    #         if(chrom %in% GENOME_CHROMS) canonicalChroms <<- c(canonicalChroms, chrom)
    #     } 
    #     sapply(c(1:90,'X','Y','M','EBV'), addChrom)
    # } else {
    #     canonicalChroms <<- GENOME_CHROMS
    # }
    maxI <- length(canonicalChroms)
    for(i in 1:maxI){
        chrom <- canonicalChroms[i]
        chromIndex[[chrom]] <<- i # special handling of unmapped reads
        revChromIndex[[i]] <<- chrom
    }
    chromIndex[['*']] <<- 99 # special handling of unmapped reads
    revChromIndex[[99]] <<- '*'
}

# construct a map of all concatenated, canonical chromosomes from a genome fai file
# must be called after setCanonicalChroms
# requires library(bit64)
loadChromSizes <- function(windowSize = 1000){ 
    message("loading chrom sizes")
    index_file <- paste(env$GENOME_FASTA, 'fai', sep = ".")
    x <- fread(
        index_file, 
        header = FALSE, 
        sep = "\t", 
        stringsAsFactors = FALSE,
        col.names  = c('chrom',     'nChromBases', 'offset',  'lineSeqLen', 'lineLen'), 
        colClasses = c('character', 'integer64',   'numeric', 'integer',    'integer')
    )[chrom %in% names(chromIndex)] # canonical chroms only
    x[, chromIndex := unlist(chromIndex[chrom])]
    x <- x[order(chromIndex)]
    x[, nChromWindows := as.integer64(floor((nChromBases - 1) / windowSize) + 1)]
    getRunningCum <- function(n){
        x <- cumsum(n)
        c(as.integer64(0), x[1:(length(x) - 1)])        
    }
    x[, ":="(
        nBasesBefore   = getRunningCum(nChromBases),
        nWindowsBefore = getRunningCum(nChromWindows)
    )]
    x[, ":="(
        nBasesThrough   = nBasesBefore   + nChromBases,
        nWindowsThrough = nWindowsBefore + nChromWindows
    )]
    x <- x[, .SD, .SDcols = c(
        "chromIndex","chrom",
        "nChromBases",  "nBasesBefore",  "nBasesThrough",
        "nChromWindows","nWindowsBefore","nWindowsThrough"
    )]
    setkey(x, chromIndex)
    x
}
expandChromWindows <- function(chromSizes){
    message("expanding chrom windows for faster lookup")
    chromSizes[, .(
        chrom = chrom,
        windowIndex = 1:as.integer(nChromWindows) - 1
    ), keyby = .(chromIndex)]
}

# nodes are codified into an integer64 for streamlined comparison
# this function expands integer64 nodes out to chrom/strand/pos, returned as a data.table
parseSignedNodes <- function(chromSizes, nodes, side = NULL, canonical = FALSE) {
    genomeIs <- abs(nodes) # 1-referenced, per initialize_windows.pl
    chromIs <- Vectorize(function(i) which(chromSizes$nBasesThrough >= i)[1])(genomeIs) # sapply does not work with integer64!
    dt <- chromSizes[chromIs][, .(
        chrom = chrom, 
        chromIndex = chromIndex,
        refPos = as.integer(genomeIs - nBasesBefore),
        strand = ifelse(nodes > 0, "+", "-")
    )]    
    if(canonical) setnames(dt, c("cChrom","cChromIndex","cRefPos","cStrand"))
    if(!is.null(side)) setnames(dt, paste0(names(dt), side))
    dt
}
# this function collapes chrom/strand/pos to integer64 nodes, returned as a vector of integer64
getSignedNode <- function(chromSizes, chromIndex, coordinate, strand, add = 0) {
    chromSizes[chromIndex, nBasesBefore] + coordinate + add * ifelse(strand == "+", 1, -1)
}

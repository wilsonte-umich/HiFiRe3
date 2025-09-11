#----------------------------------------------------------
# manipulations related to genome target regions
#----------------------------------------------------------

# load the target regions as a table of all positions by region
targetRegions <- NULL
loadTargetRegions <- function(){
    if(is.null(env$TARGETS_BED) || 
       env$TARGETS_BED == "" || 
       env$TARGETS_BED == "NA" || 
       env$TARGETS_BED == "null") return(NULL)
    message("initializing target regions")
    if(is.null(env$REGION_PADDING)) env$REGION_PADDING <- 1
    env$REGION_PADDING <- as.integer(env$REGION_PADDING)
    targets <- fread(env$TARGETS_BED, sep = "\t", stringsAsFactors = FALSE)
    if(ncol(targets) < 4) targets[[4]] <- "NA"
    targets <- targets[, 1:4]
    names(targets) <- c('chrom', 'start', 'end', 'regionName')
    targets[, ':='(
        chromI = unlist(chromIndex[chrom]),
        start = start + 1
    )]
    x <- rbind(
        targets[, .(
            chromI,
            pos = (start - env$REGION_PADDING):(start - 1),
            type = "A"
        ), by = regionName],
        targets[, .( 
            chromI, 
            pos = start:end,
            type = "T"
        ), by = regionName],
        targets[, .(
            chromI, 
            pos = (end + 1):(end + env$REGION_PADDING),
            type = "A"
        ), by = regionName]
    )
    x <- dcast(x, chromI + pos ~ regionName, value.var = 'type', fill = "-")
    regionNames <- names(x)[-(1:2)]
    x[, I := .I]
    setcolorder(x, c('chromI', 'pos', 'I', regionNames))
    targetRegions <<- list(
        N = targets[, .N],
        names = regionNames,
        positions = x
    )
    targets[, ':='(
        start = start - 1, # convert back to bed for printing
        paddedStartI = getTargetRegionsI(chromI, start - env$REGION_PADDING),
        startI       = getTargetRegionsI(chromI, start),
        centerI      = getTargetRegionsI(chromI, start + floor((end - start) / 2)),
        endI         = getTargetRegionsI(chromI, end),
        paddedEndI   = getTargetRegionsI(chromI, end + env$REGION_PADDING),
        size         = end - start + 1
    )][order(chromI, start, end)]
    targetRegions$bed <<- targets
}

# get the index of a position in an ordered, collapsed concatenation of all padded regions
getTargetRegionI <- function(chromI_, pos_){
    if(is.null(targetRegions)) return(0L)
    x <- targetRegions$positions[chromI == chromI_ & pos == pos_]
    if(nrow(x) != 1) return(0L)
    x[, I]
}
getTargetRegionsI <- Vectorize(getTargetRegionI)

# get a vector of the names of all padded target regions a position corresponds to
# typically will be one region name, unless user settings lead to overlapping padded regions
getTargetRegionName <- function(chromI_, pos_){
    if(is.null(targetRegions)) return("*")
    x <- targetRegions$positions[chromI == chromI_ & pos == pos_]
    if(nrow(x) != 1) return("*")
    x <- unlist(x[, -(1:3)])
    x <- x[x != "-"]
    if(length(x) == 0) return("*")
    names(x)
}
getTargetRegionsName <- Vectorize(getTargetRegionName)

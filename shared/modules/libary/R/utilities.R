# load RE filtering sites previously established by prepare/locate
loadFilteringSites <- function(){
    message(paste("loading RE filtering sites"))
    x <- fread(
        env$RE_SITES_TABLE,
        sep = "\t",
        quote = ""
    )
    x[, ":="(
        inSilico = as.logical(inSilico),
        rNameIndex = unlist(chromIndex[chrom])
    )]
    x
}

# load extract input data as edge nodes for RE site matching
loadObservedEndpoints <- function() {
    message(paste("loading observed endpoint nodes"))
    x <- fread(
        paste(env$EXTRACT_SVS_PREFIX, "observed_endpoints.txt.gz", sep = "."),
        col.names  = nodeCols,
        colClasses = nodeColClasses,
        sep = "\t",
        quote = ""
    )
    x[, ":="( # adjust all values to use as positionOffsets array indices
        strand     = strand + 1,    # 1/2 for +/-, top/bottom
        isJunction = isJunction + 1 # 1/2 logical, i.e., 1 = alignment, 2 = junction
    )]
    x
}

# load extract input data as edge tables
loadSvEdges <- function() {
    message(paste("loading edges from SV-containing reads"))
    x <- fread(
        paste(env$EXTRACT_SVS_PREFIX, "edges.sv.txt.gz", sep = "."),
        col.names  = edgeCols,
        colClasses = edgeColClasses,
        sep = "\t",
        quote = ""
    )
    x[, ":="(
        isFoldback = as.logical(isFoldback)
    )]
    x
}

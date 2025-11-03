
# RE filtering_sites table columns
# this table has a header on disk
filteringSiteCols_ <- list( 
    chrom       = "integer", # i.e., rName, e.g., chr13
    sitePos1    = "integer",
    inSilico    = "integer", # 0/1 logical
    nObserved   = "integer"  # sitePos1 count in the training data set (not necessarily the current sample)
)
filteringSiteCols <- names(filteringSiteCols_)
filteringSiteColClasses <- unname(unlist(filteringSiteCols_))

# extract observed_endpoints columns
nodeCols_ <- list( 
    rNameIndex  = "integer", # e.g., 13
    rPos1       = "integer",
    strand      = "integer", # 0/1 for +/-, top/bottom
    nodeN       = "integer", # 1 or 2
    isJunction  = "integer", # 0/1 logical
    nObserved   = "integer"  # rPos1 count in the current sample
)
nodeCols <- names(nodeCols_)
nodeColClasses <- unname(unlist(nodeCols_))

# extract edges columns
edgeCols_ <- list( 
    #-------------  read level
    "qName"         = "character", # qName:readN, unsplit, as provided by the sequencing platform appended by readN for uniqueness of unmerged reads
    "channel"       = "integer",
    "trim5"         = "integer",
    "trim3"         = "integer",
    "mergeLevel"    = "integer",
    "nRead1"        = "integer",
    "nRead2"        = "integer",
    "isFoldback"    = "integer",
    #-------------  block level
    "blockN"        = "integer",
    #-------------  edge level
    "edgeN"         = "integer",
    "edgeType"      = "character",
    "nodePos1_1"    = "integer64",
    "qPos1_1"       = "integer",
    "nodePos1_2"    = "integer64",
    "qPos1_2"       = "integer",
    "eventSize"     = "integer",
    "mapQ"          = "integer",
    "alnQ"          = "double", # gapCompressedIdentity
    #------------- alignment level
    "cigar"         = "character", 
    #------------- junction level
    "flankLen"      = "integer",
    "insertSize"    = "integer",
    "jxnBases"      = "character"
)
edgeCols <- names(edgeCols_)
edgeColClasses <- unname(unlist(edgeCols_))

# edge types
edgeTypes <- list(
    ALIGNMENT     = "A", # the single type for a contiguous aligned segment
    TRANSLOCATION = "T", # edge/junction types (might be several per source molecule)
    INVERSION     = "V",
    DUPLICATION   = "U",
    DELETION      = "D",
    UNKNOWN       = "?",
    INSERTION     = "I", 
    PROPER        = "P"
)

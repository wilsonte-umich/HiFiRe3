# generic support for data formats and indexed data retrieval using tabix
# these functions expect that files utilized chromIndex1, not chrom

# types
af_strands <- list(
    top    = 0,
    bottom = 1
)
af_junctions <- list(
    indexToType = c(
        "deletion",
        "duplication",
        "inversion",
        "translocation",
        "mixed" # when a floating alignment has two junction types
    ),
    indexToTypeLabel = c(
        "del    ",
        "dup    ",
        "inv    ",
        "trans  ",
        "int-gen"
    ),
    typeToIndex = list(
        proper          = 0L,
        deletion        = 1L,
        duplication     = 2L,
        inversion       = 3L,
        translocation   = 4L,
        mixed           = 5L
    ),
    indexToRowOffset = c(
        deletion        = -1/3,
        duplication     = 0,
        inversion       = 1/3,
        translocation   = 0,
        mixed           = NA
    ),
    typeToColor = list(
        proper          = CONSTANTS$plotlyColors$black,
        deletion        = CONSTANTS$plotlyColors$blue,
        duplication     = CONSTANTS$plotlyColors$green,
        inversion       = CONSTANTS$plotlyColors$red,
        translocation   = CONSTANTS$plotlyColors$orange,
        mixed           = CONSTANTS$plotlyColors$purple
    ),
    barColors = function(){
        af_junctions$typeToColor[1:5 + 1]
    },
    getIndexColor = function(indices){
        sapply(indices, function(i){
            if(i == 0) CONSTANTS$plotlyColors$black
            else af_junctions$typeToColor[[af_junctions$indexToType[i]]]
        })
    },
    getJxnTypeIndices = function(jxnTypes){
        sapply(jxnTypes, function(jxnType){
            af_junctions$typeToIndex[[jxnType]]
        })
    },
    getJxnTypeNames = function(jxnIndices){
        sapply(jxnIndices, function(jxnIndex){
            af_junctions$indexToType[jxnIndex]
        })
    },
    getJxnTypeLabels = function(jxnIndices){
        sapply(jxnIndices, function(jxnIndex){
            af_junctions$indexToTypeLabel[jxnIndex]
        })
    }
)
af_junctionClasses <- list(
    indexToClass = c(
        "intergenome",
        "artifact",
        "validated",
        "expected",
        "unexpected2",
        "unexpected1"
    ),
    indexToClassLabel = c(
        "intergenome, all n",
        "artifact, n=1",
        "validated, n>=3",
        "expected, n<3",
        "unexpected, n=2",
        "unexpected, n=1"
    ),
    sizedLabels = c(
        "artifact, n=1",
        "validated, n>=3",
        "expected, n<3",
        "unexpected, n=2",
        "unexpected, n=1"
    ),
    classToIndex = list(
        intergenome = 1L,
        artifact    = 2L,
        validated   = 3L,
        expected    = 4L,
        unexpected2 = 5L,
        unexpected1 = 6L
    ),
    classToColor = list(
        intergenome = CONSTANTS$plotlyColors$purple,
        artifact    = CONSTANTS$plotlyColors$orange,
        validated   = CONSTANTS$plotlyColors$green,
        expected    = CONSTANTS$plotlyColors$blue,
        unexpected2 = CONSTANTS$plotlyColors$red,
        unexpected1 = CONSTANTS$plotlyColors$red
    ),
    getClassColor = function(indices){
        sapply(indices, function(i){
            if(i == 0) CONSTANTS$plotlyColors$black
            else af_junctionClasses$classToColor[[af_junctionClasses$indexToType[i]]]
        })
    },
    getJxnClassLabels = function(is){
        sapply(is, function(i){
            af_junctionClasses$indexToClassLabel[i]
        })
    }
)
af_junctionStrata <- list(
    indexToStratum = c(
        "unexpected1",
        "unexpected2",
        "expected",
        "validated"
    ),
    indexToStratumLabel = c(
        "unexpected, n=1",
        "unexpected, n=2",
        "expected, n<3",
        "validated, n>=3"
    ),
    stratumToIndex = list(
        unexpected1 = 1L,
        unexpected2 = 2L,
        expected    = 3L,
        validated   = 4L
    ),
    getJxnStratumLabels = function(is){
        sapply(is, function(i){
            af_junctionStrata$indexToStratumLabel[i]
        })
    }
)
af_alignments <- list(
    typeToColor = list(
        alignment       = CONSTANTS$plotlyColors$black,
        projection      = CONSTANTS$plotlyColors$grey
    ) 
)
af_haplotypes <- list(
    haplotypesToColor = c(
        CONSTANTS$plotlyColors$black,   # 1 = reference only
        CONSTANTS$plotlyColors$blue,    # 2 = haplotype 1
        CONSTANTS$plotlyColors$blue,    # 3
        CONSTANTS$plotlyColors$red,     # 4 = haplotype 2
        CONSTANTS$plotlyColors$red,     # 5
        CONSTANTS$plotlyColors$purple,  # 6 = haplotype 1+2 
        CONSTANTS$plotlyColors$purple   # 7
    ),
    getHaplotypeColor = function(haplotypes){
        af_haplotypes$haplotypesToColor[max(haplotypes, 1L)]
    }
)

# column definitions
af_bgzColumns <- list(
    filteringSitesBgz = c(
        chromIndex1     = "integer",
        sitePos1        = "integer",
        inSilico        = "integer",
        nObserved       = "integer"
    ),
    svAlignmentsBgz = c(
        chromIndex1     = "integer",
        refPos1_5       = "integer",
        refPos1_2       = "integer",
        refPos1_3       = "integer",
        strandIndex0_5  = "integer",
        jxnType         = "character", # middle floating alignments have two comma-delimited types
        pathN           = "integer",
        nJunctions      = "integer",
        alnN            = "integer",
        nObserved       = "integer"
    ),
    svJunctions1Bgz = c(
        chromIndex1_1   = "integer",
        refPos1_1       = "integer",
        strandIndex0_1  = "integer",
        chromIndex1_2   = "integer",
        refPos1_2       = "integer",
        strandIndex0_2  = "integer",
        jxnType         = "integer",
        nObserved       = "integer",
        alnOffset       = "character", # not integer since may be * (missing) when SEQ was dropped for known SVs
        jxnBases        = "character",
        paths           = "character",
        nPathJxns       = "integer",
        mapQ            = "integer",
        deTag           = "double",
        siteDist        = "integer",
        alnFailFlag     = "character",
        jxnFailFlag     = "character",
        qNames          = "character", # comma-delimited lists, one per observed read
        insertSizes     = "character", 
        isAllowedSizes  = "character", 
        stemLengths     = "character",
        passedStems     = "character",
        seqs            = "character",
        quals           = "character",
        cigars          = "character", # one per flanking alignment as CIGAR_1::CIGAR_2
        orientations    = "character",
        isExpected      = "integer",
        hasAltAlignment = "integer",
        svSize          = "integer",
        isIntergenome   = "integer",
        target1         = "character",
        targetDist1     = "integer",
        target2         = "character",
        targetDist2     = "integer",
        genes1          = "character",
        geneDist1       = "character", # comma-delimited list of distances
        genes2          = "character",
        geneDist2       = "character",
        isExcluded_1    = "integer",
        isExcluded_2    = "integer",
        bkptCoverage_1  = "integer",
        bkptCoverage_2  = "integer",
        nSamples        = "integer",
        samples         = "character"  # ,sample1,sample2,...
    ),
    svReadPaths = c(
        qName       = "character",
        qLen        = "integer",
        channel     = "integer",
        readHasJxn  = "integer",
        chroms      = "character",
        pos1s       = "character",
        strand0s    = "character",
        mapqs       = "character",
        deTags      = "character",
        blockNs     = "character",
        nRefBases   = "character",
        alnFailFlags= "character",
        jxnFailFlags= "character",
        qryStart0s  = "character",
        qryEnd1s    = "character",
        cigars      = "character",
        seq         = "character",
        qual        = "character",
        pseudoRef   = "character"
    )
)
af_bgzColumns$svJunctions2Bgz <- af_bgzColumns$svJunctions1Bgz
af_bgzColumns$svJunctions2Bgz_multi <- c(af_bgzColumns$svJunctions1Bgz, sample  = "integer")
af_bgzColumns$svJunctions2Bgz_multi <- af_bgzColumns$svJunctions1Bgz_multi # for multiple-sample rbind

# column display definitions
af_bgzColumns_display <- list(
    svJunctions1Bgz = c(
        sample          = "sample",
        chromIndex1_1   = "chrom1",
        refPos1_1       = "bkptPos1",
        strandIndex0_1  = "strand1",
        chromIndex1_2   = "chrom2",
        refPos1_2       = "bkptPos2",
        strandIndex0_2  = "strand2",
        svSize          = "svSize",
        jxnType         = "jxnType",
        isIntergenome   = "interGen",
        nObserved       = "nObs",
        isExpected      = "expected",
        mapQ            = "mapQ",
        deTag           = "deTag",
        siteDist        = "siteDist",
        alnOffset       = "alnOffset", # not integer since may be * (missing)
        # jxnBases        = "character",
        target1         = "target1",
        # targetDist1     = "integer",
        target2         = "target2",
        # isExcluded_1    = "excl1",
        # isExcluded_2    = "excl2",
        # hasAltAlignment = "hasAlt",
        # targetDist2     = "integer",
        # genes1          = "character",
        # geneDist1       = "character", # comma-delimited list of distances
        # genes2          = "character",
        # geneDist2       = "character",
        bkptCoverage_1  = "bkptCov1",
        bkptCoverage_2  = "bkptCov2",
        # orientations    = "character",
        insertSizes     = "insSizes", 
        stemLengths     = "stemLens",
        NULL
        # paths           = "character",
        # nPathJxns       = "integer",
        # qNames          = "character",
        # seqs            = "character",
        # quals           = "character",
        # cigars          = "character"
    )
)

# chrom to chromIndex conversion for bgz queries
# simple cache, sessionCache not used
af_chromIndex <- list()  
af_getChromIndex <- function(sourceId, chrom){
    if(is.null(af_chromIndex[[sourceId]])) {
        af_chromIndex[[sourceId]] <<- fread(getSourceFilePath(sourceId, "chromsFile"))
    }
    ci <- af_chromIndex[[sourceId]]
    ci[ci[[1]] == chrom | startsWith(ci[[1]], paste0(chrom, "_"))][[2]]
}
af_getChromNames <- Vectorize(function(sourceId, chromIndex){
    if(is.null(af_chromIndex[[sourceId]])) {
        af_chromIndex[[sourceId]] <<- fread(getSourceFilePath(sourceId, "chromsFile"))
    }
    af_chromIndex[[sourceId]][af_chromIndex[[sourceId]][[2]] == chromIndex][[1]]
})
af_getChromSize <- function(sourceId, chromIndex){
    if(is.null(af_chromIndex[[sourceId]])) {
        af_chromIndex[[sourceId]] <<- fread(getSourceFilePath(sourceId, "chromsFile"))
    }
    af_chromIndex[[sourceId]][af_chromIndex[[sourceId]][[2]] == chromIndex][[3]]
}

# genome excluded regions
af_excludedRegions <- list()
af_getExcludedRegions <- function(sourceId){
    if(is.null(af_excludedRegions[[sourceId]])) {
        af_excludedRegions[[sourceId]] <<- fread(getSourceFilePath(sourceId, "exclusionsBed"), fill = TRUE)[, 1:3]
        setnames(af_excludedRegions[[sourceId]], c("chrom", "start0", "end1"))
    }
    af_excludedRegions[[sourceId]]
}

# general track data retrieval, uses tabix random access not sessionCache
af_getTrackData_bgz <- function(sourceId, fileType, coord){
    coord$chromosome <- af_getChromIndex(sourceId, coord$chromosome)
    bgzFile <- getSourceFilePath(sourceId, fileType)
    if(!isTruthy(bgzFile) || !file.exists(bgzFile)) return(data.table()) # no data of this type
    getCachedTabix(bgzFile) %>% # create = TRUE, force = TRUE
    getTabixRangeData(
        coord, 
        col.names  =  names(af_bgzColumns[[fileType]]), 
        colClasses = unname(af_bgzColumns[[fileType]]), 
        skipChromCheck = TRUE
    ) 
    # coord$chromosome <- af_getChromIndex(sourceId, coord$chromosome)
    # getSourceFilePath(sourceId, fileType) %>%
    # getCachedTabix() %>%
    # getTabixRangeData(
    #     coord, 
    #     col.names  =  names(af_bgzColumns[[fileType]]), 
    #     colClasses = unname(af_bgzColumns[[fileType]]), 
    #     skipChromCheck = TRUE
    # )
}
# get and parse RE site, with padding for plotting
af_getSites_padded <- function(sourceId, coord){
    coord$start <- coord$start - coord$width
    coord$end   <- coord$end   + coord$width
    if(coord$start < 1) coord$start <- 1L
    af_getTrackData_bgz(sourceId, "filteringSitesBgz", coord)
}
# get and parse read alignments
af_getAlignments <- function(sourceId, coord){
    x <- af_getTrackData_bgz(sourceId, "svAlignmentsBgz", coord)
    if(nrow(x) == 0) return(x)
    x[
        refPos1_5 > 1 & 
        refPos1_2 > 0
    ]
    x[, jxnType := ifelse(
        grepl(",", jxnType),
        af_junctions$typeToIndex$mixed,
        as.integer(jxnType)
    )]
    x
}
# get and parse unique junction nodes in region ...
af_getJunctions <- function(sourceId, coord){
    # svJunctions1Bgz and svJunctions2Bgz are the same row data, just sorted differently
    # having both files allows for recover of all unique SVs with either junction within coord
    # importantly, svJunctions1Bgz and svJunctions2Bgz do NOT! have reversed nodes!
    # svJunctions1Bgz node1 == svJunctions2Bgz node1 NOT(svJunctions2Bgz node2)
    rbind(
        af_getTrackData_bgz(sourceId, "svJunctions1Bgz", coord),
        af_getTrackData_bgz(sourceId, "svJunctions2Bgz", coord) # same file as svJunctions1Bgz, just sorted differently
    ) %>% unique()
}
# ... additionally filtered the same way as is active the in exploreJunctions step
af_getFilteredJunctions <- function(sourceId, coord){
    af_applyJunctionFilters(
        af_getJunctions(sourceId, coord), 
        app$exploreJunctions$settingsObject
    )
}

# loading entire junction and tally tables, cached in sessionCache
af_jxnClassNames <- c("all_jxns", "validated", "expected", "artifact", "intergenome", "candidate")
af_jxnClasses <- 1:length(af_jxnClassNames)
names(af_jxnClasses) <- af_jxnClassNames
af_jxnCreate <- "asNeeded"
af_getJunctions_all_source <- function(sourceId){
    fileType <- "svJunctions1Bgz"
    filePath <- loadPersistentFile(
        sourceId = sourceId, # ... or source
        contentFileType = fileType,
        header = FALSE,
        force = FALSE,
        colClasses = unname(af_bgzColumns[[fileType]]), 
        col.names = names(af_bgzColumns[[fileType]]), # additional argument passed to fread
        quote = "", # additional argument passed to fread
        postProcess = function(jxns){
            startSpinner(session, message = "post-process junctions")
            jxns[, ":="(
                jxnI = 1:.N,
                sourceId = sourceId,
                alnOffset = sapply(alnOffset, function(x) if(x == "*") NA_integer_ else as.integer(x)),
                sample = sub(".analyze.SVs", "", getSourceFilePackageName(sourceId)),
                isValidated = nObserved >= 3,
                isExpected  = isExpected == 1
            )][, ":="(
                isArtifact = nObserved == 1 & !isExpected & (
                    jxnType == af_junctions$typeToIndex$translocation |
                    svSize > 1e6
                )
            )][, ":="(
                jxnClass = ifelse(
                    isIntergenome, af_junctionClasses$classToIndex$intergenome,
                    ifelse(
                        isArtifact, af_junctionClasses$classToIndex$artifact,
                        ifelse(
                            isValidated, af_junctionClasses$classToIndex$validated,
                            ifelse(
                                isExpected, af_junctionClasses$classToIndex$expected,
                                ifelse(
                                    nObserved == 2, af_junctionClasses$classToIndex$unexpected2,
                                    af_junctionClasses$classToIndex$unexpected1
                                )
                            )
                        )
                    )
                ),
                jxnStratum = ifelse(
                    isValidated, af_junctionStrata$stratumToIndex$validated,
                    ifelse(
                        isExpected, af_junctionStrata$stratumToIndex$expected,
                        ifelse(
                            nObserved == 2, af_junctionStrata$stratumToIndex$unexpected2,
                            af_junctionStrata$stratumToIndex$unexpected1
                        )
                    )
                )
            )]
            jxns
        }
    )
    persistentCache[[filePath]]$data
}
af_getJunctions_all_sources <- function(sourceIds){
    if(length(sourceIds) == 1){
        af_getJunctions_all_source(sourceIds)
    } else {
        fileType <- "svJunctions1Bgz"
        sessionCache$get(
            fileType, 
            keyObject = list(sourceIds, "af_getJunctions_all_sources"), 
            permanent = TRUE,
            from = "ram",
            create = af_jxnCreate,
            createFn = function(...) {
                startSpinner(session, message = "parsing junctions")
                x <- do.call(rbind, lapply(sourceIds, af_getJunctions_all_source))
                x[, jxnI := 1:.N]
                x
            }  
        )$value
    }
}

# # SV read paths
# af_rpCreate <- "asNeeded"
# af_getReadPaths_source <- function(sourceId){
#     fileType <- "svReadPaths"
#     filePath <- loadPersistentFile(
#         sourceId = sourceId, # ... or source
#         contentFileType = fileType,
#         header = FALSE,
#         force = FALSE,
#         colClasses = unname(af_bgzColumns[[fileType]]), 
#         col.names = names(af_bgzColumns[[fileType]]), # additional argument passed to fread
#         quote = "", # additional argument passed to fread
#         postProcess = function(rps){
#             startSpinner(session, message = "post-process read paths")
#             rps[, ":="(
#                 sourceId = sourceId
#             )]
#         }
#     )
#     persistentCache[[filePath]]$data
# }
# af_getReadPaths_sources <- function(sourceIds){
#     if(length(sourceIds) == 1){
#         af_getReadPaths_source(sourceIds)
#     } else {
#         fileType <- "svReadPaths"
#         sessionCache$get(
#             fileType, 
#             keyObject = list(sourceIds, "af_getReadPaths_sources"), 
#             permanent = TRUE,
#             from = "ram",
#             create = af_rpCreate,
#             createFn = function(...) {
#                 startSpinner(session, message = "parsing read paths")
#                 x <- do.call(rbind, lapply(sourceIds, af_getReadPaths_source))
#                 # x[, rpI := 1:.N]
#                 x
#             }  
#         )$value
#     }
# }

# get a single SV-containing read path by qName
getTaskFile <- function(sourceId, suffix){
    outputDir <- getSourcePackageOption(sourceId, "output", "output-dir")
    dataName  <- getSourcePackageOption(sourceId, "output", "data-name")
    taskDir   <- file.path(outputDir, dataName)
    fileName  <- paste0(dataName, suffix)
    file.path(taskDir, fileName) %>% Sys.glob # wildcard expansion
}
af_getReadPath <- function(sourceId, qName){
    req(sourceId, qName)
    bgzFile <- getTaskFile(sourceId, '.*.read_paths.txt.bgz')
    tabix   <- getCachedTabix(bgzFile)
    gRange <- GenomicRanges::GRanges(
        seqnames = qName, # qName masquerades as rName for this purpose
        ranges = IRanges::IRanges(start = 0, end = 1e6) # coordinates are irrelevant
    )
    lines <- Rsamtools::scanTabix(tabix, param = gRange)[[1]]
    if(length(lines) == 1) lines <- paste0(lines, "\n")
    fread(
        text = lines,
        col.names  =  names(af_bgzColumns$svReadPaths), 
        colClasses = unname(af_bgzColumns$svReadPaths), 
        header = FALSE,
        sep = "\t"
    )
}

# get a single junction with SEQ and QUAL
af_getFullJunction <- function(jxn){
    req(jxn)
    bgzFile <- getTaskFile(jxn$sourceId, '.*.final_junctions_1.txt.bgz')
    tabix <- getCachedTabix(bgzFile)
    gRange <- GenomicRanges::GRanges(
        seqnames = jxn$chromIndex1_1,
        ranges = IRanges::IRanges(start = jxn$refPos1_1 - 1, end = jxn$refPos1_1 + 1)
    )
    lines <- Rsamtools::scanTabix(tabix, param = gRange)[[1]]
    if(length(lines) == 1) lines <- paste0(lines, "\n")
    fjxns <- fread(
        text = lines,
        col.names  =  names(af_bgzColumns$svJunctions1Bgz), 
        colClasses = unname(af_bgzColumns$svJunctions1Bgz), 
        header = FALSE,
        sep = "\t"
    )
    fjxns[
        chromIndex1_1 == jxn$chromIndex1_1 & 
        refPos1_1     == jxn$refPos1_1 & 
        chromIndex1_2 == jxn$chromIndex1_2 & 
        refPos1_2     == jxn$refPos1_2
    ]
}


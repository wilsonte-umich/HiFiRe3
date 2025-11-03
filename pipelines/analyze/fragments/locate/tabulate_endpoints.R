# action:
#     compare observed sequence endpoints to each other and to the 
#       in silico digested genome to identify presumptive RE sites
#       for later use during error correction
# steps:
#     compare observed and in silico sitePos1 values to each other in parallel by chromosome
#       use ACCEPT_ENDPOINT_DISTANCE to allow sitePos1 values to match if differing by a few bp
#       accounts for inaccuracies in observed sitePos1 resulting from base errors, indels, etc.
#     declare whether observed sequence endpoints define a cluster and thus presumptive RE sites
#     identify whether presumptive RE sites are known in silico sites or apparent RFLPs
#       require --min-rflp-evidence matching endpoints to accept an unexpected position cluster as an RFLP
# output: 
#     headered "filteringSites" table (in RDS) and file (txt.gz) of the unique sitePos1 where isFilteringSite == TRUE
#           chrom,sitePos1,inSilico,nObserved
#       where:
#           chrom is the named chrom, e.g., chr1
#           sitePos1 is the assigned consensus RE site position of the sequence endpoint position cluster
#           inSilico is a logical flag whether sitePos1 matches an RE site in the reference genome
#           nObserved is the total number of informative sequence endpoints assigned to sitePos1
#     headered "sampleSites" table (in RDS) of site positions observed in the sample data (i.e., nObserved > 0)
#           chrom,sitePos1,nObserved,inSilico,isFilteringSite
#       where, in addition to definitions above:
#           isFilteringSite is a logical flag whether sitePos1 will be used for endpoint-to-site matching in pipeline
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
        'GENOME_SITES_GZ',
        'GENOME',
        'BLUNT_RE_TABLE',
        'ENZYME_NAME',
        'OBSERVED_ENDPOINTS_FILE',
        'PLOT_PREFIX',
        'APP_ENDPOINTS_FILE',
        'FILTERING_SITES_FILE'
    ),
    integer = c(
        'ACCEPT_ENDPOINT_DISTANCE',
        'MIN_RFLP_EVIDENCE',
        'N_CPU'
    )
))
#-------------------------------------------------------------------------------------
# source R scripts
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
sourceScripts(file.path(rUtilDir, 'genome'), c('chroms'))
setCanonicalChroms()
#-------------------------------------------------------------------------------------
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#=====================================================================================

#=====================================================================================
# feedback functions
#-------------------------------------------------------------------------------------
plotIt <- function(plotName, plotFn, ...){
    png(filename = paste(env$PLOT_PREFIX, plotName, "png", sep = "."),
        width = 3, height = 3, units = "in", pointsize = 8,
        bg = "white",  res = 600, type = "cairo")
    plotFn(plotName, ...)
    dev.off()
}
reportIt <- function(label, value) message(paste("", value, label, sep = "\t"))
percent <- function(num, denom, digits = 2) round(num / denom * 100, digits)
#-------------------------------------------------------------------------------------
plotSitePosEvidence <- function(plotName, sitePosCounts, log){
    peakCount <- sitePosCounts[!is.na(reference) & reference > 0, weighted.mean(nObserved, reference)]
    maxX <- peakCount * 3
    maxY <- max(sitePosCounts$reference)
    getX <- function(x, y) if(log)       x[y > 0]  else x
    getY <- function(y)    if(log) log10(y[y > 0]) else y
    plot(
        x = NA,
        y = NA,
        xlab = "Site Endpoint Count",
        ylab = paste(if(log) "log10 " else "", "# of Observed Sites"),
        xlim = c(0, maxX),
        ylim = range(getY(0:maxY)) * 1.05,
        main = paste(sitePosCounts[nObserved == 1, rflp], "non-reference singletons")
    )
    abline(v = env$MIN_RFLP_EVIDENCE)
    y <- sitePosCounts$reference
    points(
        x = getX(sitePosCounts$nObserved, y),
        y = getY(y),
        col = "blue",
        pch = 19,
        cex = 0.5
    )
    y <- sitePosCounts$rflp
    points(
        x = getX(sitePosCounts$nObserved, y),
        y = getY(y),
        col = "red3",
        pch = 19,
        cex = 0.5
    )
}
#=====================================================================================

#=====================================================================================
# load the RE sites expected from in silico digestion of the genome
#-------------------------------------------------------------------------------------
sitePosCols  <- c("chrom","sitePos1")

message("gathering cut site positions from in silico digestion")
# enzyme,strand,cut_site,regex,offset,CpG_priority,high_fidelity,site_length,degenerate,effective_length,CpG_sensitive,CpG_level,units_ul,units,cost,cents_unit,temperature,buffer,heat_inactivation,star_activity,one_hour_limit,HiFiRe3_compatible,initial_priority,comments
res <- fread(env$BLUNT_RE_TABLE)
res[, reKey := paste(gsub("\\s", "", enzyme), toupper(cut_site), sep = "_")]
reKey <- res[enzyme == env$ENZYME_NAME, reKey]
sitesFile <- env$GENOME_SITES_GZ
inSilicoPositions <- fread(sitesFile)
names(inSilicoPositions) <- sitePosCols
inSilicoPositions[, inSilico := TRUE]
#=====================================================================================

#=====================================================================================
# load the sequence endpoints observed in actual alignment data
#-------------------------------------------------------------------------------------
message("gathering sequence endpoints")
observedPositions <- fread(env$OBSERVED_ENDPOINTS_FILE)
names(observedPositions) <- c("chrom", "sitePos1", "nObserved")
#=====================================================================================

#=====================================================================================
# report summary statistics to log
#-------------------------------------------------------------------------------------
message("in silico summary statistics")
nInSilicoPostions <- inSilicoPositions[, .N]
reportIt("genome", env$GENOME)
reportIt("restriction enzyme", env$ENZYME_NAME)
reportIt("# in silico positions", nInSilicoPostions)

message("observed endpoint summary statistics")
nObservedPositions <- observedPositions[, .N]
nEndpoints   <- observedPositions[, sum(nObserved)]
reportIt("# input positions", nObservedPositions)
reportIt("# input endpoints", nEndpoints)
#=====================================================================================

#=====================================================================================
# attempt to fuzzy-match observed to in silico position within ACCEPT_ENDPOINT_DISTANCE
# or, if not possible, commit adjusted positions as candidate RE sites
#-------------------------------------------------------------------------------------
message("fuzzy matching observed endpoint and in-silico RE site positions")
OBS_TYPES <- list(inSilico = 0, observed = 1)

# assign closely spaced observed positions to the most likely candidate RE position
commitNonBreakingGroup <- function(ordPos){
    nInSilico <- ordPos[, sum(obsType == OBS_TYPES$inSilico)]

    # only unmatched in silico sites remain, ignore them (will add these back later)
    if(nInSilico == ordPos[, .N]) { 
        NULL

    # no matching in silico positions to one or more observed positions
    # commit as the most frequent observed position
    } else if (nInSilico == 0){
        ordPos[, .(
            sitePos1 = ordPos[, sitePos1[which.max(nObserved)]],
            inSilico = FALSE,
            nObserved = nObserved
        )]

    # exactly one matching in silico position to one or more observed positions
    # commit using the in silico position as the candidate RE position
    } else if (nInSilico == 1){
        ordPos[obsType == OBS_TYPES$observed, .(
            sitePos1 = ordPos[obsType == OBS_TYPES$inSilico, sitePos1],
            inSilico = TRUE,
            nObserved = nObserved
        )]

    # exactly one observed position with multiple matching in silico positions
    # commit using the best matching in silico position
    } else if (ordPos[, sum(obsType == OBS_TYPES$observed) == 1]){
        inSilPos <- ordPos[obsType == OBS_TYPES$inSilico]
        obsPos   <- ordPos[obsType == OBS_TYPES$observed]
        dist <- inSilPos[, abs(obsPos$sitePos1 - sitePos1)]
        obsPos[, .(
            sitePos1 = inSilPos[, sitePos1[which.min(dist)]],
            inSilico = TRUE,
            nObserved = nObserved
        )]

    # multiple matching in silico positions to multiple observed positions
    # hopefully uncommon but can occur, especially with closely spaced 4-base cutters
    # match each observed position to its best matching in silico position
    } else {
        inSilPos <- ordPos[obsType == OBS_TYPES$inSilico]
        obsPos   <- ordPos[obsType == OBS_TYPES$observed]
        distanceMatrix <- sapply(inSilPos$sitePos1, function(sitePos1) abs(obsPos$sitePos1 - sitePos1))
        do.call(rbind, lapply(1:nrow(obsPos), function(i){
            j <- which.min(distanceMatrix[i,])
            obsPos[i, .(
                sitePos1 = inSilPos[j, sitePos1],
                inSilico = TRUE,
                nObserved = nObserved
            )]
        }))
    }
}

# split apart groups that didn't break break sequentially but still need to break
# use kmeans to iteratively break each overly large group into two smaller groups until all subgroups resolve
splitPositionGroup <- function(ordPos){

    # all remaining positions are a coherent non-breaking group
    if(ordPos[, diff(range(sitePos1)) <= env$ACCEPT_ENDPOINT_DISTANCE]){
        commitNonBreakingGroup(ordPos)

    # group still needs breaking; do it and try each part again
    } else {
        cluster <- kmeans(ordPos$sitePos1, 2)$cluster
        rbind(
            splitPositionGroup(ordPos[cluster == 1]),
            splitPositionGroup(ordPos[cluster == 2])
        )
    }
}

# work in parellel by chromosome
matchResults <- do.call(rbind, mclapply(unique(observedPositions$chrom), function(chrom_){
# matchResults <- do.call(rbind, lapply(unique(observedPositions$chrom), function(chrom_){
# matchResults <- do.call(rbind, mclapply(c("chr21","chr22"), function(chrom_){

    message(paste(" ", chrom_))
    orderedPos <- rbind( # generate an ordered list of all candidate RE site positions on the chromosome ...
        inSilicoPositions[chrom == chrom_, .(  # ... from in silico digest
            sitePos1,
            obsType = OBS_TYPES$inSilico,
            nObserved = NA_integer_
        )],
        observedPositions[chrom == chrom_, .( # ... and from observed alignment ends
            sitePos1,
            obsType = OBS_TYPES$observed,
            nObserved = nObserved
        )]
    )[order(sitePos1)]
    interPosDistances <- orderedPos[, .( # calculate the coordinate distances between all positions
        leftPosI  = 1:(.N - 1),
        rightPosI = 2:.N,
        distance = diff(sitePos1))
    ]
    interPosDistances[, breaking := distance > env$ACCEPT_ENDPOINT_DISTANCE] # determine if position pairs are within the site tolerance 
    breakingPosGroups <- interPosDistances[, .( # declare runs of non-breaking and breaking positions
        leftPosI  = min(leftPosI),
        rightPosI = max(rightPosI),
        nPos = .N + 1 # i.e., one more than the number of distances
    ), by = .(breaking, rleid(breaking))]
    breakingPosGroups <- breakingPosGroups[breaking == FALSE | nPos > 2] # adjust groups to find positions that broke on both sides
    breakingPosGroups[breaking == TRUE, ":="(
        leftPosI  = leftPosI  + 1,
        rightPosI = rightPosI - 1
    )]
    breakingPosGroups[, {
        ordPos <- orderedPos[leftPosI:rightPosI]

        # one or more ordered positions that all broke relative to each other
        if(breaking){ 
            ordPos[obsType == OBS_TYPES$observed, .( # ignore in silico sites, they had no matching observed
                inSilico = FALSE,
                nObserved = nObserved
            ), by = .(sitePos1)] # commit observed positions individually

        # one or more ordered positions that did not break sequentially
        # but that may still break between first and last positions
        } else {
            splitPositionGroup(ordPos)
        }
    }, by = .(rleid)][, .(
        chrom = chrom_,
        sitePos1 = sitePos1,
        inSilico = inSilico,
        nObserved = nObserved
    )]
}, mc.cores = env$N_CPU))
# }))
#=====================================================================================

#=====================================================================================
# more reporting
#-------------------------------------------------------------------------------------
message("plotting RE site position evidence counts by in silico site match")
sitePosCounts <- dcast(
    matchResults[
        !(chrom %in% c("chrY","chrM","chrEBV")), 
        .(nObserved = sum(nObserved)), 
        by = .(chrom, sitePos1, inSilico)
    ], 
    nObserved ~ inSilico, 
    fun.aggregate = length, 
    value.var = "sitePos1"
)
names(sitePosCounts) <- c("nObserved","rflp","reference")
invisible(plotIt(
    "siteEvidenceCounts_linear", 
    plotSitePosEvidence, 
    sitePosCounts,
    FALSE
))
invisible(plotIt(
    "siteEvidenceCounts_log", 
    plotSitePosEvidence, 
    sitePosCounts,
    TRUE
))
#=====================================================================================

#=====================================================================================
# construct the output tables
#-------------------------------------------------------------------------------------
message("assembling table of RE site positions seen as sequence endpoint clusters")
sampleSites <- matchResults[, # output has one row per assigned endpoint cluster position
    .(                        # i.e., fuzzy matches are now collapsed to their consensus position
        nObserved = sum(nObserved),
        inSilico = any(inSilico)
    ), 
    keyby = sitePosCols
]
sampleSites[,
    isFilteringSite := inSilico == TRUE |  # keep any cluster that matched a reference RE site
        nObserved >= env$MIN_RFLP_EVIDENCE # or any non-reference position with a significant endpoint cluster (an RFLP)
]                                          # drop low-coverage positions not at expected sites, nearly always random single fragments or low quality bases

message("assembling table of RE site positions for use in SV filtering (all in silico sites + observed RFLPs)")
filteringSites <- merge( # output has one row per assigned endpoint cluster position plus unmatched in silico sites
    inSilicoPositions,
    sampleSites[isFilteringSite == TRUE, .(chrom, sitePos1, nObserved)],
    all = TRUE,
    by = sitePosCols
)[order(chrom, sitePos1)]
filteringSites[is.na(inSilico),  ":="(inSilico = FALSE)]
filteringSites[is.na(nObserved), nObserved := 0L]

message("final site tallies (chr1 to chrX only)")
allSites <- merge( # output has one row per genome position of interest
    inSilicoPositions,
    sampleSites[, .(chrom, sitePos1, nObserved, isFilteringSite)],
    all = TRUE,
    by = sitePosCols
)[order(chrom, sitePos1)]
allSites[is.na(inSilico),  ":="(inSilico = FALSE)]
allSites[is.na(nObserved), nObserved := 0L]
allSites[is.na(isFilteringSite), isFilteringSite := inSilico]
allSites[, observed  := nObserved > 0]
allSites[, passedCov := nObserved >= env$MIN_RFLP_EVIDENCE]
print(allSites[
    !(chrom %in% c("chrY","chrM","chrEBV")), 
    .(
        nSites = .N,
        medianNObs = median(nObserved)
    ), 
    keyby = .(inSilico, observed, passedCov, isFilteringSite)
])
#=====================================================================================

#=====================================================================================
# write output files
#-------------------------------------------------------------------------------------

# endpoints as RDS for plotting in app
saveRDS(
    list(
        env = env,
        sampleSites = sampleSites[nObserved > 0],
        filteringSites = filteringSites
    ), 
    file = env$APP_ENDPOINTS_FILE
)

# filtering sites as txt.gz for use in various types of downstream code
fwrite(
    filteringSites, 
    file = env$FILTERING_SITES_FILE, 
    quote = FALSE,
    sep = "\t",
    logical01 = TRUE,
    scipen = 999L,
    nThread = env$N_CPU,
    compress = "gzip"
)
#=====================================================================================

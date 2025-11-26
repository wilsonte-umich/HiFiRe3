#=====================================================================================
# script initialization
#-------------------------------------------------------------------------------------
message("initializing")
suppressPackageStartupMessages(suppressWarnings({
    library(data.table)
}))
#-------------------------------------------------------------------------------------
# load, parse and save environment variables
env <- as.list(Sys.getenv())
rUtilDir <- file.path(env$MODULES_DIR, 'utilities', 'R')
source(file.path(rUtilDir, 'workflow.R'))
checkEnvVars(list(
    string = c(
        'GENOME_REMAPS_DIR',
        'GENOME',
        'BLUNT_RE_TABLE',
        'RE_SUMMARIES_RDS',
        'RE_DISTRIBUTIONS_RDS',
        'RE_FRAGMENTS_RDS'
    ),
    integer = c(
        'MIN_PRIORITY_LEVEL'
    )
))
#-------------------------------------------------------------------------------------
# source R scripts
sourceScripts(file.path(rUtilDir, 'genome'), c('stats'))
#-------------------------------------------------------------------------------------
# set some options
options(scipen = 999) # prevent 1+e6 in printed, which leads to read.table error when integer expected
options(warn = 2) 
#=====================================================================================

#=====================================================================================
# load the REs
#-------------------------------------------------------------------------------------
message("gathering restriction enzymes")
# enzyme,strand,cut_site,regex,offset,CpG_priority,high_fidelity,site_length,degenerate,effective_length,CpG_sensitive,CpG_level,units_ul,units,cost,cents_unit,temperature,buffer,heat_inactivation,star_activity,one_hour_limit,HiFiRe3_compatible,initial_priority,comments
res <- fread(env$BLUNT_RE_TABLE)[CpG_priority >= env$MIN_PRIORITY_LEVEL]
res[, reKey := paste(gsub("\\s", "", enzyme), toupper(cut_site), sep = "_")]
#=====================================================================================

#=====================================================================================
# turn site positions into fragment spans
#-------------------------------------------------------------------------------------
message("parsing restriction fragments")
frags <- sapply(res$reKey, function(reKey_){
    message(paste("", reKey_, sep = "\t"))
    filePrefix <- paste(env$GENOME, "digest", reKey_, sep = ".")
    sitesFile <- paste(filePrefix, "txt.gz", sep = ".")
    sitesFile <- file.path(env$GENOME_REMAPS_DIR, sitesFile)
    rePositions <- fread(sitesFile)
    names(rePositions) <- c("chrom","pos1")
    rePositions[
        !grepl("_.+_", chrom) # for composite genomes, only consider the primary chromosomes (o/w overruns memory)
    ][,      
        if(.N > 1) .(  
            start0 = pos1[1:(.N - 1)] - 1L, # known limitation, does not score fragments from cut site to chromosome ends 
            end1   = pos1[2:.N]             # only significant for small genomes with rare cut sites
        ) else .(
            start0 = integer(), 
            end1   = integer()
        ), 
        by = .(chrom)
    ]
}, simplify = FALSE, USE.NAMES = TRUE)
#=====================================================================================

#=====================================================================================
message("summarizing fragments")
distributions <- list()
lengths_ <- sapply(res$reKey, function(reKey_){
    frags[[reKey_]][, end1 - start0]
}, simplify = FALSE, USE.NAMES = TRUE)
avg_MW_g_mol_bp <- 615.96
one_ug_in_g <- 1e-6
distributions <- sapply(res$reKey, function(reKey_){
    x <- data.table(
        length = lengths_[[reKey_]]
    )[, .(
            N = .N,
            bp = length * .N,
            avg_MW_g_mol = length * avg_MW_g_mol_bp
        ), 
        keyby = .(length)
    ][, ":="(
        molar_fraction  = N  / sum(N),
        mass_fraction   = bp / sum(bp),
        avg_fmol_per_ug = one_ug_in_g / avg_MW_g_mol * 1e15
    )][, ":="(
        avg_fmol_per_ug = avg_fmol_per_ug * mass_fraction
    )][, ":="(
        bp = NULL,
        avg_MW_g_mol = NULL
    )]
}, simplify = FALSE, USE.NAMES = TRUE)
siteSummary <- data.table(
    reKey   = res$reKey,
    nSites  = sapply(res$reKey, function(k) length(lengths_[[k]])),
    totalGb = round(sapply(res$reKey, function(k) sum(lengths_[[k]])) / 1e9, 2),
    minLen  = sapply(res$reKey, function(k) min(lengths_[[k]])),
    maxLen  = sapply(res$reKey, function(k) max(lengths_[[k]]))
)
for(percent in c(1,2.5,5,10,50,90,95,97.5,99)){
    fraction <- percent / 100
    suffix <- as.character(percent)
    siteSummary[[paste0("q", suffix)]] <- sapply(res$reKey, function(reKey_){
        quantile(lengths_[[reKey_]], fraction)
    })
    siteSummary[[paste0("N", suffix)]] <- sapply(res$reKey, function(reKey_){
        N50(lengths_[[reKey_]], fraction)
    })
}
for(minFragSize_kb in c(0.145, 1:15)){
    minFragSize <- minFragSize_kb * 1000
    maxFragSize <- minFragSize * 2
    siteSummary[[paste("molarFrac", minFragSize, maxFragSize, sep = "_")]] <- sapply(res$reKey, function(reKey_){
        distributions[[reKey_]][length >= minFragSize & length < maxFragSize, sum(molar_fraction)]
    })
    siteSummary[[paste("massFrac", minFragSize, maxFragSize, sep = "_")]] <- sapply(res$reKey, function(reKey_){
        distributions[[reKey_]][length >= minFragSize & length < maxFragSize, sum(mass_fraction)]
    })
    siteSummary[[paste("N", minFragSize, maxFragSize, sep = "_")]] <- sapply(res$reKey, function(reKey_){
        distributions[[reKey_]][length >= minFragSize & length < maxFragSize, sum(N)]
    })
}
siteSummaries <- merge(
    res,
    siteSummary,
    by = "reKey",
    all = TRUE
)
#=====================================================================================

#=====================================================================================
# save results for app
message("saving results for app")
saveRDS(
    list(
        enzymes = res,
        siteSummaries = siteSummaries,
        env = env
    ),
    file = env$RE_SUMMARIES_RDS
)
saveRDS(
    distributions,
    file = env$RE_DISTRIBUTIONS_RDS
)
saveRDS(
    frags,
    file = env$RE_FRAGMENTS_RDS
)
#=====================================================================================

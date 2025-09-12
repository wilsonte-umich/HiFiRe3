# utility functions for loading RE site data and related
# caller must call stopSpinner()
snvData_create <- "asNeeded"

snvSummaryColumns <- c(
    SNV_CHROM_INDEX1 = "integer",
    SNV_START0 = "integer",
    SNV_END1 = "integer",
    SNV_REF_BASES = "character",
    SNV_ALT_BASES = "character",
    SNV_N_TOP_STRAND = "integer",
    SNV_N_BOTTOM_STRAND = "integer",
    SNV_N_HIGH_QUAL = "integer",
    SNV_N_OBSERVED = "integer",
    SNV_COVERAGE_TOP_STRAND = "integer",
    SNV_COVERAGE_BOTTOM_STRAND = "integer",
    SNV_COVERAGE_UNSTRANDED = "integer",
    SNV_ZYGOSITY = "integer",
    SNV_GENOTYPE_COVERAGE_UNSTRANDED = "integer",
    SNV_GENOTYPE_ZYGOSITY = "integer"
)
zygosityTypeIs <- list(
    na = 0,
    singleton = 1,
    subclonal = 2,
    heterozygous = 3,
    homozygous = 4,
    not_in_genotype = 5
)
zygosityTypes <- names(zygosityTypeIs)
zygosityTypeLabels <- c(
    "Single",
    "Subclonal",
    "Heterozyg",
    "Homozyg",
    "Missing"
)
varTypeIs <- list(
    snp     = 1, # single-nucleotide polymorphism, e.g. A>G
    del1    = 2, # deletion of 1 base, e.g. A>-
    ins1    = 3, # insertion of 1 base, e.g. ->A
    indel0  = 4, # equal insertion + deletion of >1 base; a multi-nucleotide polymorphism, e.g. AT>CG
    delN    = 5, # deletion of >1 base, e.g. AT>-
    insN    = 6, # insertion of >1 base, e.g. ->AT
    indelX  = 7, # unequal insertion + deletion of >1 base, e.g. AT>CGT, catch-all for complex events
    # -----------
    match   = 8, # no variant, e.g. A>A
    clipped = 9, # 5'-clipped base
    lowQual = 10 # masked low quality base
)
varTypeColors <- c(
    snp = "blue",
    del1 = "red",
    ins1 = "green",
    indel0 = "purple",
    delN = "orange",
    insN = "cyan",
    indelX = "brown",
    match = "gray",
    clipped = "white",
    lowQual = "white"
)
pileupCodes <- list(
    CS_MATCH        = ":",
    CLIPPED         = "!",
    MASKED_LOW_QUAL = "q"
)
pileupCodeVarTypeIs <- c(
    ":" = varTypeIs$match,
    "!" = varTypeIs$clipped,
    "q" = varTypeIs$lowQual
)
pileupCodeVals <- names(pileupCodeVarTypeIs)
varTypes <- names(varTypeIs)
parseVarType <- function(refBases, altBases){
    nRefBases <- nchar(refBases)
    nAltBases <- nchar(altBases)
    delta <- nAltBases - nRefBases
    ifelse(
        refBases %in% pileupCodeVals,
        pileupCodeVarTypeIs[refBases],
        ifelse(
            delta == 0,
            ifelse(nRefBases == 1, varTypeIs$snp, varTypeIs$indel0),
            ifelse(
                delta == 1,
                ifelse(nRefBases == 0, varTypeIs$ins1, varTypeIs$indelX),
                ifelse(
                    delta == -1,
                    ifelse(nRefBases == 1, varTypeIs$del1, varTypeIs$indelX),
                    ifelse(
                        delta > 1,
                        ifelse(nRefBases == 0, varTypeIs$insN, varTypeIs$indelX),
                        ifelse(nAltBases == 0, varTypeIs$delN, varTypeIs$indelX)
                    )
                )
            )
        )
    )
}
parseVarType_long <- function(varType, refBases, altBases){
    nRefBases <- nchar(refBases)
    nAltBases <- nchar(altBases)
    delta <- nAltBases - nRefBases
    ifelse(
        varType %in% c(varTypeIs$snp, varTypeIs$del1, varTypeIs$ins1),
        paste(refBases, altBases, sep = ">"),
        ifelse(
            varType %in% c(varTypeIs$indel0, varTypeIs$delN),
            nRefBases,
            ifelse(
                varType == varTypeIs$insN,
                nAltBases,
                pmax(nRefBases, nAltBases)
            )
        )
    )
}

fillSnvColumns <- function(vars){
    vars[, maskedZygosity := ifelse(
        SNV_GENOTYPE_ZYGOSITY == zygosityTypeIs$not_in_genotype, 
        SNV_ZYGOSITY, 
        SNV_GENOTYPE_ZYGOSITY # known SNV/indels gety promoted to the expected genotype zygosity
    )]
    vars[, ":="(
        isValidated = maskedZygosity >= zygosityTypeIs$heterozygous,
        vaf = SNV_N_OBSERVED / SNV_COVERAGE_UNSTRANDED,
        varType = parseVarType(SNV_REF_BASES, SNV_ALT_BASES)
    )]
    vars[, ":="(
        varType_long = parseVarType_long(varType, SNV_REF_BASES, SNV_ALT_BASES)
    )]    
}
snvs_loadSnvSummary <- function(sourceId){
    sessionCache$get(
        'snvs_loadSnvSummary', 
        key = sourceId, 
        permanent = TRUE,
        from = "ram",
        create = snvData_create,
        createFn = function(...) {
            startSpinner(session, message = "loading SNV summary")
            bgzFile <- getSourceFilePath(sourceId, "snvSummaryBgz")
            gzFile <- paste(bgzFile, "gz", sep = ".")
            file.copy(bgzFile, gzFile, overwrite = FALSE)
            vars <- fread(
                gzFile, header = FALSE, sep = "\t",
                col.names = names(snvSummaryColumns),
                colClasses = unname(snvSummaryColumns)
            ) %>% fillSnvColumns()
        }
    )$value
}

snvs_loadReadPileup <- function(sourceId){
    sessionCache$get(
        'snvs_loadReadPileup', 
        key = sourceId, 
        permanent = TRUE,
        from = "ram",
        create = snvData_create,
        createFn = function(...) {
            startSpinner(session, message = "loading read pileup")
            gzFile <- getSourceFilePath(sourceId, "readPileupGz")
            rp <- fread(gzFile, header = TRUE, sep = "\t")
            rp <- rp[, 
                .(nObserved = sum(nObserved)), # agregate the chromsomes, hasn't been done yet
                keyby = .(readN, readPos1, refBases, altBases, zygosity, genotypeZygosity, genotypeCoverage)
            ]
            rp[, ":="(
                varType = parseVarType(refBases, altBases),
                baseSpec = paste(refBases, altBases, sep = ">")
            )]
            rp[, ":="(
                varType_long = parseVarType_long(varType, refBases, altBases),
                isValidated = refBases == pileupCodes$CS_MATCH | 
                              (!is.na(zygosity)         & zygosity         >= zygosityTypeIs$heterozygous) | 
                              (!is.na(genotypeZygosity) & genotypeZygosity != zygosityTypeIs$not_in_genotype)
            )]
            rp[varType <= varTypeIs$match]
        }
    )$value
}

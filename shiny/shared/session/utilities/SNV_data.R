# utility functions for loading RE site data and related
# caller must call stopSpinner()
snvData_create <- "asNeeded"

# get track data
hf3_getPileup <- function(sourceId, coord, readType){
    if(readType == "error_corrected"){
        hf3_getTrackData_bgz(sourceId, "errorCorrectedPileupBgz", coord, use_chrom = FALSE, debug = FALSE)
    } else {
        hf3_getTrackData_bgz(sourceId, "allReadsPileupBgz", coord, use_chrom = FALSE, debug = FALSE)
    }
}
hf3_getVariants <- function(sourceId, coord, readType){
    if(readType == "error_corrected"){
        hf3_getTrackData_bgz(sourceId, "errorCorrectedVariantsBgz", coord, use_chrom = FALSE, debug = FALSE)
    } else {
        hf3_getTrackData_bgz(sourceId, "allReadsVariantsBgz", coord, use_chrom = FALSE, debug = FALSE)
    }
}
hf3_getEncodings <- function(sourceId, coord, readType){
    if(readType == "error_corrected"){
        hf3_getTrackData_bgz(sourceId, "errorCorrectedEncodingsBgz", coord, use_chrom = FALSE, debug = FALSE)
    } else {
        hf3_getTrackData_bgz(sourceId, "allReadsEncodingsBgz", coord, use_chrom = FALSE, debug = FALSE)
    }
}

# get cached data, full tables
hf3_cached_create <- "asNeeded"
hf3_getVariants_cached <- function(sourceId){
    fileType <-"errorCorrectedVariantsBgz"
    sessionCache$get(
        fileType, 
        key = sourceId, 
        permanent = TRUE,
        from = "ram",
        create = hf3_cached_create,
        createFn = function(...) {
            startSpinner(session, message = "loading variants")
            dataFilePath <- getSourceFilePath(sourceId, fileType)
            d <- fread(
                cmd = paste("zcat", dataFilePath),
                col.names  =  names(hf3_bgzColumns[[fileType]]), 
                colClasses = unname(hf3_bgzColumns[[fileType]])
            )
            d[, ":="(
                n_alt_bases = nchar(alt_bases),
                vaf = count / coverage
            )]
            d[, ":="(
                is_snp = n_ref_bases == 1 & n_ref_bases == n_alt_bases,
                is_mnp = n_ref_bases > 1  & n_ref_bases == n_alt_bases,
                is_ins = n_ref_bases == 0,
                is_del = n_alt_bases == 0,
                vaf_bin = floor(vaf * 51) / 51
            )]
            d
        }  
    )$value
}
hf3_getEncodings_cached <- function(sourceId, coord){
    fileType <-"errorCorrectedEncodingsBgz"
    sessionCache$get(
        fileType, 
        key = sourceId, 
        permanent = TRUE,
        from = "ram",
        create = hf3_cached_create,
        createFn = function(...) {
            startSpinner(session, message = "loading encodings")
            dataFilePath <- getSourceFilePath(sourceId, fileType)
            d <- fread(
                cmd = paste("zcat", dataFilePath),
                col.names  =  names(hf3_bgzColumns[[fileType]]), 
                colClasses = unname(hf3_bgzColumns[[fileType]])
            )
            d[, ":="(
                n_bases = end1 - start0,
                n_bases_bin = floor((end1 - start0) / 250) * 250
            )]
            d
        }  
    )$value
}

# varTypeIs <- list(
#     snp     = 1, # single-nucleotide polymorphism, e.g. A>G
#     del1    = 2, # deletion of 1 base, e.g. A>-
#     ins1    = 3, # insertion of 1 base, e.g. ->A
#     indel0  = 4, # equal insertion + deletion of >1 base; a multi-nucleotide polymorphism, e.g. AT>CG
#     delN    = 5, # deletion of >1 base, e.g. AT>-
#     insN    = 6, # insertion of >1 base, e.g. ->AT
#     indelX  = 7, # unequal insertion + deletion of >1 base, e.g. AT>CGT, catch-all for complex events
#     # -----------
#     match   = 8, # no variant, e.g. A>A
#     clipped = 9, # 5'-clipped base
#     lowQual = 10 # masked low quality base
# )
# varTypeColors <- c(
#     snp = "blue",
#     del1 = "red",
#     ins1 = "green",
#     indel0 = "purple",
#     delN = "orange",
#     insN = "cyan",
#     indelX = "brown",
#     match = "gray",
#     clipped = "white",
#     lowQual = "white"
# )
# pileupCodes <- list(
#     CS_MATCH        = ":",
#     CLIPPED         = "!",
#     MASKED_LOW_QUAL = "q"
# )
# pileupCodeVarTypeIs <- c(
#     ":" = varTypeIs$match,
#     "!" = varTypeIs$clipped,
#     "q" = varTypeIs$lowQual
# )
# pileupCodeVals <- names(pileupCodeVarTypeIs)
# varTypes <- names(varTypeIs)
# parseVarType <- function(refBases, altBases){
#     nRefBases <- nchar(refBases)
#     nAltBases <- nchar(altBases)
#     delta <- nAltBases - nRefBases
#     ifelse(
#         refBases %in% pileupCodeVals,
#         pileupCodeVarTypeIs[refBases],
#         ifelse(
#             delta == 0,
#             ifelse(nRefBases == 1, varTypeIs$snp, varTypeIs$indel0),
#             ifelse(
#                 delta == 1,
#                 ifelse(nRefBases == 0, varTypeIs$ins1, varTypeIs$indelX),
#                 ifelse(
#                     delta == -1,
#                     ifelse(nRefBases == 1, varTypeIs$del1, varTypeIs$indelX),
#                     ifelse(
#                         delta > 1,
#                         ifelse(nRefBases == 0, varTypeIs$insN, varTypeIs$indelX),
#                         ifelse(nAltBases == 0, varTypeIs$delN, varTypeIs$indelX)
#                     )
#                 )
#             )
#         )
#     )
# }
# parseVarType_long <- function(varType, refBases, altBases){
#     nRefBases <- nchar(refBases)
#     nAltBases <- nchar(altBases)
#     delta <- nAltBases - nRefBases
#     ifelse(
#         varType %in% c(varTypeIs$snp, varTypeIs$del1, varTypeIs$ins1),
#         paste(refBases, altBases, sep = ">"),
#         ifelse(
#             varType %in% c(varTypeIs$indel0, varTypeIs$delN),
#             nRefBases,
#             ifelse(
#                 varType == varTypeIs$insN,
#                 nAltBases,
#                 pmax(nRefBases, nAltBases)
#             )
#         )
#     )
# }

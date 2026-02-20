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

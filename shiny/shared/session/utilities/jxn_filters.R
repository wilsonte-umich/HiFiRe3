# for a read, a single READ_FAILURE_FLAG is assembled sequentially as stated below
# included here for completeness; failed reads were not used during SV analysis
# however, note that a READ may be on target due to one/the first alignment 
#      but carry a distal off-target JUNCTION that is tabulated,
#      i.e., on-target filtering in the pipeline is a read, not a junction, property
hf3_readFailureBits <- list(
    Unmapped      = 1, # flag bits are enforced sequentially with short-circuiting 
    OffTarget     = 2, #     later bits are not set if earlier bits fail and are set
    SiteLookup    = 4,
    ClipTolerance = 8, # outer endpoint criteria reject entire reads, but otherwise have no bearing on junctions
    SiteDistance  = 16
)
# for an alignment, a single ALN_FAILURE_FLAG is assembled sequentially as stated below
# for a junction instance, it is a bitwise OR of two ALN_FAILURE_FLAGs from the flanking alignments
# for a final Junction, it is a bitwise OR of that value over all junction instances
hf3_alnFailureBits <- list( 
    Mapq        = 1, # flag bits are enforced sequentially with short-circuiting 
    Divergence  = 2, #     later bits are not set if earlier bits fail and are set
    FlankLen    = 4, # only set on junction reads
    BaseQual = 8  # only set on ONT junctions reads
)
# for a junction instance, a single JXN_FAILURE_FLAG is assembled sequentially as stated below
# for a final junction, it is a bitwise OR of JXN_FAILURE_FLAGs over all junction instances
# traversal failures, and only traversal failures, were not processed as candidate junctions
#     (traversal failures are not considered "real", all others are true, if artifactual, junctions in reads)
hf3_jxnFailureBits <- list( 
    Traversal = 1,  # the first group of bits are enforced sequentially with short-circuiting 
    Noncanonical   = 2,  #     later bits are not set if earlier bits fail and are set
    FoldbackInv    = 4,
    OntFollowOn    = 8, 
    HasAdapter     = 16,
    #-------------------
    SiteMatch      = 32, # the second group of bits are all checked on all junctions and OR'ed together with first group
    StemLength     = 64
    # InsertSize     = 128 # this flag if never set, kept here for legacy reasons
)
hf3_alnFailureBitsUI <- function(inputId) checkboxGroupInput(
    inputId,
    "Aln Flag Bits",
    inline = TRUE,
    choices = hf3_alnFailureBits
)
hf3_jxnFailureBitsUI <- function(inputId) checkboxGroupInput(
    inputId,
    "Jxn Flag Bits",
    inline = TRUE,
    choices = hf3_jxnFailureBits
)
hf3_flagFilterModeUI <- function(inputId, label, selected) radioButtons(
    inputId, 
    label, 
    inline = TRUE,
    choices = c("show_all", "any_pass", "all_pass", "any_fail", "all_fail"), 
    selected = selected
)

# filter junction clusters based on track or other settings
hf3_JxnFilterDefaults <- list( 
    Min_Breakpoint_Coverage = 0,
    One_End_In_Genome       = "any",
    # Enforce_1N_to_2N_Size   = FALSE,
    Allow_Excluded_Regions  = FALSE
)
hf3_anyReadFlagPassed <- function(readFlags, flagBits){
    for (readFlag in readFlags) {
        if(bitwAnd(readFlag, flagBits) == 0) return(TRUE)
    }
    FALSE
}
hf3_anyReadFlagFailed <- function(readFlags, flagBits){
    for (readFlag in readFlags) {
        if(bitwAnd(readFlag, flagBits) == flagBits) return(TRUE)
    }
    FALSE
}
hf3_allReadFlagsPassed <- function(readFlags, flagBits){
    for (readFlag in readFlags) {
        if(bitwAnd(readFlag, flagBits) != 0) return(FALSE)
    }
    TRUE
}
hf3_allReadFlagsFailed <- function(readFlags, flagBits){
    for (readFlag in readFlags) {
        if(bitwAnd(readFlag, flagBits) != flagBits) return(FALSE)
    }
    TRUE
}
hf3_enforceFlagFilters <- function(failure_flags, failure_bits, filter_mode){
    readFlagsList <- sapply(failure_flags, function(x) as.integer(strsplit(sub(",", "", x), ",")[[1]]))
    flagBits  <- sum(as.integer(failure_bits))
    fn <- switch(
        filter_mode,
        any_pass = hf3_anyReadFlagPassed,
        any_fail = hf3_anyReadFlagFailed,
        all_pass = hf3_allReadFlagsPassed,
        all_fail = hf3_allReadFlagsFailed
    )
    sapply(readFlagsList, fn, flagBits)
}

hf3_applyJunctionFilters <- function(jxns, settings, input){

    startSpinner(session, message = paste("filtering junctions"))

    # enforce filters from setting, generally less commonly used
    filters <- if(is.null(settings$Junction_Filters)) list() else settings$Junction_Filters()
    filters <- lapply(names(hf3_JxnFilterDefaults), function(filter){
        if(is.null(filters[[filter]])) hf3_JxnFilterDefaults[[filter]]
        else if(!is.null(filters[[filter]]$selected)) filters[[filter]]$selected else filters[[filter]]$value
    })
    names(filters) <- names(hf3_JxnFilterDefaults)
    if(filters$Min_Breakpoint_Coverage > 1) jxns <- jxns[
        bkpt_coverage_1 >= filters$Min_Breakpoint_Coverage | 
        bkpt_coverage_2 >= filters$Min_Breakpoint_Coverage
    ]
    oneEndInGenome <- trimws(filters$One_End_In_Genome)
    if(oneEndInGenome != "any" && oneEndInGenome != ""){
        oneEndInGenome <- paste0("_", oneEndInGenome)
        jxns <- jxns[
            endsWith(hf3_getChromNames(sourceId, chrom_index1_1), oneEndInGenome) |
            endsWith(hf3_getChromNames(sourceId, chrom_index1_2), oneEndInGenome)
        ]
    }
    if(!filters$Allow_Excluded_Regions) jxns <- jxns[
        is_excluded_1 == 0 & is_excluded_2 == 0
    ]

    # enforce filters from inputs specifying failure flag bits
    if(input$alnFlagFilterMode != "show_all" && !is.null(input$alnFailureBits)){
        keep <- hf3_enforceFlagFilters(jxns$aln_failure_flags, input$alnFailureBits, input$alnFlagFilterMode)
        jxns <- jxns[keep]
    }
    if(input$jxnFlagFilterMode != "show_all" && !is.null(input$jxnFailureBits)){
        keep <- hf3_enforceFlagFilters(jxns$jxn_failure_flags, input$jxnFailureBits, input$jxnFlagFilterMode)
        jxns <- jxns[keep]
    }
    jxns
}

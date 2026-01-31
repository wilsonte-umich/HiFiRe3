# for a read, a single READ_FAILURE_FLAG was assembled sequentially as stated below
# included here for completeness; failed reads were not used during SV analysis
# however, a READ may be on target due to one/the first alignment 
#      but carry a distal off-target JUNCTION that is tabulated,
#      i.e., on-target filtering in the pipeline is a read, not a junction-level, property
# read failure flag bits were enforced sequentially with short-circuiting 
#     i.e., later bits will never be set if earlier bits failed and were set
# outer endpoint criteria reject entire reads, but otherwise have no bearing on junctions
hf3_readFailureBits <- list(
    Unmapped      = 1,  # minimap2 failed to align the read to the reference genome
    OffTarget     = 2,  # the first (ONT) or all (other platforms) alignments were off-target
    SiteLookup    = 4,  # one or more alignment nodes failed lookup in the closest RE site list (not common)
    ClipTolerance = 8,  # one or both read outer clips exceeded the maximum allowed clip tolerance
    SiteDistance  = 16  # one or both read outer nodes were too far from a restriction enzyme site (ligFree only)
)
# for an alignment, a single ALN_FAILURE_FLAG was assembled sequentially as stated below
# for a junction instance, the value is a bitwise OR of two ALN_FAILURE_FLAGs from the flanking alignments
# for a final Junction, the value is a bitwise OR of the values over all junction instances
# alignment failure flag bits were enforced sequentially with short-circuiting 
#     i.e., later bits will never be set if earlier bits failed and were set
hf3_alnFailureBits <- list( 
    Mapq        = 1, # the minimap2 mapping quality (MAPQ) was below the minimum threshold
    Divergence  = 2, # the alignment divergence tag (de:f:) value was above the maximum threshold
    FlankLen    = 4, # the alignment was shorter than the allowed flank length
    BaseQual    = 8  # the average base quality was below the minimum threshold
)
# for a junction instance, a single JXN_FAILURE_FLAG was assembled sequentially as stated below
# for a final junction, the value is a bitwise OR of the JXN_FAILURE_FLAGs over all junction instances
# traversal failures, and only traversal failures, were not processed as candidate junctions
#     traversal failures are not considered "real", others are true, if artifactual, junctions in reads
# the first group of junction failure flag bits were enforced sequentially with short-circuiting 
#     i.e., later bits will never be set if earlier bits failed and were set
# the second group (SiteMatch and StemLength) were both checked on all junctions and OR'ed together with first group
hf3_jxnFailureBits <- list( 
    # Traversal    = 1,  # the number of base traversed on reference and read between breakpoints was too similar
    # Noncanonical = 2,  # one alignment was to a non-nuclear chromosome, e.g., an unplaced contig
    FoldbackInv    = 4,  # the junction was classified as a foldback inversion consistent with sequence two strands of one duplex
    OntFollowOn    = 8,  # the junction had inserted bases with a low quality stretch as seen in ONT follow-on artifacts; ONT only
    "HasAdapter ||"     = 16, # the junction had inserted bases that aligned to known adapter sequences
    #-------------------
    SiteMatch      = 32, # one or both junction breakpoints were too close to a filtering RE site
    StemLength     = 64  # both junction breakpoints were too far from the read end, i.e., failed the <1N filter
    # InsertSize   = 128 # the read carring the junction failed the 1N to 2N size filter; never set, listed for legacy reasons
)


# construct the filtering UI
hf3_flagFilterModeUI <- function(inputId, passFail, label = NULL) radioButtons(
    inputId, 
    label, 
    inline = TRUE,
    choices = if(passFail == "pass") c("show_all", "any_pass", "all_pass")
                               else  c("show_all", "any_fail", "all_fail"), 
    selected = "show_all"
)
hf3_alnFailureBitsUI <- function(inputId, label = FALSE) checkboxGroupInput(
    inputId,
    if(label) "Aln Failure Flag (<=1 one set, in order)" else NULL,
    inline = TRUE,
    choices = hf3_alnFailureBits
)
hf3_jxnFailureBitsUI <- function(inputId, label = FALSE) checkboxGroupInput(
    inputId,
    if(label) "Jxn Failure Flag (group 1: <=1 one set, in order; group 2: both can be set)" else NULL,
    inline = TRUE,
    choices = hf3_jxnFailureBits
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
    startSpinner(session, message = paste("filtering junctions."))
    oneEndInGenome <- trimws(filters$One_End_In_Genome)
    if(length(jxns) > 0 && oneEndInGenome != "any" && oneEndInGenome != ""){
        oneEndInGenome <- paste0("_", oneEndInGenome)
        jxns <- jxns[
            endsWith(hf3_getChromNames(sourceId, chrom_index1_1), oneEndInGenome) |
            endsWith(hf3_getChromNames(sourceId, chrom_index1_2), oneEndInGenome)
        ]
    }
    startSpinner(session, message = paste("filtering junctions..."))
    if(length(jxns) > 0 && !filters$Allow_Excluded_Regions) jxns <- jxns[
        is_excluded_1 == 0 & is_excluded_2 == 0
    ]

    # enforce filters from inputs specifying failure flag bits
    startSpinner(session, message = paste("filtering junctions...."))
    if(length(jxns) > 0 && input$alnFlagPassMode != "show_all" && !is.null(input$alnFlagPassBits)){
        keep <- hf3_enforceFlagFilters(jxns$aln_failure_flags, input$alnFlagPassBits, input$alnFlagPassMode)
        jxns <- jxns[keep]
    }
    startSpinner(session, message = paste("filtering junctions....."))
    if(length(jxns) > 0 && input$alnFlagFailMode != "show_all" && !is.null(input$alnFlagFailBits)){
        keep <- hf3_enforceFlagFilters(jxns$aln_failure_flags, input$alnFlagFailBits, input$alnFlagFailMode)
        jxns <- jxns[keep]
    }
    startSpinner(session, message = paste("filtering junctions......"))
    if(length(jxns) > 0 && input$jxnFlagPassMode != "show_all" && !is.null(input$jxnFlagPassBits)){
        keep <- hf3_enforceFlagFilters(jxns$jxn_failure_flags, input$jxnFlagPassBits, input$jxnFlagPassMode)
        jxns <- jxns[keep]
    }
    startSpinner(session, message = paste("filtering junctions......."))
    if(length(jxns) > 0 && input$jxnFlagFailMode != "show_all" && !is.null(input$jxnFlagFailBits)){
        keep <- hf3_enforceFlagFilters(jxns$jxn_failure_flags, input$jxnFlagFailBits, input$jxnFlagFailMode)
        jxns <- jxns[keep]
    }
    jxns
}

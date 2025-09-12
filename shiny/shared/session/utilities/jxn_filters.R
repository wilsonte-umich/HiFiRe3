# filter junction clusters based on track or other settings
af_JxnFilterDefaults <- list( 
    # Min_SV_Size = 0,
    # Max_SV_Size = 0,
    # Min_N_Observed = 0,
    # Max_N_Observed = 0,
    Min_Breakpoint_Coverage = 0,
    Min_Map_Quality = 0,
    Max_Frac_Base_Differences = 1,
    Min_Site_Distance = 0,
    # Is_Expected_SV = c("expected","unexpected"),
    One_End_In_Genome = "any",
    Is_On_Target = c("on_target","off_target"),
    Allow_Excluded_Regions = FALSE,
    Allow_Alternative_Alignments = FALSE
    # SV_Type = c("deletion","duplication","inversion","translocation"),
    # Mixed_Genome = c("intragenome","intergenome")
)
af_applyJunctionFilters <- function(jxns, settings){

    startSpinner(session, message = paste("filtering junctions"))

    filters <- if(is.null(settings$Junction_Filters)) list() else settings$Junction_Filters()
    filters <- lapply(names(af_JxnFilterDefaults), function(filter){
        if(is.null(filters[[filter]])) af_JxnFilterDefaults[[filter]]
        else if(!is.null(filters[[filter]]$selected)) filters[[filter]]$selected else filters[[filter]]$value
    })
    names(filters) <- names(af_JxnFilterDefaults)

    # if(filters$Min_SV_Size > 1) jxns <- jxns[svSize >= filters$Min_SV_Size]
    # if(filters$Max_SV_Size > 0) jxns <- jxns[svSize <= filters$Max_SV_Size]

    # if(filters$Min_N_Observed > 1) jxns <- jxns[nObserved >= filters$Min_N_Observed]
    # if(filters$Max_N_Observed > 0) jxns <- jxns[nObserved <= filters$Max_N_Observed]
    if(filters$Min_Breakpoint_Coverage > 1) jxns <- jxns[
        bkptCoverage_1 >= filters$Min_Breakpoint_Coverage | 
        bkptCoverage_2 >= filters$Min_Breakpoint_Coverage
    ]
    if(filters$Min_Map_Quality > 0) jxns <- jxns[
        mapQ >= filters$Min_Map_Quality
    ]
    if(filters$Max_Frac_Base_Differences < 1) jxns <- jxns[
        deTag<= filters$Max_Frac_Base_Differences
    ]

    if(filters$Min_Site_Distance > 0 && jxns[, any(siteDist !=0 )]) jxns <- jxns[
        siteDist >= filters$Min_Site_Distance
    ]

    # if(length(filters$Is_Expected_SV) == 1) 
    #     jxns <- jxns[readHasSv == (filters$Is_Expected_SV == "unexpected")]
    if(length(filters$Is_On_Target) == 1){
        jxnIsOnTarget <- jxns[, target1 != "*" | target2 != "*"]
        if (any(jxnIsOnTarget)) jxns <- jxns[
            jxnIsOnTarget == (filters$Is_On_Target == "on_target")
        ]
    }
    oneEndInGenome <- trimws(filters$One_End_In_Genome)
    if(oneEndInGenome != "any" && oneEndInGenome != ""){
        oneEndInGenome <- paste0("_", oneEndInGenome)
        jxns <- jxns[
            endsWith(af_getChromNames(sourceId, chromIndex1_1), oneEndInGenome) |
            endsWith(af_getChromNames(sourceId, chromIndex1_2), oneEndInGenome)
        ]
    }
    if(!filters$Allow_Excluded_Regions) jxns <- jxns[
        isExcluded_1 == 0 & isExcluded_2 == 0
    ]
    if(!filters$Allow_Alternative_Alignments) jxns <- jxns[
        hasAltAlignment == 0
    ]
    # if(length(filters$SV_Type) > 0 && length(filters$SV_Type) < 4) 
    #     jxns <- jxns[jxnType %in% af_junctions$getJxnTypeIndices(filters$SV_Type)]
    # if(length(filters$Mixed_Genome) == 1) 
    #     jxns <- jxns[isIntergenome == (filters$Mixed_Genome == "intergenome")]
    # req(nrow(jxns) > 0)
    jxns
}

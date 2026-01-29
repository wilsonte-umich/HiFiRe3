#----------------------------------------------------------------------
# analyze_SVs_reads trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_analyze_SVs_readsTrack <- function(...) {
    list(
        click = FALSE,
        hover = FALSE,
        brush = FALSE,
        items = TRUE,
        navigation = FALSE, 
        expand = FALSE, # TODO: enable expand on plotted spans
        expand2 = FALSE
    )
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.analyze_SVs_readsTrack <- function(...) showTrackSourcesDialog(...)

# build method for the S3 class; REQUIRED
build.analyze_SVs_readsTrack <- function(track, reference, coord, layout){

    # check and enforce max width
    maxWidth_sites <- track$settings$get("SV_Reads","Sites_Max_Width")
    if(isTruthy(coord$width > maxWidth_sites))
        return(trackInfo(track, coord, layout, "window too wide to plot SV sites and fragments"))
    maxWidth_alignments <- track$settings$get("SV_Reads","Alignments_Max_Width")

    # collect data sources
    selectedSourceIds <- names(track$settings$items())
    nSources <- length(selectedSourceIds)
    req(nSources > 0)
    if(nSources > 1)
        return(trackInfo(track, coord, layout, "analyze_SVs_reads only plots one data source (sample) at a time"))
    sourceId <- selectedSourceIds[1]
    projectName <- strsplit(getSourceFilePackageName(sourceId), "\\.")[[1]][1]

    # set plot configuration
    # maxCoverage <- track$settings$get("SV_Reads","Fragments_Max_Coverage")
    padding <- padding(track, layout)

    # get RE sites in window
    startSpinner(session, message = "building SVs reads.")
    sites <- hf3_getSites_padded(sourceId, coord)
    if(nrow(sites) == 0) return(trackInfo(track, coord, layout, "no RE sites in window"))

    # get alignments in window
    startSpinner(session, message = "building SVs reads..")
    tryCatch({ if(coord$width <= maxWidth_alignments) {
        alignments <- hf3_getAlignments(sourceId, coord)
        if(nrow(alignments) == 0) return(trackInfo(track, coord, layout, "no alignments to plot"))

        # stack the boxes
        # top strand reads
        startSpinner(session, message = "building SVs reads...")
        siteY_F <- sites[, .(refPos1 = sitePos1, ybottom = 0)]
        blocks_F <- alignments[
            strand_index0 == hf3_strands$top,
            {
                is_site_overrun <- ref_proj3 > 0 & ref_pos3 - ref_proj3 >= 1 & ref_pos3 - ref_proj3 <= 5
                xright_aln <- ifelse(is_site_overrun, ref_proj3, ref_pos3)
                xright_prj <- ifelse(ref_proj3 > 0 & xright_aln < ref_proj3, ref_proj3, NA)
                .(
                    jxn_types  = jxn_types,
                    n_observed = n_observed,
                    xleft      = ref_pos5, # matching a RE site, except for distal SV alignments on top strand
                    xright_aln = xright_aln,
                    xright_prj = xright_prj,
                    xMax       = pmax(xright_aln, xright_prj, na.rm = TRUE)
                )
            }
        ][
            order(xleft, -xMax, -xright_aln)
        ]
        for(i in seq_len(nrow(blocks_F))){
            block <- blocks_F[i]
            siteLeft <- if(block$xleft %in% siteY_F$refPos1) block$xleft
                        else siteY_F$refPos1[max(which(siteY_F$refPos1 <= block$xleft))]
            blocks_F[i, ybottom := siteY_F[refPos1 == siteLeft, ybottom]]
            siteY_F[refPos1 %between% c(siteLeft, block$xMax), ybottom := ybottom + block$n_observed]
        }

        # bottom strand reads
        siteY_R <- sites[, .(refPos1 = sitePos1 - 1, ytop = 0)]
        startSpinner(session, message = "building SVs reads....")
        blocks_R <- alignments[
            strand_index0 == hf3_strands$bottom, {
                is_site_overrun <- ref_proj3 > 0 & ref_proj3 - ref_pos3 >= 1 & ref_proj3 - ref_pos3 <= 5
                xleft_aln <- ifelse(is_site_overrun, ref_proj3, ref_pos3)
                xleft_prj <- ifelse(ref_proj3 > 0 & xleft_aln > ref_proj3, ref_proj3, NA)
                .(
                    jxn_types  = jxn_types,
                    n_observed = n_observed,
                    xright    = ref_pos5,
                    xleft_aln = ref_pos3,
                    xleft_prj = xleft_prj,
                    xMin      = pmin(xleft_aln, xleft_prj, na.rm = TRUE)
                )
            }
        ][
            order(-xright, xMin, xleft_aln)
        ]
        for(i in seq_len(nrow(blocks_R))){
            block <- blocks_R[i]
            siteRight <- if(block$xright %in% siteY_R$refPos1) block$xright
                        else siteY_R$refPos1[max(which(siteY_R$refPos1 >= block$xright))]
            blocks_R[i, ytop := siteY_R[refPos1 == siteRight, ytop]]
            siteY_R[refPos1 %between% c(block$xMin, siteRight), ytop := ytop - block$n_observed]
        }
        ylim <- c(siteY_R[, min(ytop)], siteY_F[, max(ybottom)])
        yaxt <- NULL
        alignmentPixelHeight <- track$settings$get("SV_Reads","Alignment_Pixel_Height")
        height <- (diff(ylim) + 2) * alignmentPixelHeight / layout$dpi + padding$total
    } else {
        ylim <- c(0, 1)
        yaxt <- "n"
        height <- height(track, 0.25) + padding$total
    }}, error = function(e){
        print(e)
        return(trackInfo(track, coord, layout, "error stacking fragments, probably due to too few sites, etc."))
    })
    startSpinner(session, message = "building SVs reads....")

    # initialize the plot frame
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim,  ylab = projectName, yaxt = yaxt,
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i"

        # plot the alignment rectangles
        if(coord$width <= maxWidth_alignments){
            rect(
                xleft   = blocks_F$xleft, 
                ybottom = blocks_F$ybottom, 
                xright  = blocks_F$xright_prj, 
                ytop    = blocks_F$ybottom + blocks_F$n_observed, 
                col     = hf3_alignments$typeToColor$projection, 
                border  = NA
            )
            rect(
                xleft   = blocks_F$xleft, 
                ybottom = blocks_F$ybottom, 
                xright  = blocks_F$xright_aln, 
                ytop    = blocks_F$ybottom + blocks_F$n_observed, 
                col     = hf3_junctions$getColorsFromBits_no_intergenomic(blocks_F$jxn_types), 
                border  = NA
            )
            rect(
                xleft   = blocks_R$xleft_prj, 
                ybottom = blocks_R$ytop - blocks_R$n_observed, 
                xright  = blocks_R$xright, 
                ytop    = blocks_R$ytop, 
                col     = hf3_alignments$typeToColor$projection, 
                border  = NA
            )
            rect(
                xleft   = blocks_R$xleft_aln, 
                ybottom = blocks_R$ytop - blocks_R$n_observed, 
                xright  = blocks_R$xright, 
                ytop    = blocks_R$ytop,
                col     = hf3_junctions$getColorsFromBits_no_intergenomic(blocks_R$jxn_types), 
                border  = NA
            )            
        }

        # plot the RE sites
        if(nrow(sites) > 0) {
            siteLineWidth <- track$settings$get("SV_Reads","Site_Line_Width")
            for(i in 1:nrow(sites)){
                lines(
                    rep(sites[i, sitePos1 - 0.5], 2), 
                    ylim, 
                    lwd = siteLineWidth
                )
            }
        }
    })

    list(
        ylim  = ylim,
        mai   = mai,
        image = image
    )
}

# plot interaction methods for the S3 class
# called by trackBrowser if track$click, $hover, or $brush is TRUE, above
# regionI indicates the region plot the user interacted with, the value must be passed to app$browser$jumpToCoordinates, etc.
click.analyze_SVs_readsTrack <- function(track, click, regionI){
    # custom actions, use str(click) to explore
}
hover.analyze_SVs_readsTrack <- function(track, hover, regionI){
    # custom actions, use str(hover) to explore
}

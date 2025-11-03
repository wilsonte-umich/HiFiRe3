#----------------------------------------------------------------------
# analyze_SVs_arcs trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------

# constructor for the S3 class
new_analyze_SVs_arcsTrack <- function(...) {
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
items.analyze_SVs_arcsTrack <- function(...) showTrackSourcesDialog(...)

# build method for the S3 class; REQUIRED
build.analyze_SVs_arcsTrack <- function(track, reference, coord, layout){

    # collect data sources
    selectedSourceIds <- names(track$settings$items())
    nSources <- length(selectedSourceIds)
    req(nSources > 0)
    if(nSources > 1)
        return(trackInfo(track, coord, layout, "analyze_SVs_arcs only plots one data source (sample) at a time"))
    sourceId <- selectedSourceIds[1]

    # get junctions with >=1 node in window
    startSpinner(session, message = "building SVs arcs.")
    junctions <- af_getFilteredJunctions(sourceId, coord)[isTruthy(chromIndex1_1) & isTruthy(chromIndex1_2)]
    startSpinner(session, message = "building SVs arcs..")

    # set plot configuration
    padding <- padding(track, layout)
    height <- height(track, 0.25) + padding$total # or set a known, fixed height in inches
    Upside_Down <- track$settings$get("SV_Arcs","Upside_Down")
    # Color_By  <- track$settings$get("SV_Arcs","Color_By")
    Line_Width  <- track$settings$get("SV_Arcs","Arc_Line_Width")
    Max_Y_BP    <- track$settings$get("SV_Arcs","Max_Y_BP")
    ymax <- switch(
        Max_Y_BP,
        window_width = diff(coord$range) * 0.51,
        max_sv_width = {
            intraChrom <- junctions[jxnType != af_junctions$typeToIndex$translocation]
            if(nrow(intraChrom) == 0) diff(coord$range) 
            else intraChrom[, max(abs(refPos1_2 - refPos1_1))]
        } * 0.51,
        fixed = as.integer(track$settings$get("SV_Arcs","Fixed_Y_BP"))
    )
    ylim <- c(0, ymax)

    # initialize the plot frame
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim, ylab = ylab(track, ""), yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i" 

        # plot the arcs
        x <- seq(0, pi, length.out = 25)
        if(nrow(junctions) == 0) trackNoData(coord, ylim, "no matching junctions in window") else {
            junctions <- junctions[sample.int(.N, replace = FALSE)]
            junctions[, randomI := 1:.N]
            transOffset <- coord$width / 4
            junctions[, {
                if(jxnType == af_junctions$typeToIndex$translocation){
                    if(chromIndex1_1 == af_getChromIndex(sourceId, coord$chromosome)){
                        pos <- refPos1_1
                        dir <- 1
                    } else {
                        pos <- refPos1_2
                        dir <- -1
                    }
                    lines(
                        x = c(pos, pos + dir * transOffset), 
                        y = ylim,
                        col = af_junctions$getIndexColor(jxnType),
                        lwd = Line_Width
                    )
                } else {
                    pos <- range(refPos1_1, refPos1_2)
                    halfsize <- (pos[2] - pos[1] + 1) / 2
                    center <- pos[1] + halfsize
                    y <- sin(x) * halfsize
                    if(Upside_Down) y <- ymax - y
                    lines(
                        x = cos(x) * halfsize + center, 
                        y = y,
                        col = af_junctions$getIndexColor(jxnType),
                        lwd = Line_Width
                    )
                }
            }, by = .(randomI)]
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
click.analyze_SVs_arcsTrack <- function(track, click, regionI){
    # custom actions, use str(click) to explore
}
hover.analyze_SVs_arcsTrack <- function(track, hover, regionI){
    # custom actions, use str(hover) to explore
}

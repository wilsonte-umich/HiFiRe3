#----------------------------------------------------------------------
# analyze_read_pileup trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
pileupBaseColors <- list( # generally follow IGV base color conventions
    M = rgb(0.75, 0.75, 0.75),   # any base match = light grey
    A = rgb(0,    0.8,    0),    # green  
    C = rgb(0,    0,    1),    # blue
    G = rgb(0.82, 0.43, 0),      # orange
    T = rgb(1,    0,    0),    # red
    N = rgb(0.5, 0.5, 0.5),      # N, treated as M
    D = rgb(0.1, 0.1, 0.1),      # deleted/missing = black
    I = rgb(0.75,   0,    0.75)  # insertion = purple
)

# constructor for the S3 class; REQUIRED
new_analyze_read_pileupTrack <- function(trackId) {
    list(
        click = FALSE, # whether the track type has `click`, `hover`, and/or `items` methods
        hover = FALSE,
        brush = FALSE,
        items = TRUE,
        navigation = FALSE, # whether the track offers a custom, additional row of within-track navigation inputs
        expand = FALSE,
        expand2 = FALSE
    )
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.analyze_read_pileupTrack <- function(...) showTrackSourcesDialog(...)

# build method for the S3 class; REQUIRED
build.analyze_read_pileupTrack <- function(track, reference, coord, layout){

    # collect data sources
    selectedSourceIds <- names(track$settings$items())
    nSources <- length(selectedSourceIds)
    req(nSources > 0)
    if(nSources > 1)
        return(trackInfo(track, coord, layout, "read pileup only plots one data source (sample) at a time"))
    sourceId <- selectedSourceIds[1]

    # load the pileup for the sample in the window
    Read_Type <- track$settings$get("Pileup", "Read_Type")
    pileup <- hf3_getPileup(sourceId, coord, Read_Type)
    if(nrow(pileup) == 0)
        return(trackInfo(track, coord, layout, "no read pileup data in window for this sample"))

    # select the pileup counts based on the allowed filter
    Include_Unallowed <- track$settings$get("Pileup", "Include_Unallowed")
    if(Include_Unallowed) {
        in_cols  <- c("M","A","C","G","T","D","N","I")
        out_cols <- in_cols
    } else {
        in_cols  <- c("M","A","C","G","T","D_allowed","I_allowed")
        out_cols <- c("M","A","C","G","T","D",        "I")
    }
    pos_cols <- c("start", "end")
    pileup <- pileup[, .SD, .SDcols = c(pos_cols, in_cols)]
    setnames(pileup, c(pos_cols, out_cols))
    refCols <- out_cols[out_cols != "I"]
    altCols <- refCols[refCols != "M"]

    # adjust the pileup blocks to fit the window
    pileup[, start := start + 1L]
    pileup <- pileup[end >= coord$start & start <= coord$end] 
    pileup[, ":="(
        start = pmax(start, as.integer(coord$start)), 
        end   = pmin(end,   as.integer(coord$end))
    )]

    ############################
    variants <- hf3_getVariants(sourceId, coord, Read_Type)
    dprint(pileup[, -1])
    if (nrow(variants) > 0) {
        variants[, start0 := start0 + 1L]
        variants[, alt_bases := nchar(alt_bases)]
        dprint(variants[, -1])
    }

    # set the layout
    padding <- padding(track, layout)
    height <- height(track, 1) + padding$total # or set a known, fixed height in inches
    ref_coverage <- pileup[, rowSums(.SD, na.rm = TRUE), .SDcols = refCols]
    ymin <- pileup[, max(I)]
    ymax <- max(ref_coverage)
    ylim <- c(-ymin, ymax)

    # ensure that all pileup blocks are always visible
    plotWidthPixels <- layout$plotWidth * layout$dpi
    basesPerPixel <- coord$width / plotWidthPixels
    minBarPixels <- 3
    barHalfWidthBases <- (basesPerPixel * minBarPixels) / 2
    pileup[, pixelWidth := (end - start + 1) / basesPerPixel] 
    pileup[pixelWidth < minBarPixels, ":="(
        start = floor(  (start + end) / 2 - barHalfWidthBases + 0.5),
        end   = ceiling((start + end) / 2 + barHalfWidthBases - 0.5)
    )]

    # use the mdiTrackImage helper function to create the track image
    mai <- NULL
    image <- mdiTrackImage(layout, height, function(...){
        mai <<- setMdiTrackMai(layout, padding, mar = list(top = 0, bottom = 0))
        plot(0, 0, type = "n", bty = "n",
            xlim = coord$range, xlab = "", xaxt = "n", # nearly always set `xlim`` to `coord$range`
            ylim = ylim, ylab = ylab(track, "Read Count"), #yaxt = "n",
            xaxs = "i", yaxs = "i") # always set `xaxs` and `yaxs` to "i" 

        rect(
            xleft   = pileup$start - 0.5, 
            xright  = pileup$end   + 0.5, 
            ybottom = 0, 
            ytop    = ref_coverage, 
            col     = pileupBaseColors$M, 
            border  = NA
        )
        ybottom <- rep(0, nrow(pileup))
        for(col in altCols) {
            rect(
                xleft   = pileup$start - 0.5, 
                xright  = pileup$end   + 0.5, 
                ybottom = ybottom, 
                ytop    = ybottom + pileup[[col]], 
                col     = pileupBaseColors[[col]], 
                border  = NA
            )
            ybottom <- ybottom + pileup[[col]]
        }
        if(any(pileup$I > 0)) {
            rect(
                xleft   = pileup$start, 
                xright  = pileup$end + 1,
                ybottom = -pileup$I, 
                ytop    = 0, 
                col     = pileupBaseColors$I, 
                border  = NA
            )
        }

        # add a legend
        trackLegend(
            track, coord, ylim, 
            legend = out_cols, 
            col = unlist(pileupBaseColors[out_cols]), 
            pch = 19, cex = 1.1
        )
    })

    # return the track's magick image and associated metadata
    list(ylim = ylim, mai = mai, image = image)
}

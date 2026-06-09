#----------------------------------------------------------------------
# analyze_SVs_alns trackBrowser track (i.e., a browserTrack)
#----------------------------------------------------------------------
analyze_SVs_alnsTrackBuffer <- list()

# constructor for the S3 class
new_analyze_SVs_alnsTrack <- function(...) {
    list(
        click = FALSE,
        hover = FALSE,
        brush = FALSE,
        items = TRUE,
        navigation = FALSE, 
        expand = FALSE,
        expand2 = FALSE
    )
}
 
# method for the S3 class to show a relevant trackItemsDialog or trackSamplesDialog
# used when a track can take a list of items to be plotted together and the item list icon is clicked
items.analyze_SVs_alnsTrack <- function(...) showTrackSourcesDialog(keyColumn = "Data_Package", ...)

# build method for the S3 class; REQUIRED
build.analyze_SVs_alnsTrack <- function(track, reference, coord, layout){
    req(coord, coord$chromosome)
    Max_Window_Width   <- getTrackSetting(track, "Alignments", "Max_Window_Width", 1000000)
    if(coord$width > Max_Window_Width) return(trackInfo(track, coord, layout, ": window too wide to plot alignments"))
    padding <- 100000 # avoid edge effects in alignment recovery
    dataFn <- function(track, reference, coord, sampleName, source){
        coord$start <- coord$start - padding
        coord$end   <- coord$end   + padding
        d <- hf3_getAlignments(source$Source_ID, coord)
# Classes 'data.table' and 'data.frame':  329 obs. of  9 variables:
#  $ chrom_index1 : int  1 1 1 1 1 1 1 1 1 1 ...
#  $ ref_pos5     : int  10057 10468 10476 14692 15804 20978 21095 21471 21960 24136 ...
#  $ ref_pos3     : int  10468 10001 25697 23964 27992 32986 30545 32909 31459 13867 ...
#  $ ref_proj3    : int  10468 10001 25697 23964 27992 32986 30545 32909 31459 13867 ...
#  $ strand_index0: int  0 1 0 0 0 0 0 0 0 1 ...
#  $ jxn_types    : int  8 8 0 0 0 0 0 0 0 0 ...
#  $ n_jxns       : int  6 2 0 0 0 0 0 0 0 0 ...
#  $ aln_i        : int  0 1 0 0 0 0 0 0 0 0 ...
#  $ n_observed   : int  1 1 1 1 1 1 1 1 1 1 ...
        col <- hf3_junctions$getColorsFromBits_no_intergenomic(d$jxn_types)
        x <- d[, .(
            chrom  = coord$chromosome,
            start  = pmin(ref_pos5, ref_pos3),
            end    = pmax(ref_pos5, ref_pos3),
            name   = ".",
            score  = (abs(ref_pos3 - ref_pos5) + 1) / 1000,
            ref_pos5  = ref_pos5,
            ref_pos3  = ref_pos3,
            jxnCol_5 = ifelse(jxn_types > 0 & aln_i > 0,      col, NA_integer_),
            jxnCol_3 = ifelse(jxn_types > 0 & aln_i < n_jxns, col, NA_integer_)
        )]
        x
    }
    build.genome_spans_track(
        track, reference, coord, layout, 
        dataFn = dataFn, 
        trackBuffer = NULL,
        scoreLabel = "Aln. Length (kb)",
        overplotSpansFn = function(d){
            Point_Size <- getTrackSetting(track, "Scores", "Point_Size", 1000000)
            points(d$ref_pos5, d$y, pch = 19, cex = Point_Size, col = d$jxnCol_5)
            points(d$ref_pos3, d$y, pch = 19, cex = Point_Size, col = d$jxnCol_3)
        },
        legendNames = function(track){
            sapply(track$settings$items(), function(dt){
                x <- sub(".analyze.SVs", "", dt$Data_Package[1])
                x <- sub(".compare.SVs", "", x)
                     sub(".merge.SVs",   "", x)
            })
        }
    )
}

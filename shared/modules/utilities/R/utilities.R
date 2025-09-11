
#----------------------------------------------------------------------
# miscellaneous generic support tools and functions
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# data frame tools
#----------------------------------------------------------------------

# reduce a data frame to unique rows based on queried columns
uniqueRows <- function(df, cols) df[!duplicated(df[cols]), ] 

#----------------------------------------------------------------
# string functions
#----------------------------------------------------------------

# convert first character in word to upper case
ucFirst <- function(y) { 
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1, 1)), substring(c, 2),
      sep = "", collapse = " ")
}

#----------------------------------------------------------------
# miscellaneous functions
#----------------------------------------------------------------

# shortcut for the opposite of %in%
`%notin%` <- Negate(`%in%`)

# remove the first element of an vector
# return the altered vector (NOT the shifted value)
shift <- function(v) if(length(v) <= 1) c() else v[2:length(v)]

#----------------------------------------------------------------
# generic numeric and distribution support functions
#----------------------------------------------------------------
peakValue <- function(x){ # find the peak, i.e., mode of a set of values
    x <- x[!is.na(x)]
    if(length(x) == 0) return(NA)
    d <- density(x)
    d$x[which.max(d$y)]
}
excludeOutliers <- function(v, low = 0.025, high = 0.975, min = -Inf){ # could be more sophisticated, but these are heuristics
    v <- v[!is.na(v)]
    q <- quantile(v, c(low, high))
    v[v >= q[1] & v <= q[2] & v >= min]
}

#----------------------------------------------------------------
# plot support
#----------------------------------------------------------------
addAlphaToColor <- function(color, alpha) { # expects a single color and a vector of alphas[0:1]
    if(alpha < 0 || alpha >= 1) return(color)
    rgb_ <- col2rgb(color)
    rgb(
        rgb_[1], rgb_[2], rgb_[3],
        max = 255,
        alpha = pmin(1, pmax(0, alpha)) * 255
    )
}
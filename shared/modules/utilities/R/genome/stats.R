# genome and related stats functions

# N50 of a set of contigs, reads, etc. (or any other Nxx "quantile")
N50 <- function(lengths, fraction = 0.5) {
    lengths <- rev(sort(lengths))
    lengths[cumsum(as.numeric(lengths)) >= sum(lengths) * (1 - fraction)][1]
}
N50_weighted <- function(lengths, probs, fraction = 0.5) {
    weights <- as.integer(probs / min(probs[probs > 0], na.rm = TRUE))
    i <- order(lengths)
    lengths <- rev(lengths[i])
    weights <- rev(weights[i])
    lengths <- rep(lengths, weights)
    lengths[cumsum(as.numeric(lengths)) >= sum(lengths) * (1 - fraction)][1]
}

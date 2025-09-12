# utility functions for loading RE site data and related
# caller must call stopSpinner()
reData_create <- "asNeeded"

reData_loadREData <- function(sourceId){
    sessionCache$get(
        'reData_loadREData', 
        key = sourceId, 
        permanent = FALSE, # already an RDS file 
        from = "ram",
        create = reData_create,
        createFn = function(...) {
            startSpinner(session, message = "loading RE data")
            readRDS(getSourceFilePath(sourceId, "reDataFile"))
        }  
    )$value
}
reData_getEnzymeName <- function(sourceId){
    sessionCache$get(
        'reData_getEnzymeName', 
        key = sourceId, 
        permanent = TRUE,
        from = "ram",
        create = reData_create,
        createFn = function(...) {
            startSpinner(session, message = "getting enzyme name")
            reData_loadREData(sourceId)$env$ENZYME_NAME
        }  
    )$value
}
reData_loadFilteringSites <- function(sourceId){
    sessionCache$get(
        'reData_loadFilteringSites', 
        key = sourceId, 
        permanent = TRUE,
        from = "ram",
        create = reData_create,
        createFn = function(...) {
            startSpinner(session, message = "loading filtering sites")
            reData_loadREData(sourceId)$filteringSites
        }  
    )$value
}
reData_loadEndpoints <- function(sourceIds){
    sessionCache$get(
        'reData_loadEndpoints', 
        key = sourceIds, 
        permanent = TRUE,
        from = "ram",
        create = reData_create,
        createFn = function(...) {
            startSpinner(session, message = "loading endpoints")
            x <- do.call(rbind, lapply(sourceIds, function(sourceId){
                reData_loadREData(sourceId)$endpoints
            }))
            startSpinner(session, message = "indexing endpoints")
            setkey(x, chrom, adjPos1)
            x
        }  
    )$value
}

# Classes ‘data.table’ and 'data.frame':  26361627 obs. of  7 variables:
#  $ chrom    : chr  "chr1" "chr1" "chr1" "chr1" ...
#  $ adjPos1  : num  14626 17401 49223 54341 54512 ...
#  $ refPos1  : int  14623 17398 49220 54338 54509 54669 54669 54904 54905 55577 .
# ..
#  $ refSpan  : int  150 150 150 172 172 236 142 236 150 150 ...
#  $ class    : Factor w/ 4 levels "R5","R3","L5",..: NA NA NA NA NA NA NA NA NA N
# A ...
#  $ inSilico : logi  FALSE FALSE FALSE FALSE TRUE FALSE ...
#  $ nObserved: int  1 3 1 4 4 8 2 8 1 1 ...
#  - attr(*, ".internal.selfref")=<externalptr>

##' @method [ diffAnalysisClass
##' @export
`[.diffAnalysisClass` <- function(x, i, j, asis = FALSE, ...) {
    result <- x@result
    y <- result[i,j, ...]
    if (!asis)
        return(y)
    x@result <- y
    return(x)
}

##' @method [[ diffAnalysisClass
##' @export
`[[.diffAnalysisClass` <- function(x, i) {
    result <- x@result
    if (!i %in% names(result))
        stop("input term not found...")
    result[[i]]
}

##' @method $ diffAnalysisClass
##' @export
`$.diffAnalysisClass` <-  function(x, name) {
    x <- x@result
    x[, name]
}

##' @method dim diffAnalysisClass
##' @export
dim.diffAnalysisClass <- function(x) {
    dim(x@result)
}

##' @importFrom utils tail
##' @method tail diffAnalysisClass
##' @export
tail.diffAnalysisClass <- function(x, n=6L, ...) {
    tail(x@result, n, ...)
}

#' @importFrom utils head
#' @method head diffAnalysisClass
#' @export
head.diffAnalysisClass <- function(x, n=6L, ...){
    head(as.data.frame(x), n=n, ...)
}

#' @method head alphasample
#' @export
head.alphasample <- function(x, n=6L, ...){
    head(as.data.frame(x), n=n, ...)
}

#' MPSE accessors
#' @param x MPSE object
#' @param i,j,... Indices specifying elements to extract or replace.
#' Indices are 'numeric' or 'character' vectors or empty (missing) or
#' NULL.  Numeric values are coerced to integer as by 'as.integer' 
#' (and hence truncated towards zero).  Character vectors will be matched 
#' to the 'names' of the object (or for matrices/arrays, the 'dimnames')
#' @param drop logical If 'TRUE' the result is coerced to the lowest 
#' possible dimension (see the examples).  This only works for extracting 
#' elements, not for the replacement.
#' @name MPSE-accessors
NULL

#' @rdname MPSE-accessors
#' @export
setMethod("[", signature(x="MPSE"),
          function(x, i, j, ..., drop=TRUE){
    otutree <- x@otutree
    refseq <- x@refseq
    nx <- methods::callNextMethod()
    newotus <- rownames(nx)
    if (!is.null(otutree)){
        rmotus <- setdiff(otutree@phylo$tip.label, newotus)
        otutree <- treeio::drop.tip(otutree, tip=rmotus)
    }
    if (!is.null(refseq)){
        refseq <- refseq[newotus]
    }
    nx@otutree <- otutree
    nx@refseq <- refseq
    return(nx)
})


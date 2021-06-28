#' @method as.data.frame diffAnalysisClass
#' @export
as.data.frame.diffAnalysisClass <- function(x,...){
    #efres <- tidyEffectSize(x)
    #kwres <- x@kwres
    #difftb <- merge(efres, kwres, by.x="f", by.y="f")
    #difftb <- difftb[order(difftb$pvalue),]
    #return(difftb)
    x@result
}

#' @method as.data.frame alphasample
#' @export
as.data.frame.alphasample <- function(x, ...){
    dat <- x@alpha
    if (!is.null(x@sampleda)){
        dat <- merge(dat, x@sampleda, by=0)
        rownames(dat) <- as.vector(dat$Row.names)
        dat$Row.names <- NULL
    }
    return(dat)
}

#' @method as.data.frame MPSE
#' @export
as.data.frame.MPSE <- function(x, ...){
    x <- as_tibble(x)
    return (x)
}


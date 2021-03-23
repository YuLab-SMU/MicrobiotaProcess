##' generate the result of rare curve.
##'
##' This function is designed to calculate the rare curve result of otu table
##' the result can be visualized by `ggrarecurve`.
##'
##' @title obtain the result of rare curve
##' @param obj phyloseq class or data.frame
##' shape of data.frame (nrow sample * ncol feature)
##' @param sampleda data.frame, (nrow sample * ncol factor)
##' @param chunks integer, the number of subsample in a sample, 
##' default is 400.
##' @param factorLevels list, the levels of the factors, default is NULL,
##' if you want to order the levels of factor, you can set this.
##' @param ..., additional parameters.
##' @return rarecurve class, which can be visualized by ggrarecurve
##' @author Shuangbin Xu
##' @export
##' @examples
##' data(test_otu_data)
##' set.seed(1024)
##' res <- get_rarecurve(test_otu_data, chunks=200)
##' p <- ggrarecurve(obj=res, 
##'                  indexNames=c("Observe","Chao1","ACE"),
##'                  shadow=FALSE,
##'                  factorNames="Group")
setGeneric("get_rarecurve", function(obj, ...)standardGeneric("get_rarecurve"))

#' @aliases get_rarecurve,data.frame
#' @rdname get_rarecurve
#' @export
setMethod("get_rarecurve", "data.frame", function(obj, sampleda, factorLevels=NULL, chunks=400){
    res <- stat_rare(data=obj, sampleda=sampleda, 
                     chunks=chunks, factorLevels=factorLevels, 
                     plotda=TRUE)
    res <- structure(list(data=res), class="rarecurve")
    return(res)

})

#' @aliases get_rarecurve,phyloseq
#' @rdname get_rarecurve
#' @export
setMethod("get_rarecurve", "phyloseq", function(obj, ...){
    otuda <- checkotu(obj)
    sampleda <- checksample(obj)
    res <- get_rarecurve(obj=otuda,
                         sampleda=sampleda,
                         ...)
    return(res)
})

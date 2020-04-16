#' @aliases get_vennlist,data.framet
#' @rdname get_vennlist
#' @export
setMethod("get_vennlist", "data.frame", function(obj,
    sampleinfo=NULL,
    factorNames=NULL,...){
    if (!is.null(sampleinfo) && !is.null(factorNames)){
    	sampleinfo <- sampleinfo[,match(factorNames, colnames(sampleinfo)), 
    							 drop=FALSE]
    }
    if (!is.null(sampleinfo) && is.null(factorNames)){
    	stop("when sampleinfo isn't NULL, factorNames shouldn't be NULL")
    }
    obj <- get_count(data=obj, 
                     featurelist=sampleinfo)
    vennlist <- apply(obj, 1, function(x){names(x[x>0])})
    return(vennlist)
})

#' @method get_vennlist default
#' @rdname get_vennlist
#' @export
get_vennlist.default <- function(obj,
    sampleinfo=NULL,
    factorNames=NULL,
    ...){
    if (!is.null(sampleinfo) && !is.null(factorNames)){
    	sampleinfo <- sampleinfo[,match(factorNames, colnames(sampleinfo)), 
    							 drop=FALSE]
    }
    if (!is.null(sampleinfo) && is.null(factorNames)){
    	stop("when sampleinfo isn't NULL, factorNames shouldn't be NULL")
    }
    obj <- CountOrRatios(obj, 
    			featurelist=sampleinfo, 
    			countmode=FALSE,
    			multiplenum=1,
    			rownamekeep=FALSE,
    			...)
    vennlist <- apply(obj, 1, function(x){names(x[x>0])})
    return(vennlist)
}

#' @method getvennlist default
#' @rdname getvennlist
#' @export
getvennlist.default <- function(data,
			   sampleinfo=NULL,
			   factorNames=NULL,
			   ...
			   ){
	if (!is.null(sampleinfo) && !is.null(factorNames)){
		sampleinfo <- sampleinfo[,match(factorNames, colnames(sampleinfo)), 
								 drop=FALSE]
	}
	if (!is.null(sampleinfo) && is.null(factorNames)){
		stop("when sampleinfo isn't NULL, factorNames shouldn't be NULL")
	}
	data <- CountOrRatios(data, 
				featurelist=sampleinfo, 
				countmode=FALSE,
				multiplenum=1,
   				rownamekeep=FALSE,
				...)
	vennlist <- apply(data, 1, function(x){names(x[x>0])})
	return(vennlist)
	
}

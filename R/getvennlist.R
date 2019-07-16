#' @method getvennlist default
#' @rdname getvennlist
#' @author Shuangbin Xu
#' @export
getvennlist.default <- function(da,
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
	data <- CountOrRatios(da, 
				featurelist=sampleinfo, 
				countmode=FALSE,
				multiplenum=1,
   				rownamekeep=FALSE,
				...)
	vennlist <- apply(data, 1, function(data){names(data[data>0])})
	return(vennlist)
	
}

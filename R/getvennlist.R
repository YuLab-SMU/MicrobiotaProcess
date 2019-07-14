#' @title generate a vennlist for VennDiagram
#'
#' @description
#' generate a vennlist as the input of \code{[VennDiagram]}.
#'
#' @param da dataframe; a dataframe contained one character column and the others are numeric.
#' all columns should be numeric if sampleinfo isn't NULL.
#' @param sampleinfo dataframe; a sample information, default is NULL.
#' @param factorNames character, a column name of sampleinfo,
#' when sampleinfo isn't NULL, factorNames shouldn't be NULL, default is NULL.
#' @param ...; Additional arguments passed to \code{\link[MicrobitaProcess]{CountOrRatios}}. 
#'
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

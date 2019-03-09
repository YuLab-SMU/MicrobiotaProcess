#' @title generate a vennlist for VennDiagram
#'
#' @description
#' generate a vennlist as the input of \code{[VennDiagram]}.
#'
#' @param da dataframe; a dataframe contained one character column and others is numeric 
#' if sampleinfo is NULL, else a numeirc dataframe.
#' @param sampleinfo dataframe; a sample information.
#' @param ...; Additional arguments passed to \code{\link[MicrobitaProcess]{CountOrRatios}}. 
#'
#' @author Shuangbin Xu
#' @export
#' @importFrom MicrobitaProcess CountOrRatios 

getvennlist <- function(da,
			   sampleinfo=NULL,
			   ...
			   ){
    	data <- da
	data <- CountOrRatios(data, 
				featurelist=sampleinfo, 
				countmode=FALSE,
				percentmode=FALSE,
			      	multiplenum=1,
   				rownamekeep=FALSE)
	vennlist <- apply(data, 1, function(data){names(data[data>0])})
	return(vennlist)
	
}

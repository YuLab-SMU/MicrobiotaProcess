#' @title split a dataframe contained one column
#'
#' @description
#' split a dataframe contained one column with a specify field separator character.
#' @param strdataframe dataframe; a dataframe contained one column to split.
#' @param prefix character; the result dataframe columns names prefix, default is "tax".
#' @param sep character; the field separator character, default is "; ".
#' @param extra character; See \code{[tidyr]{separate}} details.
#' @param fill character; See \code{[tidyr]{separate}} details.
#' @param ...; Additional arguments passed to \code{[tidyr]{separate}}.
#' @export
#' @author Shuangbin Xu
#' @importFrom tidyr separate
splitStrtoList <- function(strdataframe, 
			      prefix="tax", 
			      sep="; ", 
			      extra="drop",
			      fill="right",
			      ...){
    	extra <- match.arg(extra, c("drop","warn","merge"))
	fill <- match.arg(fill, c("warn", "right", "left"))
       colstr <- names(strdataframe)
	tmplength <- length(strsplit(as.character(strdataframe[1,1]), sep)[[1]])
	newcolnames <- paste(prefix, rep(1:tmplength), sep="")
	tmpdata <- tidyr::separate(strdataframe, 
				      colstr,
				      newcolnames,
				      sep=sep,
				      extra = "drop", 
				      fill = "right",...
				      )
	return(tmpdata)
}

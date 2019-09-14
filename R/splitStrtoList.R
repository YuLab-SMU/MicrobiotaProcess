#' @title split a dataframe contained one column
#'
#' @description
#' split a dataframe contained one column with a specify field separator character.
#' @param strdataframe dataframe; a dataframe contained one column to split.
#' @param prefix character; the result dataframe columns names prefix, default is "tax".
#' @param sep character; the field separator character, default is "; ".
#' @param extra character; See \code{\link[tidyr]{separate}} details.
#' @param fill character; See \code{\link[tidyr]{separate}} details.
#' @param ..., Additional arguments passed to \code{\link[tidyr]{separate}}.
#' @return data.frame of strdataframe by sep.
#' @export
#' @author Shuangbin Xu
#' @importFrom tidyr separate
#' @examples
#' otudafile <- system.file("extdata", "otu_tax_table.txt",
#'                       package="MicrobiotaProcess")
#' samplefile <- system.file("extdata",
#'                  "sample_info.txt", package="MicrobiotaProcess")
#' otuda <- read.table(otudafile, sep="\t", header=TRUE,
#'                     row.names=1, check.names=FALSE,
#'                     skip=1, comment.char="")
#' sampleda <- read.table(samplefile,
#'             sep="\t", header=TRUE, row.names=1)
#' taxdf <- otuda[!sapply(otuda, is.numeric)]
#' taxdf <- splitStrtoList(taxdf)
#' head(taxdf)
splitStrtoList <- function(strdataframe, 
    prefix="tax", 
    sep="; ", 
    extra="drop",
    fill="right",
    ...){
    extra <- match.arg(extra, c("drop","warn","merge"))
    fill <- match.arg(fill, c("warn", "right", "left"))
    colstr <- colnames(strdataframe)
    tmplength <- max(unlist(lapply(strsplit(as.vector(strdataframe[,1]), sep),
		         function(x)length(x))))
    newcolnames <- paste(prefix, rep(seq_len(tmplength)), sep="")
    tmpdata <- suppressWarnings(separate(strdataframe, 
                                         colstr,
                                         newcolnames,
                                         sep=sep,
                                         extra = "warn", 
                                         fill = "warn",
                                         ...))
    return(tmpdata)
}

#' @title Dropping Species with Few abundance and Few Occurrences
#' 
#' @description
#' Drop species or features from the feature data frame that occur fewer than or equal to a threshold number of 
#' occurrences and fewer abundance than to a threshold abundance.
#'
#' @param taxtab dataframe; a dataframe of species (or features), default is (n_sample, n_feature).
#' @param rmode boolean; whether transpose the taxtab, default is False.
#' @param minocc numeric; the threshold number of occurrences to be dropped, if < 1.0,
#' it will be the threshold ratios of occurrences, default is 0.
#' @param minabu numeric: the threshold abundance, if fewer than the threshold will be dropped, default is 0.
#' @return a list contained feature dataframe dropped, and the call, arguments.
#' @export
#' @author Shuangbin Xu
#' @examples
#' otudafile <- system.file("extdata", "otu_tax_table.txt",
#'                          package="MicrobiotaProcess")
#' otuda <- read.table(otudafile, sep="\t", 
#'                     header=TRUE, row.names=1, 
#'                     check.names=FALSE, skip=1, 
#'                     comment.char="")
#' otuda <- otuda[sapply(otuda, is.numeric)]
#' dim(otuda)
#' otudat <- droptax(otuda, minocc=0.2, minabu=10)
#' dim(otudat)
droptax <- function(taxtab, rmode=FALSE, minocc=0, minabu=0){
	if (isTRUE(rmode)){
		taxtab <- data.frame(t(taxtab), check.names=FALSE)
	}
	if (minocc < 1.0){
		minocc <- round(dim(taxtab)[1]*minocc, 0)
	}
	taxtab <- taxtab[,apply(taxtab>minabu,2,sum)>=minocc]
	return (taxtab)
}


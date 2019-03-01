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
#' 
#' library("microproccess")
#' data <- read.csv(system.file("data", package="microproccess", "Baxter_16_crc_genera_group.csv.gz"))
#' data$Group <- NULL
#' dim(data)
#' head(data)
#' newdat <- droptax(data, rmode=FALSE, minocc=0.2, minabu=0)
#' dim(newdat)
#' head(newdat)

droptax <- function(taxtab, rmode=FALSE, minocc=0, minabu=0){
	if (isTRUE(rmode)){
		taxtab <- data.frame(t(taxtab), check.names=FALSE)
	}
	if (minocc < 1.0){
		minocc <- round(dim(taxtab)[1]*minocc, 0)
	}
	taxtab <- taxtab[,apply(taxtab>minabu,2,sum)>=minocc]
	attr(taxtab, "call") <- match.call()
	attr(taxtab, "rmode") <- rmode
	attr(taxtab, "minocc") <- minocc
	attr(taxtab, "minabu") <- minabu
	return (taxtab)
}


#' @title Two sample Test
#' 
#' @description
#' Performs two sample test with `kruskal.test`, `wilcox.test`, `t.test`.
#' @param measure character; must be the names of feature.
#' @param methodnames character; the method of tow sample test, default is
#' `kruskal.test`, one of the `t.test`, `kruskal.test` and `wilcox.test`.
#' @param data dataframe; the data contaiend the features and the group for
#' the corresponding elements of `features`, shape(`row_sample * col_feature`).
#' @param ...; Additional arguments passed to the method choosed.
#' @author Shuangbin Xu
#' @export
TwoGroupTest <- function(measure, methodnames="kruskal.test", data, ...){
	methodnames <- match.arg(methodnames, c("t.test", "kruskal.test", "wilcox.test"))
	colnames(data) <- make.names(colnames(data))
	measure <- make.names(measure)
	numsindex <- !unlist(lapply(data, is.numeric))
	groupname <- names(data[,numsindex,drop=FALSE])
	formulatmp <- as.formula(paste(measure, groupname, sep=" ~ "))
	testrt <- do.call(methodnames, list(formulatmp, data=data, ...))
	return(testrt)
}

#' @title multiple comparison for multiple measure with tow sample test
#'
#' @description
#' Performs two sample test to multiple measure with 
#' `kruskal.test`, `wilcox.test`, `t.test`.
#' @param methodnames character; the method of tow sample test, default is 
#' `kruskal.test`, one of the `t.test`, `kruskal.test` and `wilcox.test`.
#' @param measurelist vector; must be the names of the features, default is NULL,
#' meaning the total features.
#' @param data dataframe; the data contaiend the features and the group for
#' the corresponding elements of `features`, shape(`row_sample * col_feature`).
#' @param ...; Additional arguments passed to the method choosed.
#' @author Shuangbin Xu
#' @export 
mapplyTwoGroupTest <- function(methodnames="kruskal.test", measurelist=NULL, data, ...){
	methodnames <- match.arg(methodnames, c("t.test", "kruskal.test", "wilcox.test"))
	colnames(data) <- make.names(colnames(data))
	numsindex <- !unlist(lapply(data, is.numeric))
	groupname <- names(data[, numsindex, drop=FALSE])
	if (is.null(measurelist)){
		nums <- unlist(lapply(data, is.numeric)) 
		measurelist <- names(data[, nums, drop=FALSE])
	}
	testresult <- mapply(TwoGroupTest,
						 measurelist,
						 MoreArgs=list(data=data,
									   methodnames=methodnames,
									   ...),
						 SIMPLIFY=FALSE)
	return(testresult)
}


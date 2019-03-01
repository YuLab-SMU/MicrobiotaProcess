#' @title caculate the count or relative abundance of replicate element with a speficify column
#' 
#' @description
#' Caculate the count or relative abundance of replicate element with a speficify columns
#' 
#' @param da dataframe; a dataframe contained one character column and others is numeric 
#' if featurelist is NULL, else a numeirc dataframe.
#' @param featurelist dataframe; a dataframe contained one chatacter column, default is NULL.
#' @param countmode boolean; whether return the count results (FALSE),or relative abundance (TRUE).
#' @param percentmode boolean; whether return the 100% percent (FALSE),or (0.1) percent(TRUE).
#' @param multiplenum numeric; the multiple you want to increase,default is 1.
#' @param rownamekeep boolean; whether you return a dataframe contained the rownames,default is FALSE.
#' @export 
#' @author Shuangbin Xu
#' @importFrom plyr ddply numcolwise
CountOrRatios <- function(da, 
			     featurelist=NULL, 
			     countmode=FALSE, 
			     percentmode=FALSE,
			     multiplenum=1, 
			     rownamekeep=FALSE){
	data <- da
       if (!is.null(featurelist)){
		data <- merge(da, featurelist, by=0)
		rownames(data) <- data$Row.names
              data$Row.names <- NULL
	}
	nums <- !unlist(lapply(data, is.numeric))
	group <- names(data[,nums,drop=FALSE])
       data <- data.frame(plyr::ddply(data, group, plyr::numcolwise(sum)), 
			     check.names=F, stringsAsFactors=FALSE)   
       rownames(data) <- as.vector(data[[group]])
       data[[group]] <- NULL
	if (!isTRUE(countmode)){
       	data <- data.frame(prop.table(as.matrix(data), 2), check.names=F, stringsAsFactors=FALSE)
		data <- data*multiplenum
	}
      	if (!isTRUE(countmode) && isTRUE(percentmode)){
       	data <- data*100
	}
	if (isTRUE(rownamekeep)){
		data <- data.frame(cbind(feature=rownames(data), data), check.names=F, stringsAsFactors=FALSE)
	}
       return (data)
}


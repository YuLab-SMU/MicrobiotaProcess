#' @keywords internal
multi.compare <- function(fun = wilcox.test, 
						data, 
						feature, 
						factorNames, ...){ 
	sapply(simplify = FALSE,
		   feature,
		   function(x){
		   tmpformula <- as.formula(paste(x, factorNames, sep="~"))
		   do.call(fun,list(tmpformula,data=data, ...))})
}	 

#' @keywords internal
classlevels <- function(sampleda, class){
	levelstmp <- levels(sampleda[,match(class, colnames(sampleda))])
	return(levelstmp)
}

#' @keywords internal
getclass2sub <- function(sampleda, class, subclass){
	tmpsplit <- sampleda[,match(class, colnames(sampleda))]
	if (!missing(subclass)){
		samplelist <- split(sampleda, tmpsplit) 
		lapply(samplelist, function(x)unique(as.vector(x[,match(subclass, colnames(x))])))
	}
}


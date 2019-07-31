#' @title Differential expression analysis based on kruskal.test, generalized fold change and wilcox.test
#' @author Shuangbin Xu
#' @export
diffAnalysis <- function(x, ...){
	UseMethod("diffAnalysis")
}

#' @method diffAnalysis data.frame
#' @rdname diffAnalysis
#' @export
diffAnalysis.data.frame <- function(data, 
									sampleda, 
									class, 
									subclass=NULL, 
									padjust="fdr",
									...){
	vars <- colnames(data)
	datameta <- merge(data, sampleda, by=0)
	kwres <- multi.compare(fun=kruskal.test,
						   data=datameta,
						   feature=vars,
						   factorNames=class)
	kwres <- lapply(kwre, function(x)x$p.value)
	kwres <- do.call("rbind", kwres)
	kwres$fdr <- p.adjust(kwres, method=padjust)
    classlevels <- getclasslevels(sampleda, class)
	compareclass <- getcompareclass(classlevels)
	if (!is.null(subclass)){
		class2sub <- getclass2sub(sampleda, class, subclass)
		comsubclass <- apply(compareclass,1,
							 function(x)getcomparesubclass(x[1],x[2],class2sub))
	}

}


###' @method diffAnalysis phyloseq
##' @rdname diffAnalysis
##' @export
#diffAnalysis.phyloseq <- function(obj,...){
#
#
#
#}

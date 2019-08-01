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
									firstcomfun='kruskal.test',
									padjust="fdr",
									filtermod="pvalue",
									firstalpha=0.05,
									strictmod=TRUE,
									fcfun="generalizedFC",
									secondcomfun="wilcox.test",
									clmin=5,
									clwilc=TRUE,
									submin=5,
									subclwilc=TRUE,
									secondalpha=0.05,
									...){
	vars <- colnames(data)
	datameta <- merge(data, sampleda, by=0)
	kwres <- multi.compare(fun=firstcomfun,
						   data=datameta,
						   feature=vars,
						   factorNames=class)
	kwres <- lapply(kwre, function(x)x$p.value)
	kwres <- do.call("rbind", kwres)
	kwres$fdr <- p.adjust(kwres, method=padjust)
	if (!filtermod=="pvalue"){
		varsfirst <- rownames(kwres[kwres$fdr <= firstalpha,])
	}else{
		varsfirst <- rownames(kwres[kwres[,1]<=0.05,,drop=FALSE])
	}
    classlevels <- getclasslevels(sampleda, class)
	compareclass <- getcompareclass(classlevels)
	if (!is.null(subclass) && strictmod){
		class2sub <- getclass2sub(sampleda, class, subclass)
		comsubclass <- apply(compareclass,1,
							 function(x)getcomparesubclass(x[1],x[2],class2sub))
		secondvars <- diffsubclass(datasample=datameta,
					 	features=varsfirst,
					 	comsubclass=comsubclass,
					 	class=class,
					 	subclass=subclass,
					 	fcfun=fcfun,
					 	secondcomfun=secondcomfun,
					 	submin=submin,
					 	subclwilc=subclwilc,
					 	pfold=secondalpha)

	}else{
		secondvars <- diffclass(datasample=datameta,
								features=varsfirst,
								comclass=compareclass,
								class=class,
								fcfun=fcfun,
								secondcomfun=secondcomfun,
								classmin=clmin,
								clwilc=clwilc,
								pfold=secondalpha,
								)
	}
	secondvars <- getconsistentfeatures(diffsubclassfeature=secondvars,
										class=class,
										classlevels=classlevels)
	
}


###' @method diffAnalysis phyloseq
##' @rdname diffAnalysis
##' @export
#diffAnalysis.phyloseq <- function(obj,...){
#
#
#
#}

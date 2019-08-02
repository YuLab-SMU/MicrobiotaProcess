#' @title Differential expression analysis based on kruskal.test, generalized fold change and wilcox.test
#' @param data data.frame, (nrow sample * ncol features)
#' @param obj object,a phyloseq class contained otu_table, sample_data, taxda.
#' @param sampleda data.frame, nrow sample * ncol factor, the sample names of sampleda and data should be the same.
#' @param class character, the factor name in sampleda.
#' @param subclass character, the factor name in sampleda, default is NULL, 
#' meaning no subclass compare.
#' @param taxda data.frame, the classification of the feature in data. default is NULL.
#' @param alltax logical, whether to set all classification as features if taxda is not NULL, 
#' default is TRUE.
#' @param mlfun character, the method for calculating the effect size of features, 
#' choose "lda" or "rf", default is "lda".
#' @param ratio numeric, range from 0 to 1, the proportion of samples for calculating the effect 
#' size of features, default is 0.7. 
#' @param firstcomfun character, the method for first test, "oneway.test" for normal distributions, 
#' suggested choosing "kruskal.test" for uneven distributions, default is "kruskal.test".
#' @param padjust character, the correction method, default is "fdr".
#' @param filtermod character, the method to filter, default is "pvalue".
#' @param firstalpha numeric, the alpha value for the first test, default is 0.05.
#' @param strictmod logical, whether to performed in one-against-one, default is TRUE (strict).
#' @param fcfun character, default is "generalizedFC", it can't be set another at the present time.
#' @param secondcomfun character, the method for one-against-one, default is "wilcox.test" for uneven distributions.
#' @param clmin integer, the minimum number of samples per class or suclass for performing test, default is 5.
#' @param clwilc logical, whether to perform test of per class or subclass, default is TRUE.
#' @param secondalpha numeric, the alpha value for the second test, default is 0.05.
#' @param normalization integer, set the normalization value, set a big number if to get more meaningful values 
#' for the LDA score, or you can set NULL for no normalization, default is 1000000.
#' @param bootnums integer, set the number of bootstrap iteration for lda or rf, default is 30.
#' @author Shuangbin Xu
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr select
#' @importFrom magrittr %>%
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
									taxda=NULL,
									alltax=TRUE,
									mlfun="lda",
									ratio=0.7,
									firstcomfun='kruskal.test',
									padjust="fdr",
									filtermod="pvalue",
									firstalpha=0.05,
									strictmod=TRUE,
									fcfun="generalizedFC",
									secondcomfun="wilcox.test",
									clmin=5,
									clwilc=TRUE,
									secondalpha=0.05,
									normalization=1000000,
									bootnums=30,
									...){
	if (!is.null(taxda) && alltax){
		taxda <- fillNAtax(taxda)
		data <- getalltaxdf(data, taxda)
	}
	vars <- colnames(data)
	datameta <- merge(data, sampleda, by=0)
	kwres <- multi.compare(fun=firstcomfun,
						   data=datameta,
						   feature=vars,
						   factorNames=class)
	kwres <- lapply(kwres, function(x)x$p.value)
	kwres <- do.call("rbind", kwres)
	kwres <- data.frame(f=rownames(kwres),pvalue=kwres[,1])
	kwres$fdr <- p.adjust(kwres$pvalue, method=padjust)
	if (!filtermod=="pvalue"){
		varsfirst <- as.vector(kwres[kwres$fdr <= firstalpha,,drop=FALSE]$f)
	}else{
		varsfirst <- as.vector(kwres[kwres$pvalue<=firstalpha,,drop=FALSE]$f)
	}
	if (!length(varsfirst)>0){stop("There are not significantly discriminative features before internal wilcoxon!")}
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
					 	submin=clmin,
					 	subclwilc=clwilc,
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
	if (!length(secondvars)>0){stop("There are not significantly discriminative features after internal wilcoxon!")}
	leaveclasslevels <- unlist(lapply(names(secondvars), 
									 function(x){unlist(strsplit(x,"-vs-"))}))
	secondvars <- getconsistentfeatures(diffsubclassfeature=secondvars,
										class=class,
										classlevels=leaveclasslevels)
	secondvarsvectors <- getsecondvarlist(secondvars=secondvars)
	if (!is.null(normalization)){
		data <- data * normalization
	}
	dameta <- merge(data, sampleda, by=0) %>% column_to_rownames(var="Row.names")
	dameta <- dameta %>% select(c(secondvarsvectors, class))
	dameta <- split(dameta, dameta[[class]])
	dameta <- sampledflist(dameta, bootnums=bootnums, ratio=ratio)
	if (mlfun=="lda"){
		mlres <- LDAeffectsize(dameta, compareclass, class, bootnums=bootnums) 
	}
	if (mlfun=="rf"){
		mlres <- rfimportance(dameta, class, bootnums=bootnums)
	}
	res <- new("diffAnalysisClass",
		originalD=data,
		sampleda=sampleda,
		taxda=taxda,
		kwres=kwres,
		secondvars=secondvars,
		mlres=mlres)
	return(res)
}


#' @method diffAnalysis phyloseq
#' @rdname diffAnalysis
#' @importFrom phyloseq tax_table
#' @export
diffAnalysis.phyloseq <- function(obj, class, subclass=NULL,...){
	otuda <- checkotu(obj)
	sampleda <- checksample(obj)
	taxda <- tax_table(obj)
	res <- diffAnalysis.data.frame(data=otuda, 
							sampleda=sampleda, 
							class=class,
							subclass=subclass,
							...)
	return(res)
}


#' @importFrom dplyr bind_rows
#' @importFrom magrittr %>%
#' @importFrom tibble column_to_rownames
#' @keywords internal
getalltaxdf <- function(data, taxda){
	data <- data.frame(t(data), check.names=FALSE)
	dat <- suppressWarnings(mapply(CountOrRatios, 
							   taxda, 
							   MoreArgs=list(da=otuda,
											 countmode=FALSE,
											 rownamekeep=TRUE),
							SIMPLIFY=FALSE) %>% bind_rows()) %>%
	       column_to_rownames(var="feature")
   return(dat)

}


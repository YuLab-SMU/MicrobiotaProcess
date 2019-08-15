#' @title Differential expression analysis based on kruskal.test, generalized fold change and wilcox.test
#' @param obj object,a phyloseq class contained otu_table, sample_data, taxda, 
#' or data.frame, nrow sample * ncol features.
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
#' @param clmin integer, the minimum number of samples per class for performing test, default is 5.
#' @param clwilc logical, whether to perform test of per class, default is TRUE.
#' @param secondalpha numeric, the alpha value for the second test, default is 0.05.
#' @param subclmin integer, the minimum number of samples per suclass for performing test, default is 3.
#' @param subclwilc logical, whether to perform test of per subclass, default is TRUE, meaning more strict. 
#' @param ldascore numeric, the threshold on the absolute value of the logarithmic LDA score, default is 2.
#' @param normalization integer, set the normalization value, set a big number if to get more meaningful values 
#' for the LDA score, or you can set NULL for no normalization, default is 1000000.
#' @param bootnums integer, set the number of bootstrap iteration for lda or rf, default is 30.
#' @param ..., additional parameters.
#' @return diffAnalysis class.
#' @author Shuangbin Xu
#' @importFrom tibble column_to_rownames
#' @importFrom dplyr select
#' @importFrom magrittr %>%
#' @export
#' @examples
#' data(kostic2012crc)
#' kostic2012crc
#' head(phyloseq::sample_data(kostic2012crc),3)
#' kostic2012crc <- phyloseq::rarefy_even_depth(kostic2012crc,rngseed=1024)
#' table(phyloseq::sample_data(kostic2012crc)$DIAGNOSIS)
#' set.seed(1024)
#' diffres <- diffAnalysis(kostic2012crc, class="DIAGNOSIS",
#'                         mlfun="lda", filtermod="fdr",
#'                         firstcomfun = "kruskal.test",
#'                         firstalpha=0.05, strictmod=TRUE,
#'                         secondcomfun = "wilcox.test",
#'                         submin=3, subclwilc=TRUE,
#'                         secondalpha=0.01, lda=3)
diffAnalysis <- function(obj, ...){
	UseMethod("diffAnalysis")
}

#' @method diffAnalysis data.frame
#' @rdname diffAnalysis
#' @importFrom tibble column_to_rownames
#' @importFrom stats p.adjust
#' @export
diffAnalysis.data.frame <- function(obj, sampleda, class, subclass=NULL, taxda=NULL,alltax=TRUE, mlfun="lda", 
    ratio=0.7,	firstcomfun='kruskal.test',	padjust="fdr",filtermod="pvalue",
    firstalpha=0.05, strictmod=TRUE, fcfun="generalizedFC",	secondcomfun="wilcox.test",
    clmin=5, clwilc=TRUE, secondalpha=0.05,	subclmin=3,	subclwilc=TRUE,	ldascore=2,
    normalization=1000000, bootnums=30,	...){
    if (!is.null(taxda)){
    	taxda <- fillNAtax(taxda)
    	if (alltax){obj <- getalltaxdf(obj, taxda)}
    }
    sampleda <- sampleda %>% select(c(class, subclass))
    if (ncol(sampleda)>1){sampleda <- duplicatedtaxcheck(sampleda) %>% column_to_rownames(var="rowname")}
    vars <- colnames(obj)
    datameta <- merge(obj, sampleda, by=0)
    kwres <- multi.compare(fun=firstcomfun, data=datameta, feature=vars, factorNames=class)
    kwres <- lapply(kwres, function(x)x$p.value)
    kwres <- do.call("rbind", kwres)
    rownames(kwres) <- vars
    kwres <- data.frame(f=rownames(kwres),pvalue=kwres[,1])
    kwres$fdr <- p.adjust(kwres$pvalue, method=padjust)
    if (!filtermod=="pvalue"){varsfirst <- as.vector(kwres[kwres$fdr<=firstalpha& !is.na(kwres$fdr),,drop=FALSE]$f)
	}else{varsfirst <- as.vector(kwres[kwres$pvalue<=firstalpha&!is.na(kwres$pvalue),,drop=FALSE]$f)}
    if (!length(varsfirst)>0){stop("There are not significantly discriminative features before internal wilcoxon!")}
    classlevels <- getclasslevels(sampleda, class)
    compareclass <- getcompareclass(classlevels)
    if (!is.null(subclass) && strictmod){
    	class2sub <- getclass2sub(sampleda, class, subclass)
    	comsubclass <- apply(compareclass,1,function(x)getcomparesubclass(x[1],x[2],class2sub))
    	secondvars <- diffsubclass(datasample=datameta, features=varsfirst, comsubclass=comsubclass,class=class,subclass=subclass,
    				 	fcfun=fcfun, secondcomfun=secondcomfun, submin=subclmin, subclwilc=subclwilc, pfold=secondalpha)
    }else{
    	secondvars <- diffclass(datasample=datameta, features=varsfirst, comclass=compareclass, class=class, fcfun=fcfun,
    							secondcomfun=secondcomfun,classmin=clmin,clwilc=clwilc,pfold=secondalpha)
	}
    if (!length(secondvars)>0){stop("There are not significantly discriminative features after internal wilcoxon!")}
    leaveclasslevels <- unlist(lapply(names(secondvars), function(x){unlist(strsplit(x,"-vs-"))}))
    secondvars <- getconsistentfeatures(diffsubclassfeature=secondvars,	class=class,classlevels=leaveclasslevels)
    secondvarsvectors <- getsecondvarlist(secondvars=secondvars)
    if (!is.null(normalization)){obj <- obj * normalization}
    dameta <- merge(obj, sampleda, by=0) %>% column_to_rownames(var="Row.names")
    dameta <- dameta %>% select(c(secondvarsvectors, class))
    dameta <- split(dameta, dameta[[class]])
    dameta <- sampledflist(dameta, bootnums=bootnums, ratio=ratio)#, randomSeed=1024)
    if (mlfun=="lda"){mlres <- LDAeffectsize(dameta, compareclass, class, bootnums=bootnums, LDA=ldascore)}
    if (mlfun=="rf"){mlres <- rfimportance(dameta, class, bootnums=bootnums)}
    res <- new("diffAnalysisClass",originalD=obj,sampleda=sampleda,taxda=taxda,kwres=kwres,
    		   secondvars=secondvars,mlres=mlres,classname=class,normalization=normalization)
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
    res <- diffAnalysis.data.frame(obj=otuda, 
    						sampleda=sampleda, 
    						class=class,
    						subclass=subclass,
    						taxda=taxda,
    						...)
    return(res)
}


#' @importFrom dplyr bind_rows
#' @importFrom magrittr %>%
#' @keywords internal
getalltaxdf <- function(data, taxda){
    data <- data.frame(t(data), check.names=FALSE)
    dt <- list()
    for (i in seq_len(ncol(taxda))){
    	dat <- CountOrRatios(data, taxda[,i,drop=FALSE], countmode=FALSE, rownamekeep=FALSE)
    	dt[[i]] <- dat
    }
    dt <- do.call("rbind", dt)
    dt <- data.frame(t(dt), check.names=FALSE)
    return(dt)
}

#' @title get the table of diffAnalysisClass
#' @param obj object, diffAnalysisClass
#' @return a data.frame contained results of diffAnalysis
#' @export
#' @examples
#' data(kostic2012crc)
#' kostic2012crc
#' head(phyloseq::sample_data(kostic2012crc),3)
#' kostic2012crc <- phyloseq::rarefy_even_depth(kostic2012crc,rngseed=1024)
#' table(phyloseq::sample_data(kostic2012crc)$DIAGNOSIS)
#' set.seed(1024)
#' diffres <- diffAnalysis(kostic2012crc, class="DIAGNOSIS",
#'                         mlfun="lda", filtermod="fdr",
#'                         firstcomfun = "kruskal.test",
#'                         firstalpha=0.05, strictmod=TRUE,
#'                         secondcomfun = "wilcox.test",
#'                         submin=3, subclwilc=TRUE,
#'                         secondalpha=0.01, lda=3)
#' restab <-tidydiffAnalysis(diffres)
#' head(restab)
tidydiffAnalysis <- function(obj){
    efres <- tidyEffectSize(obj)
    kwres <- obj@kwres
    difftb <- merge(efres, kwres, by.x="f", by.y="f")
    return(difftb)
}

#' @keywords internal
tidyEffectSize <- function(obj){
    f <- LDA <- MeanDecreaseAccuracy <- NULL
    secondvars <- getsecondTRUEvar(obj)
    efres <- merge(obj@mlres, secondvars, by.x="f", by.y="f") %>%
    		select (-c("gfc", "Freq"))
    if ("LDA" %in% colnames(efres)){
    	efres <- efres %>% mutate(f = factor(f, levels=f[order(eval(parse(text=obj@classname)), LDA)]))
    }else{
    	efres <- efres %>% mutate(f = factor(f, levels=f[order(eval(parse(text=obj@classname)),
    														   MeanDecreaseAccuracy)]))
    }
    return(efres)
}

#' @keywords internal
getsecondTRUEvar <- function(obj){
    secondvars <- do.call("rbind",c(obj@secondvars,
    								make.row.names=FALSE))
    secondvars <- secondvars %>% filter(eval(parse(text="gfc"))%in%"TRUE")
    return(secondvars)
}

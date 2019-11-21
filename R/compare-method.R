#' @title Differential expression analysis
#' @param obj object,a phyloseq class contained otu_table, sample_data, taxda, 
#' or data.frame, nrow sample * ncol features.
#' @param sampleda data.frame, nrow sample * ncol factor, the sample names of 
#' sampleda and data should be the same.
#' @param class character, the factor name in sampleda.
#' @param subclass character, the factor name in sampleda, default is NULL, 
#' meaning no subclass compare.
#' @param taxda data.frame, the classification of the feature in data. 
#' default is NULL.
#' @param alltax logical, whether to set all classification as features if taxda is not NULL, 
#' default is TRUE.
#' @param standard_method character, the method of standardization, 
#' see also \code{\link[vegan]{decostand}}, default is NULL, 
#' it represents that the relative abundance of taxonomy will be used.
#' @param mlfun character, the method for calculating the effect size of features, 
#' choose "lda" or "rf", default is "lda".
#' @param ratio numeric, range from 0 to 1, the proportion of samples for calculating the effect 
#' size of features, default is 0.7. 
#' @param firstcomfun character, the method for first test, "oneway.test" for normal distributions, 
#' suggested choosing "kruskal.test" for uneven distributions, default is "kruskal.test", or
#' you can use lm, glm, or glm.nb (for negative binomial distribution), or `kruskal_test`, 
#' `oneway_test` of `coin`.
#' @param padjust character, the correction method, default is "fdr".
#' @param filtermod character, the method to filter, default is "pvalue".
#' @param firstalpha numeric, the alpha value for the first test, default is 0.05.
#' @param strictmod logical, whether to performed in one-against-one, default is TRUE (strict).
#' @param fcfun character, default is "generalizedFC", it can't be set another at the present time.
#' @param secondcomfun character, the method for one-against-one, default is "wilcox.test" for 
#' uneven distributions, or `wilcox_test` of `coin`, or you can also use `lm`,
#' `glm`, `glm.nb`(for negative binomial distribution in `MASS`).
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
#' @return diff_analysis class.
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
#' diffres <- diff_analysis(kostic2012crc, class="DIAGNOSIS",
#'                         mlfun="lda", filtermod="fdr",
#'                         firstcomfun = "kruskal.test",
#'                         firstalpha=0.05, strictmod=TRUE,
#'                         secondcomfun = "wilcox.test",
#'                         subclmin=3, subclwilc=TRUE,
#'                         secondalpha=0.01, ldascore=3)
diff_analysis <- function(obj, ...){
	UseMethod("diff_analysis")
}

#' @method diff_analysis data.frame
#' @rdname diff_analysis
#' @importFrom tibble column_to_rownames
#' @importFrom stats p.adjust
#' @export
diff_analysis.data.frame <- function(obj, sampleda, class, subclass=NULL, taxda=NULL,alltax=TRUE, standard_method=NULL, mlfun="lda", 
    ratio=0.7, firstcomfun='kruskal.test', padjust="fdr",filtermod="pvalue",
    firstalpha=0.05, strictmod=TRUE, fcfun="generalizedFC", secondcomfun="wilcox.test",
    clmin=5, clwilc=TRUE, secondalpha=0.05, subclmin=3, subclwilc=TRUE,	ldascore=2,
    normalization=1000000, bootnums=30,	...){
    if (!is.null(taxda)){taxda <- fillNAtax(taxda)
        if (alltax){obj <- getalltaxdf(obj, taxda, method=standard_method)}
    }
    sampleda <- sampleda %>% select(c(class, subclass))
    if (ncol(sampleda)>1){sampleda <- duplicatedtaxcheck(sampleda) %>% column_to_rownames(var="rowname")}
    vars <- colnames(obj)
    datameta <- merge(obj, sampleda, by=0)
    kwres <- multi.compare(fun=firstcomfun, data=datameta, feature=vars, factorNames=class)
    kwres <- lapply(kwres, function(x)getpvalue(x))
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
                                   fcfun=fcfun, secondcomfun=secondcomfun, submin=subclmin, subclwilc=subclwilc, pfold=secondalpha, ...)
    }else{
        secondvars <- diffclass(datasample=datameta, features=varsfirst, comclass=compareclass, class=class, fcfun=fcfun,
                                secondcomfun=secondcomfun,classmin=clmin,clwilc=clwilc,pfold=secondalpha, ...)
    }
    if (!length(secondvars)>0){stop("There are not significantly discriminative features after internal wilcoxon!")}
    leaveclasslevels <- unlist(lapply(names(secondvars), function(x){unlist(strsplit(x,"-vs-"))}))
    secondvars <- getconsistentfeatures(diffsubclassfeature=secondvars,	class=class,classlevels=leaveclasslevels)
    secondvarsvectors <- getsecondvarlist(secondvars=secondvars)
    if (!is.null(normalization)){obj <- obj * normalization}
    dameta <- merge(obj, sampleda, by=0) %>% column_to_rownames(var="Row.names")
    dameta <- dameta %>% select(c(secondvarsvectors, class))
    dameta <- split(dameta, dameta[,match(class,colnames(dameta))])
    dameta <- sampledflist(dameta, bootnums=bootnums, ratio=ratio)
    dameta <- removeconstant(dameta) 
    if (mlfun=="lda"){mlres <- LDAeffectsize(dameta, compareclass, class, bootnums=bootnums, LDA=ldascore)}
    if (mlfun=="rf"){mlres <- rfimportance(dameta, class, bootnums=bootnums, effsize=ldascore)}
    tmpfun <- ifelse(!"funname" %in% names(match.call()),NA,"diff_analysis.data.frame")
    res <- new("diffAnalysisClass",originalD=obj,sampleda=sampleda,taxda=taxda,kwres=kwres,
               secondvars=secondvars,mlres=mlres,call=match.call.defaults(fun=tmpfun))
    return(res)
}


#' @method diff_analysis phyloseq
#' @rdname diff_analysis
#' @importFrom phyloseq tax_table
#' @export
diff_analysis.phyloseq <- function(obj, ...){
    otuda <- checkotu(obj)
    sampleda <- checksample(obj)
    taxda <- tax_table(obj)
    call <- match.call()
    res <- diff_analysis.data.frame(obj=otuda,
                                    sampleda=sampleda,
                                    taxda=taxda,
                                    funname=TRUE,
                                    ...)
    for (i in setdiff(names(as.list(res@call)), names(call))){
        call[i] <- list(as.list(res@call)[[i]])
    }
    res@call <- call
    return(res)
}

#' @title get the table of abundance of all level taxonomy 
#'
#' @description
#' This function was designed to get the abundance of all level taxonomy,
#' the input can be phyloseq object or data.frame.
#' @param obj object, phyloseq or data.frame
#' @param method character, the normalization method, 
#' see also \code{\link[vegan]{decostand}}, default is NULL, the relative abundance 
#' will return.
#' @param taxda data.frame, the taxonomy table.
#' @param taxa_are_rows logical, if the obj is data.frame, and the features are rownames,
#' the taxa_are_rows should be set TRUE, default FALSE, meaning the features are colnames. 
#' @param ..., additional parameters, see also \code{\link[vegan]{decostand}}.
#' @return the all taxonomy abundance table
#' @author Shuangbin Xu
#' @export
#' @examples
#' data(test_otu_data)
#' alltaxatab <- get_alltaxadf(test_otu_data)
#' head(alltaxatab[,1:10])
setGeneric("get_alltaxadf", function(obj, ...){standardGeneric("get_alltaxadf")})

#' @aliases get_alltaxadf,phyloseq
#' @rdname get_alltaxadf
#' @importFrom phyloseq tax_table
#' @export
setMethod("get_alltaxadf", "phyloseq",function(obj, ...){
    otuda <- checkotu(obj)
    if (is.null(obj@tax_table)){
        stop("The taxaonomy table is empty!")
    }else{
        taxa <- fillNAtax(tax_table(obj)) 
    }
    data <- getalltaxdf(data=otuda, taxda=taxa, ...)
    return(data)
})

#' @aliases get_alltaxadf,data.frame
#' @rdname get_alltaxadf
#' @export
setMethod("get_alltaxadf", "data.frame", function(obj, taxda, taxa_are_rows=FALSE, ...){
    if (!taxa_are_rows){
        obj <- data.frame(t(obj), check.names=FALSE)
    }
    data <- getalltaxdf(data=obj, taxda=taxda, ...)
    return(data)
})

#' @importFrom magrittr %>%
#' @keywords internal
getalltaxdf <- function(data, taxda, method=NULL, ...){
    data <- data.frame(t(data), check.names=FALSE)
    dt <- list()
    for (i in seq_len(ncol(taxda))){
        if (is.null(method)){
            dat <- CountOrRatios(data, taxda[,i,drop=FALSE],
                                 countmode=FALSE, rownamekeep=FALSE)
        }else{
            dat <- CountOrRatios(data, taxda[,i,drop=FALSE], 
                                 countmode=TRUE, rownamekeep=FALSE)
            if(method=="count"){
                dat <- dat
            }else{
                dat <- transformdf(data=dat, method=method, ...)
            }
        }
        dt[[i]] <- dat
    }
    dt <- do.call("rbind", dt)
    dt <- data.frame(t(dt), check.names=FALSE)
    return(dt)
}

#' @importFrom vegan decostand
#' @keywords internal
transformdf <- function(data, method, ...){
    data <- data.frame(t(data), check.names=FALSE)
    data <- decostand(data, method=method, ...)
    data <- data.frame(t(data), check.names=FALSE)
    return(data)
}

###' @title get the table of diffAnalysisClass
###' @param x object, diffAnalysisClass
###' @param ..., additional parameters
###' @return a data.frame contained results of diff_analysis
###' @export
###' @examples
###' data(kostic2012crc)
###' kostic2012crc
###' head(phyloseq::sample_data(kostic2012crc),3)
###' kostic2012crc <- phyloseq::rarefy_even_depth(kostic2012crc,rngseed=1024)
###' table(phyloseq::sample_data(kostic2012crc)$DIAGNOSIS)
###' set.seed(1024)
###' diffres <- diff_analysis(kostic2012crc, class="DIAGNOSIS",
###'                         mlfun="lda", filtermod="fdr",
###'                         firstcomfun = "kruskal.test",
###'                         firstalpha=0.05, strictmod=TRUE,
###'                         secondcomfun = "wilcox.test",
###'                         subclmin=3, subclwilc=TRUE,
###'                         secondalpha=0.01, lda=3)
###' restab <- as.data.frame(diffres)
###' head(restab)
##base::as.data.frame

#' @method as.data.frame diffAnalysisClass
#' @rdname as.data.frame
#' @export
as.data.frame.diffAnalysisClass <- function(x,...){
    efres <- tidyEffectSize(x)
    kwres <- x@kwres
    difftb <- merge(efres, kwres, by.x="f", by.y="f")
    difftb <- difftb[order(difftb$pvalue),]
    return(difftb)
}

#' @method as.data.frame alphasample
#' @rdname as.data.frame
#' @export
as.data.frame.alphasample <- function(x, ...){
    dat <- x@alpha
    if (!is.null(x@sampleda)){
        dat <- merge(dat, x@sampleda, by=0)
        rownames(dat) <- as.vector(dat$Row.names)
        dat$Row.names <- NULL
    }
    return(dat)
}

#' @keywords internal
tidyEffectSize <- function(obj){
    f <- LDA <- MeanDecreaseAccuracy <- NULL
    secondvars <- getsecondTRUEvar(obj)
    classname <- getcall(obj, "class")
    efres <- merge(obj@mlres, secondvars, by.x="f", by.y="f") %>%
             select (-c("gfc", "Freq"))
    if ("LDA" %in% colnames(efres)){
        efres <- efres %>% mutate(f = factor(f, levels=f[order(eval(parse(text=classname)), LDA)]))
    }else{
        efres <- efres %>% mutate(f = factor(f, levels=f[order(eval(parse(text=classname)),
                                                         MeanDecreaseAccuracy)]))
    }
    return(efres)
}

#' @keywords internal
getsecondTRUEvar <- function(obj){
    secondvars <- do.call("rbind",c(obj@secondvars,make.row.names=FALSE))
    secondvars <- secondvars %>% filter(eval(parse(text="gfc"))%in%"TRUE")
    return(secondvars)
}

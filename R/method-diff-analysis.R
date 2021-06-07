#' @title Differential expression analysis
#' @param obj object,a phyloseq class contained otu_table, sample_data, taxda, 
#' or data.frame, nrow sample * ncol features.
#' @param sampleda data.frame, nrow sample * ncol factor, the sample names of 
#' sampleda and data should be the same.
#' @param classgroup character, the factor name in sampleda.
#' @param subclass character, the factor name in sampleda, default is NULL, 
#' meaning no subclass compare.
#' @param taxda data.frame, the classification of the feature in data. 
#' default is NULL.
#' @param alltax logical, whether to set all classification as features if taxda is not NULL, 
#' default is TRUE.
#' @param standard_method character, the method of standardization, 
#' see also \code{\link[vegan]{decostand}}, default is NULL, 
#' it represents that the relative abundance of taxonomy will be used. If count was set,
#' it represents the count reads of taxonomy will be used.
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
#' @param clmin integer, the minimum number of samples per classgroup for performing test, default is 5.
#' @param clwilc logical, whether to perform test of per classgroup, default is TRUE.
#' @param secondalpha numeric, the alpha value for the second test, default is 0.05.
#' @param subclmin integer, the minimum number of samples per suclass for performing test, default is 3.
#' @param subclwilc logical, whether to perform test of per subclass, default is TRUE, meaning more strict. 
#' @param ldascore numeric, the threshold on the absolute value of the logarithmic LDA score, default is 2.
#' @param normalization integer, set the normalization value, set a big number if to get more meaningful values 
#' for the LDA score, or you can set NULL for no normalization, default is 1000000.
#' @param bootnums integer, set the number of bootstrap iteration for lda or rf, default is 30.
#' @param ci numeric, the confidence interval of effect size (LDA or MDA), default is 0.95.
#' @param type character, the type of datasets, default is "species", if the dataset is not about species, 
#' such as dataset of kegg function, you should set it to "others".
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
#' diffres <- diff_analysis(kostic2012crc, classgroup="DIAGNOSIS",
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
diff_analysis.data.frame <- function(obj, sampleda, classgroup, subclass=NULL, taxda=NULL,alltax=TRUE, standard_method=NULL, mlfun="lda", 
    ratio=0.7, firstcomfun='kruskal.test', padjust="fdr",filtermod="pvalue",
    firstalpha=0.05, strictmod=TRUE, fcfun="generalizedFC", secondcomfun="wilcox.test",
    clmin=5, clwilc=TRUE, secondalpha=0.05, subclmin=3, subclwilc=TRUE,	ldascore=2,
    normalization=1000000, bootnums=30, ci=0.95, type="species", ...){
    params <- list(...)
    if (!is.null(params$class) && inherits(class,"character")){
        message("The class argument has been deprecated. Please use `classgroup` instead!")
        classgroup <- params$class
    }
    if (!is.null(taxda)){
        if (!"fillNA" %in% names(attributes(taxda))){
            taxda <- fillNAtax(taxda, type=type)
        }
        if (alltax){obj <- get_alltaxdf(obj, taxda, method=standard_method)}
    }
    sampleda <- sampleda %>% select(c(classgroup, subclass))
    if (ncol(sampleda)>1){sampleda <- duplicatedtaxcheck(sampleda) %>% column_to_rownames(var="rowname")}
    vars <- colnames(obj)
    datameta <- merge(obj, sampleda, by=0)
    if (!is.factor(datameta[,match(x=classgroup, colnames(datameta))])){
        datameta[,match(x=classgroup, colnames(datameta))]<- as.factor(datameta[,match(x=classgroup, colnames(datameta))])}
    kwres <- multi_compare(fun=firstcomfun, data=datameta, feature=vars, factorNames=classgroup)
    kwres <- lapply(kwres, function(x)get_pvalue(x))
    kwres <- do.call("rbind", kwres)
    rownames(kwres) <- vars
    kwres <- data.frame(f=rownames(kwres),pvalue=kwres[,1])
    kwres$fdr <- p.adjust(kwres$pvalue, method=padjust)
    if (!filtermod=="pvalue"){varsfirst <- as.vector(kwres[kwres$fdr<=firstalpha& !is.na(kwres$fdr),,drop=FALSE]$f)
    }else{varsfirst <- as.vector(kwres[kwres$pvalue<=firstalpha&!is.na(kwres$pvalue),,drop=FALSE]$f)}
    if (!length(varsfirst)>0){stop(paste0("There are not significantly discriminative features before internal", secondcomfun," !"))}
    classlevels <- get_classlevels(sampleda, classgroup)
    compareclass <- get_compareclass(classlevels)
    if (!is.null(subclass) && strictmod){
        class2sub <- get_class2sub(sampleda, classgroup, subclass)
        comsubclass <- apply(compareclass,1,function(x)get_comparesubclass(x[1],x[2],class2sub))
        secondvars <- diffsubclass(datasample=datameta, features=varsfirst, comsubclass=comsubclass,classgroup=classgroup,subclass=subclass,
                                   fcfun=fcfun, secondcomfun=secondcomfun, submin=subclmin, subclwilc=subclwilc, pfold=secondalpha, ...)
    }else{
        secondvars <- diffclass(datasample=datameta, features=varsfirst, comclass=compareclass, classgroup=classgroup, fcfun=fcfun,
                                secondcomfun=secondcomfun,classmin=clmin,clwilc=clwilc,pfold=secondalpha, ...)
    }
    if (!length(secondvars)>0){stop(paste0("There are not significantly discriminative features after internal", secondcomfun," !"))}
    leaveclasslevels <- unlist(lapply(names(secondvars), function(x){unlist(strsplit(x,"-vs-"))}))
    secondvars <- get_consistentfeatures(diffsubclassfeature=secondvars, classgroup=classgroup,classlevels=leaveclasslevels)
    secondvarsvectors <- get_secondvarlist(secondvars=secondvars)
    if (!is.null(normalization)){obj <- obj * normalization}
    dameta <- merge(obj, sampleda, by=0) %>% column_to_rownames(var="Row.names")
    dameta <- dameta %>% select(c(secondvarsvectors, classgroup))
    dameta <- split(dameta, dameta[,match(classgroup,colnames(dameta))])
    dameta <- get_sampledflist(dameta, bootnums=bootnums, ratio=ratio)
    dameta <- remove_constant(dameta)
    if (mlfun=="lda"){mlres <- LDAeffectsize(dameta, compareclass, classgroup, bootnums=bootnums, LDA=ldascore, ci=ci)}
    if (mlfun=="rf"){mlres <- rfimportance(dameta, classgroup, bootnums=bootnums, effsize=ldascore, ci=ci)}
    params <- list(mlfun=mlfun, firstcomfun=firstcomfun, secondcomfun=secondcomfun,
                   firstalpha=firstalpha, filtermod=filtermod, classgroup=classgroup,
                   normalization=normalization, type=type, standard_method=standard_method,
                   subclass=subclass, strictmod=strictmod, fcfun=fcfun, clmin=clmin, clwilc=clwilc,
                   subclmin=subclmin, subclwilc=subclwilc)
    result <- merge_total_res(kwres=kwres, secondvars=secondvars, mlres=mlres, params=params)
    res <- new("diffAnalysisClass", originalD=obj, sampleda=sampleda, taxda=taxda, result=result, kwres=kwres,
               secondvars=secondvars, mlres=mlres, someparams=params)
    return(res)
}


#' @method diff_analysis phyloseq
#' @rdname diff_analysis
#' @importFrom phyloseq tax_table
#' @export
diff_analysis.phyloseq <- function(obj, ...){
    otuda <- checkotu(obj)
    sampleda <- checksample(obj)
    taxda <- obj@tax_table
    res <- diff_analysis.data.frame(obj=otuda,
                         sampleda=sampleda,
                         taxda=taxda,
                         ...)
    return(res)
}

#
#' @keywords internal 
#normalize_da <- function(data, method){
#    if (all.equal(data, round(data))){
#        if (is.null(method)){
#	    data <- apply(data, 1, function(x)100*x/sum(x))
#        }else{
#            data <- decostand(data, method=method, ...)
#        }
#        data <- data.frame(data, check.names=FALSE)
#    }
#    return(data)
#}

#' @title get the table of abundance of all level taxonomy 
#'
#' @description
#' This function was designed to get the abundance of all level taxonomy,
#' the input can be phyloseq object or data.frame.
#' @param obj object, phyloseq or data.frame
#' @param method character, the normalization method, 
#' see also \code{\link[vegan]{decostand}}, default is NULL, the relative abundance 
#' will be return, if it set `count`, the count table will be return.
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
setMethod("get_alltaxadf", "phyloseq",function(obj, method=NULL, type="species", ...){
    otuda <- checkotu(obj)
    if (is.null(obj@tax_table)){
        stop("The taxaonomy table is empty!")
    }else{
        taxa <- tax_table(obj)
        if (!"fillNA" %in% names(attributes(taxa))){
            taxa <- fillNAtax(taxa, type = type)
        }
    }
    data <- get_alltaxdf(data=otuda, taxda=taxa, taxa_are_rows=FALSE, method=method, ...)
    return(data)
})

#' @aliases get_alltaxadf,data.frame
#' @rdname get_alltaxadf
#' @export
setMethod("get_alltaxadf", "data.frame", function(obj, taxda, taxa_are_rows=FALSE, method=NULL, type="species", ...){
    #if (taxa_are_rows){
    #    obj <- data.frame(t(obj), check.names=FALSE)
    #}
    if (!"fillNA" %in% names(attributes(taxda))){
        taxda <- fillNAtax(taxda, type=type)
    }
    data <- get_alltaxdf(data=obj, taxda=taxda, taxa_are_rows=taxa_are_rows, method=method, ...)
    return(data)
})

#' @importFrom magrittr %>%
#' @keywords internal
get_alltaxdf <- function(data, taxda, method=NULL, taxa_are_rows=FALSE, ...){
    if (!taxa_are_rows){
        data <- data.frame(t(data), check.names=FALSE)
    }
    dt <- list()
    for (i in seq_len(ncol(taxda))){
        if (is.null(method)){
            dat <- get_ratio(data, taxda[,i,drop=FALSE])
        }else{
            dat <- get_count(data, taxda[,i,drop=FALSE])
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
###' diffres <- diff_analysis(kostic2012crc, classgroup="DIAGNOSIS",
###'                         mlfun="lda", filtermod="fdr",
###'                         firstcomfun = "kruskal.test",
###'                         firstalpha=0.05, strictmod=TRUE,
###'                         secondcomfun = "wilcox.test",
###'                         subclmin=3, subclwilc=TRUE,
###'                         secondalpha=0.01, lda=3)
###' restab <- as.data.frame(diffres)
###' head(restab)
##base::as.data.frame

## #' @method as.data.frame diffAnalysisClass
## #' @rdname as.data.frame
## #' @export
## as.data.frame.diffAnalysisClass <- function(x,...){
##     efres <- tidyEffectSize(x)
##     kwres <- x@kwres
##     difftb <- merge(efres, kwres, by.x="f", by.y="f")
##     difftb <- difftb[order(difftb$pvalue),]
##     return(difftb)
## }
## 
## #' @method as.data.frame alphasample
## #' @rdname as.data.frame
## #' @export
## as.data.frame.alphasample <- function(x, ...){
##     dat <- x@alpha
##     if (!is.null(x@sampleda)){
##         dat <- merge(dat, x@sampleda, by=0)
##         rownames(dat) <- as.vector(dat$Row.names)
##         dat$Row.names <- NULL
##     }
##     return(dat)
## }
## 
## #' @keywords internal
## tidyEffectSize <- function(obj){
##     f <- LDAmean <- MDAmean <- NULL
##     secondvars <- get_second_true_var(obj)
##     classname <- extract_args(obj, "classgroup")
##     efres <- merge(obj@mlres, secondvars, by.x="f", by.y="f") %>%
##              select (-c("gfc", "Freq"))
##     if ("LDAmean" %in% colnames(efres)){
##         efres <- efres %>% mutate(f = factor(f, levels=f[order(eval(parse(text=classname)), LDAmean)]))
##     }else{
##         efres <- efres %>% mutate(f = factor(f, levels=f[order(eval(parse(text=classname)),
##                                                          MDAmean)]))
##     }
##     return(efres)
## }
## 
## #' @keywords internal
## get_second_true_var <- function(obj){
##     secondvars <- do.call("rbind",c(obj@secondvars,make.row.names=FALSE))
##     secondvars <- secondvars %>% dplyr::filter(eval(parse(text="gfc"))%in%"TRUE")
##     return(secondvars)
## }

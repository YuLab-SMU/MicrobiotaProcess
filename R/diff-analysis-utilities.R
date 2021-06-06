#' @title a container for performing two or more sample test.
#' @param fun character, the method for test, optional ""
#' @param data data.frame, nrow sample * ncol feature+factorNames.
#' @param feature vector, the features wanted to test.
#' @param factorNames character, the name of a factor giving the corresponding groups.
#' @param subgroup vector, the names of groups, default is NULL.
#' @param ..., additional arguments for fun.
#' @return the result of fun, if fun is wilcox.test, it 
#' will return the list with class "htest".
#' @importFrom rlang new_formula
#' @importFrom stats wilcox.test
#' @author Shuangbin Xu
#' @export
#' @examples
#' datest <- data.frame(A=rnorm(1:10,mean=5),
#'                      B=rnorm(2:11, mean=6), 
#'                      group=c(rep("case",5),rep("control",5)))
#' head(datest)
#' multi_compare(fun=wilcox.test,data=datest,
#'               feature=c("A", "B"),factorNames="group")
#' da2 <- data.frame(A=rnorm(1:15,mean=5),
#'                   B=rnorm(2:16,mean=6),
#'                   group=c(rep("case1",5),rep("case2",5),rep("control",5)))
#' multi_compare(fun=wilcox.test,data=da2,
#'               feature=c("A", "B"),factorNames="group",
#'               subgroup=c("case1", "case2"))
multi_compare <- function(fun = wilcox.test, 
                          data, 
                          feature, 
                          factorNames, 
                          subgroup=NULL, ...){ 
    if (!is.null(subgroup)){
        data <- data[data[[factorNames]] %in% subgroup, ,drop=FALSE]
        data[[factorNames]] <- factor(data[[factorNames]], levels=subgroup)
    }
    lapply(feature,
           function(x){
           tmpformula <- new_formula(as.name(x), as.name(factorNames))
           suppressWarnings(do.call(fun,list(tmpformula,data=data, ...)))}) 
}

#' @title Methods for computation of the p-value
#' @param obj object, such as htest, lm, negbin
#' ScalarIndependenceTest class.
#' @return pvalue.
#' @author Shuangbin Xu
#' @export
#' @examples
#' library(nlme)
#' lmeres <- lme(distance ~ Sex,data=Orthodont)
#' pvalue <- get_pvalue(lmeres)
get_pvalue <- function(obj){
    UseMethod("get_pvalue")
}

#' @method get_pvalue htest
#' @rdname get_pvalue
#' @export
get_pvalue.htest <- function(obj){
    pvalue <- obj$p.value
    return(pvalue)
}

#' @method get_pvalue lme
#' @rdname get_pvalue
#' @export
get_pvalue.lme <- function(obj){
    anres <- get_anova(obj)
    pvalue <- anres[2,4]
    return (pvalue)
}

#' @method get_pvalue negbin
#' @rdname get_pvalue
#' @export
get_pvalue.negbin <- function(obj){
    anres <- get_anova(obj)
    pvalue <- anres[2,5]
    return (pvalue)
}

#' @method get_pvalue ScalarIndependenceTest
#' @rdname get_pvalue
#' @importFrom coin pvalue
#' @export
get_pvalue.ScalarIndependenceTest <- function(obj){
    pvalue(obj)
}

#' @method get_pvalue QuadTypeIndependenceTest
#' @rdname get_pvalue
#' @export
get_pvalue.QuadTypeIndependenceTest <- function(obj){
    get_pvalue.ScalarIndependenceTest(obj)
}

#' @method get_pvalue lm
#' @rdname get_pvalue
#' @export
get_pvalue.lm <- function(obj){
    anres <- get_anova(obj)
    pvalue <- anres[1,5]
    return(pvalue)
}

#' @method get_pvalue glm
#' @rdname get_pvalue
#' @importFrom stats pnorm
#' @export
get_pvalue.glm <- function(obj){
    anres <- summary(obj)
    pvalue <- 2*pnorm(abs(anres$coeff[2,3]),
					  lower.tail=FALSE)
	return(pvalue)
}

#' @importFrom stats anova
#' @keywords internal 
get_anova <- function(obj){
    anres <- suppressWarnings(anova(obj))
    return (anres)
}

#' @keywords internal
get_classlevels <- function(sampleda, classgroup){
    levelstmp <- unique(as.vector(sampleda[,match(classgroup, colnames(sampleda))]))
    return(levelstmp)
}

#' @keywords internal
get_class2sub <- function(sampleda, classgroup, subclass){
    tmpsplit <- sampleda[,match(classgroup, colnames(sampleda))]
    if (!missing(subclass)){
    	samplelist <- split(sampleda, tmpsplit) 
    	lapply(samplelist, function(x)unique(as.vector(x[,match(subclass, colnames(x))])))
    }
}

#' @importFrom gtools combinations
#' @keywords internal
get_compareclass <- function(classlevels){
    combinations(n=length(classlevels), r=2, v=classlevels,
    						  repeats.allowed=FALSE)
}

#' @keywords internal
get_comparesubclass <- function(xlevel, ylevel, class2sublist){
    res <- expand.grid(class2sublist[[xlevel]],
    				   class2sublist[[ylevel]])
    colnames(res) <- c(xlevel, ylevel)
    return(res)
}

#' @keywords internal
diffclass <- function(datasample,
                      features,
                      comclass,
                      classgroup,
                      fcfun="generalizedFC",
                      secondcomfun="wilcox.test",
                      classmin=3,
                      clwilc=TRUE,
                      pfold=0.05,
                      ...){
    keepfeature <- list()
    for (i in seq_len(nrow(comclass))){
        classtmp <- as.vector(comclass[i,])
        #clsize <- min(table(datasample[[class]]))
        datatmp <- datasample %>% dplyr::filter(eval(parse(text=classgroup)) %in% classtmp)
        datatmp[[match(classgroup, colnames(datatmp))]] <- factor(datatmp[[match(classgroup,colnames(datatmp))]], 
                                                             levels=classtmp)
        clsize <- min(table(datatmp[[match(classgroup,colnames(datatmp))]]))
        resgfoldC <- get_gfc_wilc(datasample=datatmp,
                                fun1=fcfun,
                                classlevelsnum=clsize,
                                vars=features,
                                classname=classgroup,
                                minnum=classmin,
                                fun2=secondcomfun, 
                                wilc=clwilc,
                                ...)
        resgfoldC <- get_compareres(resgfoldC, pfold=pfold)
        keepfeature[[paste(classtmp, collapse="-vs-")]] <- resgfoldC[resgfoldC$Freq==1,]
    }
    return(keepfeature)
}

#' @keywords internal
diffsubclass <- function(datasample,
                         features,
                         comsubclass, 
                         classgroup, 
                         subclass, 
                         fcfun="generalizedFC",
                         secondcomfun="wilcox.test",
                         submin=3,
                         subclwilc=TRUE,
                         pfold=0.05,
                         ...){
    keepfeature <- list()
    for (i in seq_len(length(comsubclass))){
        classtmp <- colnames(comsubclass[[i]])
        reslist <- list()
        for (j in seq_len(nrow(comsubclass[[i]]))){
            subclasstmp <- as.vector(unlist(comsubclass[[i]][j,])) 
            datatmp <- datasample %>% 
                       dplyr::filter(eval(parse(text=classgroup)) %in% classtmp & eval(parse(text=subclass)) %in%subclasstmp)
            datatmp[[match(subclass, colnames(datatmp))]] <- factor(datatmp[[match(subclass,colnames(datatmp))]],
                                                                             levels=subclasstmp)
            subclsize <- min(table(datatmp[[match(subclass,colnames(datatmp))]]))
            resgfoldC <- get_gfc_wilc(datasample=datatmp, 
                                    classlevelsnum=subclsize,
                                    fun1=fcfun,
                                    vars=features,
                                    classname=subclass,
                                    minnum=submin,
                                    fun2=secondcomfun,
                                    wilc=subclwilc,
                                    ...)
            reslist[[j]] <- resgfoldC
        }
        reslist <- do.call("rbind", reslist)
        reslist <- get_compareres(reslist, pfold=pfold)
        reslist <- reslist[reslist$Freq==nrow(comsubclass[[i]]),]
        if (nrow(reslist)>0){
            keepfeature[[paste(classtmp, collapse="-vs-")]] <- reslist
        }
    }
    return(keepfeature)
}

#' @keywords internal
get_gfc_wilc <- function(datasample, classlevelsnum, fun1='generalizedFC', 
                       vars, classname, minnum, fun2='wilcox.test', wilc,...){
    resgfoldC <- multi_compare(fun=fun1, data=datasample,
                               feature=vars, factorNames=classname)
    resgfoldC <- switch(fun1,
                        generalizedFC=lapply(resgfoldC, function(x)x$gfc),
                        compare_mean = lapply(resgfoldC, function(x)x$diffmean),
                        compare_median = lapply(resgfoldC, function(x)x$diffmedian)
                        )
    resgfoldC <- do.call("rbind", resgfoldC)
    rownames(resgfoldC) <- vars
    if (classlevelsnum>= minnum &&  wilc){
        tmpres <- multi_compare(fun=fun2, data=datasample,
                                feature=vars, factorNames=classname, ...)
        pvaluetmp <- lapply(tmpres, function(x)get_pvalue(x))
        pvaluetmp <- do.call("rbind", pvaluetmp)
        rownames(pvaluetmp) <- vars
        resgfoldC <- merge(resgfoldC, pvaluetmp, by=0)
        colnames(resgfoldC) <- c("f", "gfc", "pvalue")
    }else{
        resgfoldC <- data.frame(f=rownames(resgfoldC), gfc=resgfoldC[,1], pvalue=0)
    }
    return(resgfoldC)
}

#' @keywords internal
get_compareres <- function(reslist, pfold){
    if (ncol(reslist)<3){
        reslist <- data.frame(f=rownames(reslist),gfc=reslist[,1], stringsAsFactors =FALSE)
        reslist$gfc <- reslist$gfc >0
        reslist <- data.frame(table(reslist), stringsAsFactors =FALSE)
    }else{
        reslist <- data.frame(reslist, stringsAsFactors =FALSE)
        colnames(reslist) <- c("f", "gfc", "pvalue")
        reslist <- reslist[reslist$pvalue <= pfold &!is.na(reslist$pvalue),,drop=FALSE]
        reslist$pvalue <- NULL
        reslist$gfc <- reslist$gfc > 0 
        reslist <- data.frame(table(reslist),stringsAsFactors =FALSE)
    }
    return(reslist)
}

#' @keywords internal
get_consistentfeatures <- function(diffsubclassfeature, 
                                  classgroup,
                                  classlevels, ...){
    leftfeature <- list()
    for (i in classlevels){
        tmpindex <- grep(i, names(diffsubclassfeature))
        tmpkeepfeature <- diffsubclassfeature[tmpindex]
        checkflag <- grepl(paste0("^",i, "-vs-"), names(tmpkeepfeature))
        falseindex <- which(!checkflag) 
        if (length(falseindex)>0){
            for (j in falseindex){
                tmpkeepfeature[[j]]$gfc <- !as.logical(tmpkeepfeature[[j]]$gfc)
            }
        }
        tmpkeepfeature <- do.call("rbind", tmpkeepfeature) 
        tmpkeepfeature <- data.frame(table(tmpkeepfeature[,colnames(tmpkeepfeature) %in% c("f", "gfc")]))
        tmpkeepfeature <- tmpkeepfeature[tmpkeepfeature$Freq==length(tmpindex),]
        if(nrow(tmpkeepfeature)>0){
            tmpkeepfeature[[classgroup]] <- i
            leftfeature[[i]] <- tmpkeepfeature
        }
    }
    return(leftfeature)
}

#' @keywords internal
get_secondvarlist <- function(secondvars){
    vars <- as.vector(unique(unlist(lapply(secondvars, 
                     function(x){as.vector(x$f)}))))
    return(vars)
}


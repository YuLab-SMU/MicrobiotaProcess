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
#' multi.compare(fun=wilcox.test,data=datest,
#'               feature=c("A", "B"),factorNames="group")
#' da2 <- data.frame(A=rnorm(1:15,mean=5),
#'                   B=rnorm(2:16,mean=6),
#'                   group=c(rep("case1",5),rep("case2",5),rep("control",5)))
#' multi.compare(fun=wilcox.test,data=da2,
#'               feature=c("A", "B"),factorNames="group",
#'               subgroup=c("case1", "case2"))
multi.compare <- function(fun = wilcox.test, 
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
#' pvalue <- getpvalue(lmeres)
getpvalue <- function(obj){
    UseMethod("getpvalue")
}

#' @method getpvalue htest
#' @rdname getpvalue
#' @export
getpvalue.htest <- function(obj){
    pvalue <- obj$p.value
    return(pvalue)
}

#' @method getpvalue lme
#' @rdname getpvalue
#' @export
getpvalue.lme <- function(obj){
    anres <- getanova(obj)
    pvalue <- anres[2,4]
    return (pvalue)
}

#' @method getpvalue negbin
#' @rdname getpvalue
#' @export
getpvalue.negbin <- function(obj){
    anres <- getanova(obj)
    pvalue <- anres[2,5]
    return (pvalue)
}

#' @method getpvalue ScalarIndependenceTest
#' @rdname getpvalue
#' @importFrom coin pvalue
#' @export
getpvalue.ScalarIndependenceTest <- function(obj){
    pvalue(obj)
}

#' @method getpvalue QuadTypeIndependenceTest
#' @rdname getpvalue
#' @export
getpvalue.QuadTypeIndependenceTest <- function(obj){
    getpvalue.ScalarIndependenceTest(obj)
}

#' @method getpvalue lm
#' @rdname getpvalue
#' @export
getpvalue.lm <- function(obj){
    anres <- getanova(obj)
    pvalue <- anres[1,5]
    return(pvalue)
}

#' @method getpvalue glm
#' @rdname getpvalue
#' @importFrom stats pnorm
#' @export
getpvalue.glm <- function(obj){
    anres <- summary(obj)
    pvalue <- 2*pnorm(abs(anres$coeff[2,3]),
					  lower.tail=FALSE)
	return(pvalue)
}

#' @importFrom stats anova
#' @keywords internal 
getanova <- function(obj){
    anres <- suppressWarnings(anova(obj))
    return (anres)
}

#' @keywords internal
getclasslevels <- function(sampleda, class){
    levelstmp <- unique(as.vector(sampleda[,match(class, colnames(sampleda))]))
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

#' @importFrom gtools combinations
#' @keywords internal
getcompareclass <- function(classlevels){
    combinations(n=length(classlevels), r=2, v=classlevels,
    						  repeats.allowed=FALSE)
}

#' @keywords internal
getcomparesubclass <- function(xlevel, ylevel, class2sublist){
    res <- expand.grid(class2sublist[[xlevel]],
    				   class2sublist[[ylevel]])
    colnames(res) <- c(xlevel, ylevel)
    return(res)
}

#' @keywords internal
diffclass <- function(datasample,
                      features,
                      comclass,
                      class,
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
        datatmp <- datasample %>% filter(eval(parse(text=class)) %in% classtmp)
        datatmp[[match(class, colnames(datatmp))]] <- factor(datatmp[[match(class,colnames(datatmp))]], 
                                                             levels=classtmp)
        clsize <- min(table(datatmp[[match(class,colnames(datatmp))]]))
        resgfoldC <- getgfcwilc(datasample=datatmp,
                                fun1=fcfun,
                                classlevelsnum=clsize,
                                vars=features,
                                classname=class,
                                minnum=classmin,
                                fun2=secondcomfun, 
                                wilc=clwilc,
                                ...)
        resgfoldC <- getcompareres(resgfoldC, pfold=pfold)
        keepfeature[[paste(classtmp, collapse="-vs-")]] <- resgfoldC[resgfoldC$Freq==1,]
    }
    return(keepfeature)
}

#' @keywords internal
diffsubclass <- function(datasample,
                         features,
                         comsubclass, 
                         class, 
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
                       filter(eval(parse(text=class)) %in% classtmp & eval(parse(text=subclass)) %in%subclasstmp)
            datatmp[[match(subclass, colnames(datatmp))]] <- factor(datatmp[[match(subclass,colnames(datatmp))]],
                                                                             levels=subclasstmp)
            subclsize <- min(table(datatmp[[match(subclass,colnames(datatmp))]]))
            resgfoldC <- getgfcwilc(datasample=datatmp, 
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
        reslist <- getcompareres(reslist, pfold=pfold)
        reslist <- reslist[reslist$Freq==nrow(comsubclass[[i]]),]
        if (nrow(reslist)>0){
            keepfeature[[paste(classtmp, collapse="-vs-")]] <- reslist
        }
    }
    return(keepfeature)
}

#' @keywords internal
getgfcwilc <- function(datasample, classlevelsnum, fun1='generalizedFC', 
                       vars, classname, minnum, fun2='wilcox.test', wilc,...){
    resgfoldC <- multi.compare(fun=fun1, data=datasample,
                               feature=vars, factorNames=classname)
    resgfoldC <- lapply(resgfoldC, function(x)x$gfc)
    resgfoldC <- do.call("rbind", resgfoldC)
    rownames(resgfoldC) <- vars
    if (classlevelsnum>= minnum &&  wilc){
        tmpres <- multi.compare(fun=fun2, data=datasample,
                                feature=vars, factorNames=classname, ...)
        pvaluetmp <- lapply(tmpres, function(x)getpvalue(x))
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
getcompareres <- function(reslist, pfold){
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
getconsistentfeatures <- function(diffsubclassfeature, 
                                  class,
                                  classlevels, ...){
    leftfeature <- list()
    for (i in classlevels){
        tmpindex <- grep(i, names(diffsubclassfeature))
        tmpkeepfeature <- diffsubclassfeature[tmpindex]
        checkflag <- grepl(paste0(i, "-vs-"), names(tmpkeepfeature))
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
            tmpkeepfeature[[class]] <- i
            leftfeature[[i]] <- tmpkeepfeature
        }
    }
    return(leftfeature)
}

#' @keywords internal
getsecondvarlist <- function(secondvars){
    vars <- as.vector(unique(unlist(lapply(secondvars, 
                     function(x){as.vector(x$f)}))))
    return(vars)
}


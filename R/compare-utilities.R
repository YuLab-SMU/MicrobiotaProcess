#' @title a container for performing two or more sample test.
#' @param fun character, the method for test, optional ""
#' @param data data.frame, nrow sample * ncol feature+factorNames.
#' @param feature vector, the features wanted to test.
#' @param factorNames character, the name of a factor giving the corresponding groups.
#' @param ..., additional arguments for fun.
#' @author ShuangbinXu
#' @export
multi.compare <- function(fun = wilcox.test, 
						data, 
						feature, 
						factorNames, ...){ 
	sapply(feature,
		   function(x){
		   tmpformula <- as.formula(paste(x, factorNames, sep="~"))
		   do.call(fun,list(tmpformula,data=data, ...))}, 
		   simplify = FALSE)
}	 

#' @keywords internal
getclasslevels <- function(sampleda, class){
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

#' @importFrom gtools combinations
#' @keywords internal
getcompareclass <- function(classlevels){
	combinations(length(classlevels), 2, classlevels,
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
					  classmin=5,
					  clwil=TRUE,
					  pfold=0.05,
					  ...){
	keepfeature <- list()
	for (i in 1:nrow(comclass)){
		classtmp <- as.vector(comclass[i,])
		datatmp <- datasample %>% filter(class %in% classtmp)
		datatmp[[class]] <- factor(datatmp[[class]], levels=classtmp)
		clsize <- table(datatmp[[class]])
		resgfoldC <- getgfcwilc(datasample=datatmp,
								classlevelsnum=clsize,
								vars=features,
								classname=class,
								minnum=classmin,
								wilc=clwil)
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
						 submin=5,
						 subclwil=TRUE,
						 pfold=0.05,
						 ...){
	keepfeature <- list()
	for (i in 1:length(comsubclass)){
		classtmp <- colnames(comsubclass[[i]])
		reslist <- list()
		for (j in 1:nrow(comsubclass[[i]])){
			subclasstmp <- as.vector(unlist(comsubclass[[i]][j,])) 
			datatmp <- datasample %>% filter(class %in% classtmp & subclass%in%subclasstmp)
			datatmp[[subclass]] <- factor(datatmp[[subclass]], levels=subclasstmp)
			subclsize <- table(datatmp[[subclass]])
			resgfoldC <- getgfcwilc(datasample=datatmp, 
									classlevelsnum=subclsize,
									vars=features,
									classname=subclass,
									minnum=submin,
									wilc=subclwil)
			}
			reslist[[j]] <- resgfoldC
		}
		reslist <- do.call("rbind", reslist)
		reslist <- getcompareres(retslist, pfold=pflod)
		keepfeature[[paste(classtmp, collapse="-vs-")]] <- reslist[reslist$Freq==nrow(comsubclass[[i]]),]
		return(keepfeature)
}

#' @keywords internal
getgfcwilc <- function(datasample, classlevelsnum, vars, classname, minnum, wilc){
	resgfoldC <- multi.compare(fun=generalizedFC, data=datasample,
							   feature=vars, factorNames=classname)
	resgfoldC <- lapply(resgfoldC, function(x)x$gfc)
	resgfoldC <- do.call("rbind", resgfoldC)
	if (min(classlevelsnum)>= minnum &&  subclwil){
		tmpres <- multi.compare(fun=wilcox.test, data=datasample,
								feature=features, factorNames=classname)
		pvaluetmp <- lapply(tmpres, function(x)x$p.value)
		pvaluetmp <- do.call("rbind", pvaluetmp)
		resgfoldC <- merge(resgfoldC, pvaluetmp, by=0)
	}
	return(resgfoldC)
}

#' @keywords internal
getcompareres <- function(relist, pfold){
	if (ncol(reslist)<3){
		reslist <- data.frame(f=rownames(reslist),gfc=reslist[,1], stringsAsFactors =FALSE)
		reslist$gfc <- reslist$gfc >0
		reslist <- data.frame(table(reslist), stringsAsFactors =FALSE)
	}else{
		reslist <- data.frame(reslist, stringsAsFactors =FALSE)
		colnames(reslist) <- c("f", "gfc", "pvalue")
		reslist <- reslist[reslist$pvalue < pflod,]
		reslist$pvalue <- NULL
		reslist$gfc <- reslist$gfc > 0 
		reslist <- data.frame(table(reslist),stringsAsFactors =FALSE)
	}
	return(reslist)
}

#' @keywords internal
getconsistentfeatures <- function(diffsubclassfeature, 
								  classlevels, ...){
	leftfeature <- list()
	for (i in classlevels){
		tmpindex <- grep(i, names(diffsubclassfeature))
		tmpkeepfeature <- keepfeature[tmpindex]
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
		leftfeature[[i]] <- tmpkeepfeature
	}
	return(leftfeature)
}


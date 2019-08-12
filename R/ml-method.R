#' @importFrom dplyr select
#' @importFrom MASS lda
#' @keywords internal
LDAeffectsize <- function(datalist, compareclass, class, bootnums=30, LDA=2){
	res <- list()
	tmpformula <- as.formula(paste0(class, " ~. "))
	for (i in seq_len(bootnums)){
		compareres <- list()
		df <- datalist[[i]]
		for (p in seq_len(nrow(compareclass))){
			tmppairs <- compareclass[p,]
			z <- suppressWarnings(lda(tmpformula,data=df, tol=1e-6))
			w <- z$scaling[,1]
			w.unit <- w/sqrt(sum(w^2))
			ss <- df %>% select(-c(class)) %>% as.matrix()
			LD <- ss %*% w.unit
			tmpp1 <- df[[match(class, colnames(df))]]==tmppairs[1]
			tmpp2 <- df[[match(class, colnames(df))]]==tmppairs[2]
			effect.size <- abs(mean(LD[tmpp1,])-mean(LD[tmpp2,]))
			wfinal <- w.unit * effect.size
			coeff <- abs(wfinal)
			mm <- z$means
			gm <- abs(mm[match(tmppairs[1], rownames(mm)),,drop=FALSE]- mm[match(tmppairs[2], rownames(mm)),,drop=FALSE])
			tmpres <- (gm+coeff) *0.5
			compareres[[p]] <- tmpres
		}
	compareres <- do.call("rbind", compareres)
	res[[i]] <- compareres
	}
	res <- Reduce("+", res)
	res <- apply(res/bootnums, 2, function(x)max(x))
	res <- log(1+res,10)
	res <- data.frame(f=names(res), LDA=res)
	res$f <- gsub("`", "", as.vector(res$f))
	res <- res[res$LDA>=LDA,]
	if (!nrow(res)>0){stop("No features with significant differences after LDA analysis.")}
	rownames(res) <- NULL
	return(res)
}

#' @importFrom randomForest randomForest
#' @keywords internal
rfimportance <- function(datalist, class, bootnums){
	rfres <- list()
	tmpformula <- as.formula(paste0(class, " ~. "))
	for (i in seq_len(bootnums)){
		df <- datalist[[i]]
		dfres <- randomForest(tmpformula, data=df, importance=TRUE, proximity=TRUE)
		imres <- dfres$importance	
		rfres[[i]] <- imres
	}
	rfres <- Reduce("+", rfres)
	rfres <- data.frame(rfres/bootnums)
	rfres <- data.frame(f=rownames(rfres), MeanDecreaseAccuracy=rfres$MeanDecreaseAccuracy)
	return(rfres)
}


#' @title Generate random data list from a original data.
#' @param dalist list, a list contained multi data.frame.
#' @param ratio numeric, the ratios of each data.frame to keep.
#' @param seednums interger, random number generation, default is 1024.
#' @param makerownames logical, whether build row.names,default is FALSE.
#' @param randomSeed integer, A single integer value passed to ‘set.seed’, 
#' which is used to fix a seed for reproducibly random number, default is 1024.
#' @author Shuangbin Xu
#' @export
sampledflist <- function(dalist, 
						 bootnums=30, 
						 ratio=0.7, 
						 makerownames=FALSE,
						 randomSeed=1024){
	datalist <- list()
	if (!is.null(randomSeed)){
		set.seed(randomSeed)
	}
	for (i in seq_len(bootnums)){
		res <- lapply(dalist,function(x){x[sample(nrow(x), trunc(nrow(x)*ratio)),,drop=FALSE]})
		res <- do.call("rbind", c(res, make.row.names=makerownames))
		datalist[[i]] <- res
	}
	return(datalist)
}


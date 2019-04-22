#' @title Rarefaction alpha index
#'
#' @description
#' building curve of alpha index in Rarefied samples.
#' @param data data.frame,(nrow sample * ncol feature (factor)) or
#' the data.frame for stat_smooth.
#' @param nrows, the nrow of facet, default is 2.
#' @param mapping, set of aesthetic mapping of ggplot2, default is NULL,
#' if the data is the data.frame for stat_smooth, the mapping should be set.
#' @param linesize integer, default is 0.5.
#' @param chunks integer, the number of subsample in a sample, 
#'  default is 400.
#' @param sampleda, data.frame, (nrow sample * ncol factor)
#' @param factorNames character, default is missing. 
#' @param factorlevels list, the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this.
#' @param facetNames vector character, default is NULL.
#' @param se boolean, default is FALSE.
#' @param method character, default is lm.
#' @param formula formula, default is `y ~ log(x)`
#' @author ShuangbinXu
#' @importFrom ggplot2 ggplot aes_ stat_smooth facet_wrap
#' @importFrom dplyr filter

ggrarecurve <- function(data,
			   nrows=2,
			   mapping=NULL,
			   linesize=0.5,	
			   chunks=400,
			   sampleda,
			   factorNames,
			   factorlevels,
			   facetNames=NULL,
			   se=FALSE,
			   method="lm",			   
			   formula=y ~ log(x),
			   ...){
	if (is.null(mapping)){
		data <- stat_rare(data,
						  chunks=chunks,
						  sampleda=sampleda,
						  factorlevels=factorlevels,
						  plotda=TRUE)
		mapping <- aes_(~readsNums, ~value, color=~sample)
		if (!missing(factorNames)){
		    	tmpcolor <- as.formula(paste0("~", factorNames))
			mapping <- modifyList(mapping,
						 aes_(~readsNums,
						      ~value, 
						      color=tmpcolor))
		}
	}
	if (!is.null(facetNames)){
		data <- data %>% filter(Alpha %in% facetNames)
	}
	p <- ggplot(data=data,
		     mapping=mapping) +
		 stat_smooth(se=se, 
		 		method = method,
		 		size=linesize,
		  		formula = formula,
				...) 
	if (!missing(nrows)){
		p <- p + facet_wrap(~ Alpha, 
				      scales="free", 
				      nrow=nrows)
	}
	return(p)
}

#' @param data data.frame,(nrow sample * ncol taxonomy 
#' (feature) or and factor)
#' @param chunks integer, the number of subsample in a sample,
#' default is 400.
#' @param sampleda data.frame, (nrow sample * ncol factor)
#' @param factorlevels list, the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this.
#' @param plotda boolean, default is TRUE, whether build the data.frame for
#' `geom_bar` of `ggplot2`.
#' @author ShuangbinXu
#' @importFrom dplyr bind_rows
#' @importFrom reshape melt
#' @importFrom magrittr %>%

stat_rare <- function(data, 
					  chunks=400, 
					  sampleda,
					  factorlevels,
					  plotda=TRUE){
	tmpfeature <- colnames(data)[sapply(data,is.numeric)]
    	tmpfactor <- colnames(data)[!sapply(data,is.numeric)]
	dat <- data[ , tmpfeature, drop=FALSE]
	out <- apply(dat, 1, samplealpha, chunks=chunks) %>% 
		bind_rows(,.id="sample")
	if (plotda){
		if (!missing(sampleda)){
			sampleda$sample <- rownames(sampleda)
			out <- merge(out, sampleda)
			out <- melt(out,id.vars=c(colnames(sampleda), "readsNums"),
						variable_name="Alpha")
		}
		if (missing(sampleda) && length(tmpfactor) > 0){
			tmpsample <- data[, tmpfactor, drop=FALSE]
			tmpsample$sample <- rownames(tmpsample)
			out <- merge(out, tmpsample)
			out <- melt(out, id.vars=c("sample", "readsNums", tmpfactor),
						variable_name="Alpha")
		}
	}else{
		if (!missing(sampleda)){
			sampleda$sample <- rownames(sampleda)
			out <- merge(out, sampleda)
		}
		if (missing(sampleda) && length(tmpfactor) >0){
			tmpsample <- data[,tmpfactor,drop=FALSE]
			tmpsample$sample <- rownames(tmpsample)
			out <- merge(out, tmpsample)
		}
	}
	if (!missing(factorlevels)){
		out <- setfactorlevels(out, factorlevels)
	}
	return(out)
}

# @keyword internal
samplealpha <- function(data, chunks=200){
	sdepth <- sum(data)
	step <- trunc(sdepth/chunks)
	n <- seq(0, sdepth, by=step)[-1]
	n <- c(n, sdepth)	
	out <- lapply(n, function(x){
			tmp <- alphaindex(data, mindepth=x)
			tmp$readsNums <- x
		    return(tmp)}) 
	out <- do.call("rbind", c(out, make.row.names=FALSE))
	out[is.na(out)] <- 0
	return (out)
}

#' @importFrom stats predict
#' @importFrom dplyr tibble
#' @export
# this is from ggplot2
predictdf.lm <- function(model, xseq, se, level) {
    pred <- predict(model, newdata = tibble(x = xseq), se.fit = se,
    level = level, interval = if (se) "confidence" else "none")
   
  if (se) {
    fit <- as.data.frame(pred$fit)
    names(fit) <- c("y", "ymin", "ymax")
    res <- data.frame(x = xseq, fit, se = pred$se.fit)
  } else {
    res <- data.frame(x = xseq, y = as.vector(pred))
  }
  # add the x=zero ,y=zero
  res <- rbind(rep(0, ncol(res)), res)
  return(res)
}  


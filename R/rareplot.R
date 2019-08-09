#' @method ggrarecurve default
#' @importFrom ggplot2 ggplot aes_ stat_smooth facet_wrap
#' @importFrom dplyr filter
#' @rdname ggrarecurve
#' @export
ggrarecurve.default <- function(data,
			   sampleda,
			   indexNames="Observe",
			   linesize=0.5,
		   	   facetnrow=1,
			   mapping=NULL,	   
			   chunks=400,
			   factorNames,
			   factorLevels,
			   se=FALSE,
			   method="lm",			   
			   formula=y ~ log(x),
			   ...){
	if (is.null(mapping)){
		data <- stat_rare(data,
						  chunks=chunks,
						  sampleda=sampleda,
						  factorLevels=factorLevels,
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
	if (!is.null(indexNames)){
		data <- data %>% filter(Alpha %in% indexNames)
	}
	p <- ggplot(data=data,
		     mapping=mapping) +
		 stat_smooth(se=se, 
		 		method = method,
		 		size=linesize,
		  		formula = formula,
				...) 
	#if (!missing(nrows)){
	p <- p + facet_wrap(~ Alpha, scales="free", nrow=facetnrow)
	#}else{
	#	p <- p + facet_wrap(~ Alpha, scales="free")
	#}
	return(p)
}

#' @title mapping data of ggrarecurve
#' @description
#' generating the data of ggrarecurve.
#' @param data data.frame,(nrow sample * ncol taxonomy 
#' (feature) or and factor)
#' @param chunks integer, the number of subsample in a sample,
#' default is 400.
#' @param sampleda data.frame, (nrow sample * ncol factor)
#' @param factorLevels list, the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this.
#' @param plotda boolean, default is TRUE, whether build the data.frame for
#' `geom_bar` of `ggplot2`.
#' @author ShuangbinXu
#' @importFrom dplyr bind_rows
#' @importFrom reshape melt
#' @importFrom magrittr %>%
#' @export
stat_rare <- function(data, 
					  chunks=400, 
					  sampleda,
					  factorLevels,
					  plotda=TRUE){
	tmpfeature <- colnames(data)[sapply(data,is.numeric)]
   	tmpfactor <- colnames(data)[!sapply(data,is.numeric)]
	dat <- data[ , match(tmpfeature, colnames(data)), drop=FALSE]
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
		if (missing(sampleda)&&length(tmpfactor) == 0){
			out <- melt(out, id.vars=c("sample", "readsNums"),
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
	if (!missing(factorLevels)){
		out <- setfactorlevels(out, factorLevels)
	}
	return(out)
}

#' @keywords internal
samplealpha <- function(data, chunks=200){
	sdepth <- sum(data)
	step <- trunc(sdepth/chunks)
	n <- seq(0, sdepth, by=step)[-1]
	n <- c(n, sdepth)
	out <- lapply(n, function(x){
			tmp <- alphaindex(data, mindepth=x)
			#tmp <- tmp$indexs
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


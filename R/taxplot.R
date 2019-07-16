#' @method ggbartax default
#' @rdname ggbartax
#' @importFrom ggplot2 ggplot aes_ geom_bar
#' @export
ggbartax.default <- function(data,
				  mapping=NULL,
				  position = "stack",
				  stat="identity",
				  width=0.7,
				  topn=30,
				  count=FALSE,
				  sampleda=NULL,
				  factorLevels=NULL,
				  settheme=TRUE,
				  facetNames=NULL,
				  setColors=TRUE,
				  ...){
	if (is.null(mapping)){
		mapping <- aes_(~sample, ~value, fill=~feature)
		data <- mappingtaxda(data, 
							 topn=topn, 
							 count=count,
							 sampleda=sampleda, 
							 factorLevels=factorLevels,
							 plotda=TRUE)
	}else{
		mapping <- mapping
	}
	p <- ggplot(data,
				mapping=mapping,
				...) + 
		geom_bar(position = position, 
				 stat=stat,
				 width=width) + 
		theme_bw() 
	tmpfactor <- setdiff(colnames(data), c("feature", "sample", "value"))
	#if (is.null(facetNames) && length(tmpfactor)>0){
	#	tmpformula <- as.formula(paste0("~ ",tmpfactor[1]))
	#	p <- p + facet_grid(tmpformula, scales="free_x", space="free_x")
	#}
	if(setColors){
		tmpn <- length(levels(data$feature))
		p <- p + scale_fill_manual(values=getCols(tmpn))
	}
	if (!is.null(facetNames)){
		tmpformula <- as.formula(paste0("~ ", facetNames))
		p <- p + facet_grid(tmpformula, scales="free_x", space="free_x")
	}
	if (settheme){
		p <- p + taxbartheme()
	}
	return(p)
}

#' @title generate the mapping data
#' @description
#' build the data.frame for `geom_bar` of `ggplot2`,
#' if the data.frame was built, you can visulize with `ggbartax`.
#' @param data data.frame, (nrow sample * ncol feature (factor))
#' @param topn integer, the top number of abundance taxonomy(feature).
#' @param count boolean, default is FALSE.
#' @param sampleda data.frame, (nrow sample * ncol factor) the 
#' sample information, if the data doesn't contain the information.
#' @param factorLevels list, the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this. 
#' @param plotda boolean, default is TRUE, whether build the data.frame for
#' `geom_bar` of `ggplot2`.
#' @author ShuangbinXu
#' @export
#' @importFrom magrittr %>%
#' @importFrom reshape melt

mappingtaxda <- function(data, topn=30, 
						 count=FALSE,
						 sampleda=NULL, 
						 factorLevels=NULL,
						 plotda=TRUE){
	tmpfeature <- colnames(data)[sapply(data,is.numeric)]
	tmpfactor <- colnames(data)[!sapply(data,is.numeric)]
	dat <- data[, tmpfeature,drop=FALSE] %>% t() %>% 
		data.frame(check.names=FALSE)
	if(!count){
		dat <- apply(dat, 2, function(x){100*x/sum(x)}) %>% 
			data.frame(check.names=FALSE)
	}
	dat$sums <- apply(dat, 1, sum)
	dat <- dat[order(dat$sum, decreasing = TRUE),,drop=FALSE]
	dat$sums <- NULL
	tmpsums <- matrix(colSums(dat),nrow=1) %>% data.frame()
	if (!missing(topn)){
		dat <- dat[1:topn,,drop=FALSE]
		if (!count){
			others <- 100 - (matrix(apply(dat,2,sum),nrow=1) %>% 
							 data.frame(check.names=FALSE))
		}
		if (count){
			others <- tmpsums - (matrix(apply(dat,2,sum),nrow=1) %>% 
							 data.frame(check.names=FALSE))
		}
		colnames(others) <- colnames(dat)
		rownames(others) <- "Others"
		dat <- rbind(dat, others)
	}
	featurelevels <- rownames(dat)
	if (plotda){
		dat <- melt(as.matrix(dat))
		colnames(dat) <- c("feature", "sample", "value")
		if (!is.null(sampleda)){
			sampleda$sample <- rownames(sampleda)
			dat <- merge(dat, sampleda)
		}
		if (is.null(sampleda) && length(tmpfactor)>0){
			tmpsample <- data[,tmpfactor,drop=FALSE]
			tmpsample$sample <- rownames(tmpsample)
			dat <- merge(dat, tmpsample)
		}
		dat$feature <- factor(dat$feature, levels=featurelevels)
	}else{
		if (!is.null(sampleda)){
			dat <- dat %>% t() %>% 
				data.frame(check.names=F)
			dat <- merge(dat, sampleda, by=0)
		}
		if (is.null(sampleda)&&length(tmpfactor)>0){
			tmpsample <- data[,tmpfactor,drop=FALSE]
			dat <- merge(dat, tmpsample, by=0)
		}
		colnames(dat)[1] <- "sample"
		
	}
	if (!is.null(factorLevels)){
		dat <- setfactorlevels(dat, factorLevels)	
	}
	return(dat)
}


#setfactorlevels <- function(data, factorlist){
#	factornames <- intersect(colnames(dat), names(factorlist))
#	if (length(factornames)>0){
#		for(i in factornames){
#				   data[[i]] <- factor(data[[i]], 
#									 levels=as.vector(factorlist[[i]]))
#		}
#	}
#	return(data)
#}


#' @importFrom ggplot2 theme
#' @keywords internal
taxbartheme <- function(){
	theme(axis.text.x = element_text(angle = -45, hjust = 0, size=7),
		  panel.grid = element_blank(),
		  legend.position = "bottom", 
		  legend.box = "horizontal", 
		  legend.text = element_text(size = 8), 
		  legend.title=element_blank(),
		  plot.margin=unit(c(0.2,1,0.2,0.2),"cm"),
	      panel.spacing = unit(0.2, "mm"),
	      strip.text.x = element_text(size=12, face="bold"),
	      strip.background = element_rect(colour="white", fill="grey"))
}

#' @title legend guides for ggbartax
#' @description
#' the legend guides for `ggbartax`
#' @param keywidth numeric, the keywidth of `guide_legend`.
#' @param keyheight numeric, the keyheight of `guide_legend`.
#' @param ncol integer, the ncol of legend.
#' @param ..., additional parameter.
#' @author ShuangbinXu
#' @importFrom ggplot2 guides
#' @export

taxbarguildes <- function(keywidth=0.4, keyheight=0.4, ncol=5, ...){
	guides(fill=guide_legend(keywidth = keywidth, 
							 keyheight = keyheight,
							 ncol=ncol),...)
}


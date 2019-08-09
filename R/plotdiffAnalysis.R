#' @title plot the clade tree results of different analysis
#' @param taxda data.frame, data.frame contained hierarchical relationship or other classes,
#' such like the tax_data of phyloseq.
#' @param obj object, diffAnalysisClass, the results of diffAnalysis \code{\link[MicrobiotaProcess]{diffAnalysis}}.
#' @param nodedf data.frame, contained the tax and the factor information and(or pvalue).
#' @param factorName character, the names of factor in nodedf.
#' @param removeUnkown logical, whether do not show unkown taxonomy, default is TRUE.
#' @param layout character, the layout of ggtree, but only "rectangular" and "circular" in here, 
#' default is circular.
#' @param size numeric, the size of segment of ggtree, default is 0.6.
#' @param skpointsize numeric, the point size of skeleton of tree, default is 0.8 .
#' @param alpha numeric, the alpha of clade, default is 0.4.
#' @param taxlevel positive integer, the full text of clade, default is 5.
#' @param cladetext numeric, the size of text of clade, default is 2.
#' @param factorLevels list, the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this.
#' @param setColors logical, whether set the color of clade, default is TRUE, 
#' or set FALSE,then use 'scale_fill_manual' setting.
#' @param settheme logical, whether set the theme of tree, default is TRUE,
#' or set FALSE, then use 'theme' setting.
#' @param setlegend logical, whether set the legend, default is TRUE,
#' or set FALSE, then use 'guides' or 'theme' setting.
#' @author Shuangbin Xu
#' @export
ggdiffclade <- function(obj,...){
	UseMethod("ggdiffclade")
}

#' @method ggdiffclade data.frame
#' @rdname ggdiffclade
#' @importFrom ggplot2 geom_point
#' @importFrom magrittr %<>%
#' @export
ggdiffclade.data.frame <- function(taxda, 
					 nodedf,
					 factorName,
					 layout="circular", 
					 size=0.6, 
					 skpointsize=0.8,
					 alpha=0.4,
					 taxlevel=6,
					 cladetext=2,
					 factorLevels=NULL,
					 setColors=TRUE,
					 settheme=TRUE,
					 setlegend=TRUE,
					   ...){
	treedata <- convert_to_treedata(taxda)
	layout %<>% match.arg(c("rectangular", "circular"))
	if (!is.null(factorLevels)){
		nodedf <- setfactorlevels(nodedf, factorLevels)
	}
	nodedf <- getnode(treedata, nodedf)
	p <- treeskeleton(treedata,
				layout=layout,
				size=size,
				pointsize=skpointsize)
	cladecoord <- getcladedf(p, nodedf$node)
	cladecoord <- merge(cladecoord, nodedf, by.x="node", by.y="node")
	labelannotcoord <- getlabeldf(p, nodedf$node)
	if (layout=="rectangular"){taxlevel=7}
	labelannotcoord <- getannotlabel(labelannotcoord, 
									 classlevel=taxlevel)
	labelcoord <- labelannotcoord$labeldf
	annotcoord <- labelannotcoord$annotdf
	if ("pvalue" %in% colnames(nodedf)){
		nodedf[["-log10(pvalue)"]] <- -log10(nodedf$pvalue)
		pointmapping <- aes_(x=~x,y=~y,
							 fill=as.formula(paste("~",factorName)), 
							 size=~-log10(pvalue))
	}else{
		pointmapping <- aes_(x=~x,y=~y,
							 fill=as.formula(paste("~",factorName)))
	}
	p <- p + geom_rect(data=cladecoord,
					   aes_(fill=as.formula(paste("~",factorName)),
						   xmin=~xmin,
						   ymin=~ymin,
						   xmax=~xmax,
						   ymax=~ymax),
					   alpha=alpha,
					   show.legend=FALSE) +
		geom_text(data=labelcoord,
				  aes(x=x,y=y,
					  label=label,
					  angle=angle),
				  size=2) +
		geom_point(data=cladecoord,
				   pointmapping,
				   shape=21)+
		geom_point(data=annotcoord, 
				   aes(x=0,y=0, color=label),
				   size=0, stroke=0) +
		scale_size_continuous(range = c(1, 3))
	if (setColors){
		tmpn <- length(unique(as.vector(nodedf[[match(factorName,colnames(nodedf))]])))
		p <- p + scale_fill_manual(values=getCols(tmpn))
	}
	if (settheme){
		p <- p + ggdiffcladetheme()
	}
	if (setlegend){
		p <- p + ggdiffcladeguides()
	}
	return(p)
}


#' @method ggdiffclade diffAnalysisClass
#' @rdname ggdiffclade
#' @importFrom magrittr %>%
#' @importFrom dplyr filter select
#' @export
ggdiffclade.diffAnalysisClass <- function(obj, removeUnkown=TRUE, ...){
	taxda <- obj@taxda
	#secondvars <- getsecondTRUEvar(obj)
	#nodedfres <- merge(obj@kwres, secondvars, by.x="f", by.y="f") %>% 
	#			select(-c("gfc", "Freq"))
	#nodedfres <- nodedfres[as.vector(nodedfres$f)%in%as.vector(obj@mlres$f),]
	#nodedfres <- merge(nodedfres, obj@mlres$f, by.x="f", by.y="f")
	nodedfres <- tidydiffAnalysis(obj)
	if (removeUnkown){
		tmpflag <- grep("__un_",as.vector(nodedfres$f))
		if (length(tmpflag)>0){
			nodedfres <- nodedfres[-tmpflag,,drop=FALSE]
		}
	}
	p <- ggdiffclade.data.frame(taxda=taxda,
								nodedf=nodedfres,
								factorName=obj@classname,
								...)
	return(p)
}


#' @title significantly discriminative feature barplot
#' @param obj object, diffAnalysisClass see also \code{\link[MicrobiotaProcess]{diffAnalysis}} 
#' or featureMeanMedian see also \code{\link[MicrobiotaProcess]{getMeanMedian}}.
#' @param filepath character, default is NULL, meaning current path. 
#' @param output character, the output dir name, default is "biomarker_barplot".
#' @param xtextsize numeric, the size of axis x label, default is 3.
#' @param removeUnkown logical, whether do not show unkown taxonomy, default is TRUE.
#' @param figwidth numeric, the width of figures, default is 6.
#' @param figheight numeric, the height of figures, default is 3.
#' @param featurename character, the feature name, contained at the objet.
#' @param class character, factor name.
#' @param subclass character, factor name. 
#' @param setColors logical, whether set the colors, default is TRUE, or FALSE,then
#' use scale_color_manual setting.
#' @param settheme logical, whether set the theme, default is TRUE, or FALSE, then
#' use theme and guides setting.
#' @author Shuangbin Xu
#' @export
ggdifftaxbar <- function(obj,...){
	UseMethod("ggdifftaxbar")
}

#' @keywords internal
setGeneric("ggdifftaxbar")

#' @aliases ggdifftaxbar,diffAnalysisClass
#' @rdname ggdifftaxbar
#' @export
setMethod("ggdifftaxbar","diffAnalysisClass",function(obj,
										   filepath=NULL,
										   output="biomarker_barplot",
										   removeUnkown=TRUE,
										   figwidth=6,
										   figheight=3,
										   ...){
	featureda <- obj@originalD
	if (!is.null(obj@normalization)){
		featureda <- featureda / obj@normalization
	}
	sampleda <- obj@sampleda
	featureda <- merge(featureda, sampleda, by=0) %>%
			column_to_rownames(var="Row.names")
	nodedfres <- tidydiffAnalysis(obj)
	featurelist <- as.vector(nodedfres$f)
	if (removeUnkown && length(grep("__un_", featurelist))>0){
			featurelist <- featurelist[-grep("__un_", featurelist)]
	}
	if (ncol(sampleda)>1){
		subclass <- colnames(sampleda)[-match(obj@classname, colnames(sampleda))]
	}else{
		subclass <- obj@classname
	}
	if(is.null(filepath)){filepath <- getwd()}
	filepath <- paste(filepath, output, sep="/")
	dir.create(filepath, showWarnings = FALSE)
	for (vars in featurelist){
		resdf <- getMeanMedian(datameta=featureda, 
							   feature=vars, 
							   subclass=subclass)
		p <- ggdifftaxbar.featureMeanMedian(resdf,
									 vars,
									 obj@classname,
									 subclass,...)
		filename <- paste(filepath, paste0(vars,".svg"), sep="/")
		ggsave(filename, p, device="svg", width = figwidth, height=figheight)
	}
})

#' @method ggdifftaxbar featureMeanMedian
#' @rdname ggdifftaxbar
#' @importFrom ggplot2 aes geom_errorbar 
#' @export
ggdifftaxbar.featureMeanMedian <- function(feMeanMedian, 
									featurename, 
									class, 
									subclass,
									xtextsize=3,
									factorLevels=NULL,
									setColors=TRUE, 
									coloslist=NULL,
									settheme=TRUE,
									...){
	data <- feMeanMedian$singlefedf
	dastatistic <- feMeanMedian$singlefestat
	if (!is.null(factorLevels)){
		data <- setfactorlevels(data, factorLevels)
		dastatistic <- setfactorlevels(dastatistic, factorLevels)
	}
	if (missing(subclass)){subclass <- class}
	p <- ggplot(data, aes_(x=~sample, 
						   y=~RelativeAbundance,
						   fill=as.formula(paste0("~",subclass))))+
		geom_bar(stat="identity") +
		geom_errorbar(data=dastatistic,
					  aes_(x=~sample,
						  ymax=~value,
						  ymin=~value,
						  linetype=~statistic),
					  size=0.5, width=1, 
					  inherit.aes=FALSE)+
		scale_linetype_manual(values=c("solid", "dotted"))+
		facet_grid(as.formula(paste0("~",class)),
				   space="free_x", scales="free_x") + 
		labs(title=featurename) + xlab(NULL)+
		scale_y_continuous(expand=c(0,0),
						   limits=c(0,max(data$RelativeAbundance)*1.05))
	if (settheme){
		p <- p + 
				theme_bw() + 
				guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5, order=1),
					   linetype=guide_legend(keywidth = 0.7, keyheight = 0.5, order=2))+
				theme(plot.title = element_text(face="bold",lineheight=25,hjust=0.5),
					  panel.grid=element_blank(),
					  legend.text = element_text(size=6.5),
					  legend.title=element_text(size=7),
					  legend.background=element_rect(fill=NA),
					  axis.text.x=element_text(angle=-45, hjust = 0, size=xtextsize),
					  panel.spacing = unit(0.2, "mm"), 
					  strip.background = element_rect(colour=NA,fill="grey"))
	}
	if (setColors){
		if (is.null(coloslist)){
			tmpn <- length(unique(as.vector(data[[match(subclass,colnames(data))]])))
			p <- p + scale_fill_manual(values=getCols(tmpn))
		}else{
			p <- p + scale_fill_manual(values=coloslist)
		}
	}
	return(p)
}


#' @title get the mean and median of specific feature.
#' @param datameta data.frame, nrow sample * ncol feature + factor.
#' @param feature character vector, the feature contained in datameta.
#' @param subclass character, factor name.
#' @importFrom dplyr rename group_by_ mutate
#' @importFrom magrittr %>%
#' @author Shuangbin Xu
#' @export
getMeanMedian <- function(datameta, feature, subclass){
	#featureda <- merge(data,sampleda,by=0) %>% rename(sample=Row.names)
	factornames <- colnames(datameta)[!unlist(sapply(datameta,is.numeric))]
	featuredatmp <- datameta %>% rownames_to_column(var="sample") %>%
			select(c("sample", feature, factornames)) %>%
			rename(RelativeAbundance=feature) %>% 
			mutate(sample = factor(sample, levels=sample[order(eval(parse(text=subclass)), -RelativeAbundance)]))
	meantmp <- featuredatmp %>% group_by_(subclass) %>% 
			mutate(value = mean(RelativeAbundance)) %>% mutate(statistic="mean")
	mediantmp <- featuredatmp %>% group_by_(subclass) %>% 
			mutate(value = median(RelativeAbundance)) %>% mutate(statistic="median")
	festatic <- rbind(meantmp, mediantmp) %>% data.frame(check.names=FALSE)
	res <- structure(list(singlefedf=featuredatmp, singlefestat=festatic), 
					 class="featureMeanMedian")
	#attr(res, "class") <- "featureMeanMedian" 
	return (res)
}

#' @title visualization of effect size by the Linear Discriminant Analysis or randomForest
#' @param obj object, diffAnalysisClass see \code{\link[MicrobiotaProcess]{diffAnalysis}}
#' @param removeUnkown logical, whether do not show unkown taxonomy, default is TRUE.
#' @param data data.frame, data.frame contained effect size and the group information. 
#' @param factorName character, the column name contained group information in data.frame.
#' @param effectsizename character, the column name contained effect size information.
#' @param setColors logical, whether set the colors, default is TRUE, or FALSE,then 
#' use scale_color_manual setting.
#' @param settheme logical, whether set the theme, default is TRUE, or FALSE, then
#' use theme and guides setting.
#' @author Shuangbin Xu
#' @export
ggeffectsize <- function(obj,...){
	UseMethod("ggeffectsize")
}

#' @method ggeffectsize data.frame
#' @rdname ggeffectsize
#' @export 
ggeffectsize.data.frame <- function(data, 
									factorName, 
									effectsizename,
									factorLevels=NULL,
									setColors=TRUE,
									settheme=TRUE,
									...){
	if (effectsizename %in% "LDA"){
		xlabtext <- bquote(paste(Log[10],"(",.("LDA"), ")"))
	}else{
		xlabtext <- effectsizename
	}
	if (!is.null(factorLevels)){
		data <- setfactorlevels(data,factorLevels)
	}
	p <- ggplot(data=data, 
				aes_(x=as.formula(paste0("~",effectsizename)),
					 y=~f)) + 
		geom_segment(aes(xend=0, yend=f), 
					 color="grey") + 
		geom_point(aes_(color=as.formula(paste0("~",factorName)))) +
		facet_grid(as.formula(paste0(factorName," ~.")),
				   scales = "free_y", space = "free_y")+
		scale_x_continuous(expand=c(0,0), 
						   limits=c(0, max(data[[effectsizename]])*1.1))+
		ylab(NULL) +
		xlab(xlabtext) 
	if (setColors){
		tmpn <- length(unique(as.vector(data[[factorName]])))
		p <- p + scale_color_manual(values=getCols(tmpn))
	}
	if (settheme){
		p <- p + theme_bw()+
				theme(axis.text.y = element_text(size=7),
					   panel.grid=element_blank(),
					   panel.spacing = unit(0.2, "mm"),
					   legend.title=element_text(size=7),
					   legend.text=element_text(size=6),
					   strip.background = element_rect(colour=NA, 
													   fill="grey"))
	}
	return(p)
}

#' @method ggeffectsize diffAnalysisClass
#' @rdname ggeffectsize
#' @export
ggeffectsize.diffAnalysisClass <- function(obj, removeUnkown=TRUE,...){
	efres <- tidyEffectSize(obj)
	if (removeUnkown && length(grep("__un_",efres$f))){
		efres <- efres[-grep("__un_",efres$f),,drop=FALSE]
	}
	if ("LDA" %in% colnames(obj@mlres)){
		effectsizename <- "LDA"
	}else{
		effectsizename <- "MeanDecreaseAccuracy" 
	}
	p <- ggeffectsize.data.frame(data=efres, 
								 factorName=obj@classname, 
								 effectsizename=effectsizename,
								 ...)
	return (p)
}


#' @importFrom dplyr rename
#' @keywords internal
getnode <- function(treedata, nodedf){
	if (is.null(treedata@data)){stop("The data slot of treedata should't be NULL.")}
	nodelist <- treedata@data[match(as.vector(nodedf[,1]), as.vector(treedata@data$labelnames)),]$node
	if (!"node" %in% colnames(nodedf)){
		nodedf$node <- nodelist
	}else{
		nodedf %>% rename(nodetmp=node)
		nodedf$node <- nodelist
	}
	return(nodedf)
}


#' @importFrom ggtree ggtree
#' @keywords internal
treeskeleton <- function(treedata, layout, size, pointsize=1){
	 p <- ggtree(treedata,
				layout=layout,
				size=size) +
		 geom_point(size=pointsize,
					shape=21, 
					fill="white",
					color="black")
	return(p)
}

#' @importFrom ggplot2 guides
#' @keywords internal
ggdiffcladeguides <- function(...){
	guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5, order=1),
		   size=guide_legend(keywidth = 0.5, keyheight = 0.5, order=2),
		   color = guide_legend(keywidth = 0.1, ncol=1, keyheight = 0.6, order = 3),...)
}

#' @importFrom ggplot2 theme
#' @keywords internal
ggdiffcladetheme <- function(...){
	theme(legend.position="right",
		  legend.margin=margin(0,0,0,0),
		  legend.spacing.y = unit(0.02, "cm"),
		  legend.box.spacing=unit(0.02,"cm"),
		  legend.text = element_text(size=5),
		  legend.title=element_text(size=6),
		  legend.background=element_rect(fill=NA))
}

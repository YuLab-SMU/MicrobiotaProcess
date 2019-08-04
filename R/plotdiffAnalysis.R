#' @title plot the clade tree results of different analysis
#' @param taxda data.frame, data.frame contained hierarchical relationship or other classes,
#' such like the tax_data of phyloseq.
#' @param obj object, diffAnalysisClass, the results of diffAnalysis \code{\link[MicrobiotaProcess]{diffAnalysis}}.
#' @param nodedf data.frame, contained the tax and the factor information and(or pvalue).
#' @param factorName character, the names of factor in nodedf.
#' @param layout character, the layout of ggtree, but only "rectangular" and "circular" in here, 
#' default is circular.
#' @param size numeric, the size of segment of ggtree, default is 0.6.
#' @param alpha numeric, the alpha of clade, default is 0.4.
#' @param taxlevel positive integer, the full text of clade, default is 5.
#' @param cladetext numeric, the size of text of clade, default is 2.
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
					 alpha=0.4,
					 taxlevel=6,
					 cladetext=2,
					 setColors=TRUE,
					 settheme=TRUE,
					 setlegend=TRUE,
					   ...){
	treedata <- as.treedata(taxda)
	layout %<>% match.arg(c("rectangular", "circular"))
	nodedf <- getnode(treedata, nodedf)
	p <- treeskeleton(treedata,
				layout=layout,
				size=size)
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
		scale_size_continuous(range = c(1, 2.5))
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
ggdiffclade.diffAnalysisClass <- function(obj,...){
	taxda <- obj@taxda
	secondvars <- obj@secondvars
	secondvars <- do.call("rbind", 
						  c(secondvars, make.row.names=FALSE))
	secondvars <- secondvars %>% filter(gfc=="TRUE")
	nodedfres <- merge(kwres, secondvars, by.x="f", by.y="f") %>% 
				select(-c("gfc", "Freq"))
	nodedfres <- nodedfres[as.vector(nodedfres$f)%in%as.vector(obj@mlres$f),]
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

#' @method ggdifftaxbar diffAnalysisClass
#' @rdname ggdifftaxbar
#' @export
ggdifftaxbar.diffAnalysisClass <- function(obj,
										   filepath=NULL,
										   output="biomarker_barplot",
										   ...){
	featureda <- obj@originalD
	if (!is.null(test@normalization)){
		featureda <- featureda / test@normalization
	}
	sampleda <- test@sampleda
	featureda <- merge(featureda, sampleda, by=0) %>%
			column_to_rownames(var="Row.names")
	featurelist <- as.vector(obj@mlres$f)
	if (ncol(obj@sampleda)>1){
		subclass <- colnames(sampleda)[-match(obj@classname, colnames(sampleda))]
	}else{
		subclass <- obj@classname
	}
	if(is.null(filepath)){filepath <- getwd()}
	filepath <- paste(filepath,output,sep="/")
	dir.create(filepath, showWarnings = FALSE)
	for (vars in featurelist){
		resdf <- getMeanMedian(datameta=featureda, 
							   feature=vars, 
							   subclass=subclass)
		p <- ggdifftaxbar.featureMeanMedian(resdf,
									 vars,
									 obj@classname,
									 subclass,...)
		filename <- paste(filepath, paste0(i,".svg"), sep="/")
		ggsave(filename, p, device="svg", width = figwidth, height=figheight)
	}
}

#' @method ggdifftaxbar featureMeanMedian
#' @rdname ggdifftaxbar
#' @export
ggdifftaxbar.featureMeanMedian <- function(feMeanMedian, 
									featurename, 
									class, 
									subclass,
									setColors=TRUE, 
									settheme=TRUE,
									...){
	data <- feMeanMedian$singlefedf
	dastatistic <- feMeanMedian$singlefestat
	if (missing(subclass)){subclass <- class}
	p <- ggplot(data, aes_(x=~sample, 
						   y=~RelativeAbundance,
						   fill=as.formula(paste0("~",subclass))))+
		geom_bar(stat="identity") +
		geom_errorbar(data=dastatistic,
					  aes(x=sample,
						  ymax=value,
						  ymin=value,
						  linetype=statistic),
					  size=0.3, width=0.8)+
		facet_wrap(as.formula(paste0("~",class)),
				   nrow=1, scale="free_x") + 
		labs(title=featurename) +
		scale_y_continuous(expand=c(0,0),
						   limits=c(0,max(data$RelativeAbundance)*1.05))
	if (settheme){
		p <- p + 
				theme_bw() + 
				guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5, order=1),
					   linetype=guide_legend(keywidth = 0.5, keyheight = 0.5, order=2))+
				theme(plot.title = element_text(face="bold",lineheight=25,hjust=0.5),
					  panel.grid=element_blank(),
					  legend.text = element_text(size=6.5),
					  legend.title=element_text(size=7),
					  legend.background=element_rect(fill=NA),
					  axis.text.x=element_text(angle=-45, hjust = 0, size=4),
					  panel.spacing = unit(0.2, "mm"), 
					  strip.background = element_rect(colour=NA,fill="grey"))
	}
	if (setColors){
		tmpn <- length(unique(as.vector(data[[match(subclass,colnames(data))]])))
		p <- p + scale_fill_manual(values=getCols(tmpn))
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
			rename(RelativeAbundance=i) %>% 
			mutate(sample = factor(sample, levels=sample[order(eval(parse(text=subclass)), -RelativeAbundance)]))
	meantmp <- featuredatmp %>% group_by_(subclass) %>% 
			mutate(value = mean(RelativeAbundance)) %>% mutate(statistic="mean")
	mediantmp <- featuredatmp %>% group_by_(subclass) %>% 
			mutate(value = median(RelativeAbundance)) %>% mutate(statistic="median")
	festatic <- rbind(meantmp, mediantmp)
	res <- list(singlefedf=featuredatmp, singlefestat=festatic)
	attr(res, "class") <- "featureMeanMedian" 
	return (res)
}

#' @title visualization of effect size by the Linear Discriminant Analysis or randomForest
#' @param obj object, diffAnalysisClass see \code{\link[MicrobiotaProcess]{diffAnalysis}}
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
									setColors=TRUE,
									settheme=TRUE,
									...){
	if (effectsizename %in% "LDA"){
		xlabtext <- bquote(paste(Log[10],"(",.("LDA"), ")"))
	}else{
		xlabtext <- effectsizename
	}
	p <- ggplot(data=efres, 
				aes_(x=as.formula(paste0("~",effectsizename)),
					 y=~f)) + 
		geom_segment(aes(xend=0, yend=f), 
					 color="grey") + 
		geom_point(aes_(color=as.formula(paste0("~",factorName)))) +
		facet_grid(as.formula(paste0(factorName," ~.")),
				   scales = "free_y", space = "free_y")+
		scale_x_continuous(expand=c(0,0), 
						   limits=c(0, max(efres[[effectsizename]])*1.1))+
		ylab(NULL) +
		xlab(xlabtext) 
	if (setColors){
		tmpn <- length(unique(as.vector(efres[[factorName]])))
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
ggeffectsize.diffAnalysisClass <- function(obj, ...){
	secondvars <- do.call("rbind",c(obj@secondvars, 
									make.row.names=FALSE)) 
	secondvars <- secondvars %>% filter(gfc=="TRUE")
	efres <- merge(obj@mlres, secondvars, by.x="f", by.y="f") %>%
			select (-c("gfc", "Freq"))
	if ("LDA" %in% colnames(efres)){
		efres <- efres %>% mutate(f = factor(f, levels=f[order(eval(parse(text=obj@classname)), LDA)]))
		effectsizename <- "LDA"
	}else{
		efres <- efres %>% mutate(f = factor(f, levels=f[order(eval(parse(text=obj@classname)), 
															   MeanDecreaseAccuracy)]))
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
	nodelist <- treedata@data[match(as.vector(nodedf[,1]), as.vector(taxtree@data$labelnames)),]$node
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
treeskeleton <- function(treedata, layout, size){
	 p <- ggtree(treedata,
				layout=layout,
				size=size) +
		 geom_point(size=1,
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
		   color = guide_legend(keywidth = 0.1, keyheight = 0.6, order = 3),...)
}

#' @importFrom ggplot2 theme
#' @keywords internal
ggdiffcladetheme <- function(...){
	theme(legend.position=c(0.9, 0.5),
		  legend.text = element_text(size=6.5),
		  legend.title=element_text(size=7),
		  legend.background=element_rect(fill=NA))
}

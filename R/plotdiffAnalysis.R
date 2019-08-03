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
		tmpn <- length(unique(as.vector(nodedf[[factorName]])))
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
	p <- ggdiffclade.data.frame(taxda=taxda,
								nodedf=nodedfres,
								factorName=obj@classname,
								...)
	return(p)
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

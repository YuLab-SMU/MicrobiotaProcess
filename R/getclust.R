#' @title Hierarchical cluster analysis for the samples
#' @param obj phyloseq, phyloseq class.
#' @param distobj dist, dist class.
#' @param data data.frame, data.frame, default is nrow samples * ncol features.
#' @param taxa_are_rows logical, if the features is in column, it should set FALSE.
#' @param distmethod character, the method of dist, when the obj is data.frame 
#' or phyloseq default is "euclidean". see also \code{\link[MicrobiotaProcess]{getdist}}.
#' @param sampleda data.frame, nrow sample * ncol factor. default is NULL.
#' @param method character, the standardization methods for community ecologists, see also
#' \code{\link[vegan]{decostand}}
#' @param hclustmethod character, the method of hierarchical cluster, default is average.
#' @param tree phylo, the phylo class, see also \code{\link[ape]{as.phylo}}.
#' @param ..., additional parameters.
#' @return clustplotClass object.
#' @author Shuangbin Xu
#' @export
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' subGlobal <- subset_samples(GlobalPatterns, 
#'          SampleType %in% c("Feces", "Mock", "Ocean", "Skin"))
#' hcsample <- getclust(subGlobal, distmethod="jaccard",
#'                   method="hellinger", hclustmethod="average")
getclust <- function(obj,...){
	UseMethod("getclust")
}

#' @method getclust dist
#' @rdname getclust
#' @importFrom ape as.phylo
#' @importFrom stats hclust
#' @export
getclust.dist <- function(distobj,
						  distmethod,
						  sampleda=NULL,
						  hclustmethod="average",
						  ...){
	if (missing(distmethod) && is.null(attr(distobj, "distmethod"))){
		stop("the method of distance should be provided!")
	}
	if (!is.null(attr(distobj, "distmethod"))){
		distmethod <- attr(distobj, "distmethod")
	}
	if (is.null(attr(distobj, "distmethod")) && !missing(distmethod)){
		distmethod <- distmethod
	}
	hclustobj <- hclust(distobj, 
						method=hclustmethod, 
						...)
	phyloobj <- as.phylo(hclustobj)
	clustplot <- new("clustplotClass",
					 hclustphylo=phyloobj,
					 sampleda=sampleda,
					 distmethod=distmethod)
	return(clustplot)
}

#' @method getclust default
#' @rdname getclust
#' @export
getclust.default <- function(data, 
							 distmethod="euclidean",
							 taxa_are_rows=FALSE,
							 sampleda=NULL,
							 tree=NULL,
							 method="hellinger",
							 hclustmethod="average",
							 ...){
	distobj <- getdist(data=data, 
			distmethod=distmethod,
			taxa_are_rows=taxa_are_rows,
			sampleda=sampleda,
			tree=tree,
			type="sample",
			method=method, ...)
	phyloobj <- getclust.dist(distobj, 
							  distmethod=distmethod,
							  sampleda=sampleda,
							  hclustmethod=hclustmethod)
	return(phyloobj)
}

#' @method getclust phyloseq
#' @rdname getclust
#' @export
getclust.phyloseq <- function(obj, 
							  distmethod="euclidean", 
							  method="hellinger",
							  hclustmethod="average",
							  ...){
	distobj <- getdist(obj,
					   distmethod=distmethod,
					   method=method,
					   type="sample",
					   ...)
	sampleda <- checksample(obj)
	phyloobj <- getclust.dist(distobj, 
							  sampleda=sampleda, 
							  hclustmethod=hclustmethod)
	return(phyloobj)
}

#' @title plot the result of hierarchical cluster analysis for the samples
#' @param obj R object, clustplotClass.
#' @param layout character, the layout of tree, see also \code{\link[ggtree]{ggtree}}.
#' @param factorNames character, default is NULL.
#' @param factorLevels list, default is NULL.
#' @param pointsize numeric, the size of point, default is 2.
#' @param fontsize numeric, the size of text of tiplabel, default is 2.6.
#' @param hjust numeric, default is -0.1
#' @param settheme logical, default is TRUE.
#' @param ..., additional params, see also \code{\link[ggtree]{geom_tippoint}}
#' @return the figures of hierarchical cluster.
#' @author Shuangbin Xu
#' @export
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' subGlobal <- subset_samples(GlobalPatterns,
#'          SampleType %in% c("Feces", "Mock", "Ocean", "Skin"))
#' hcsample <- getclust(subGlobal, distmethod="jaccard",
#'                   method="hellinger", hclustmethod="average") 
ggclust <- function(obj,...){
	UseMethod("ggclust")
}

#' @method ggclust clustplotClass
#' @rdname ggclust
#' @importFrom ggtree ggtree %<+% geom_tippoint geom_tiplab geom_tiplab2
#' @importFrom ggplot2 labs element_text
#' @export
ggclust.clustplotClass <- function(obj, 
								   layout="rectangular",
								   factorNames=NULL,
								   factorLevels=NULL,
								   pointsize=2,
								   fontsize=2.6,
								   hjust=-0.1,
								   settheme=TRUE,
								   ...){
	phyloclass <- obj@hclustphylo
	samplehcp <- ggtree(phyloclass, 
						layout=layout)
	if (!is.null(obj@sampleda)){
		sampleda <- data.frame(obj@sampleda, check.names=FALSE)
		phyloclass <- obj@hclustphylo
		sampleda <- sampleda[match(phyloclass$tip.label, rownames(sampleda)),,drop=FALSE]
		sampleda <- data.frame(sample=rownames(sampleda),sampleda, check.names=FALSE)
		rownames(sample) <- NULL
		if(!is.null(factorNames)){
			tmpfactormap <- getfactormap(factorNames)	
		}else{
			tmpfactormap <- getfactormap(colnames(sampleda)[-1])
		}
		if(!is.null(factorLevels)){
			sampleda <- setfactorlevels(sampleda, factorLevels)
		}
		samplehcp <- samplehcp %<+% sampleda + 
				geom_tippoint(tmpfactormap, size=pointsize, ...)
	}
	if (layout=="circular"){
		samplehcp <- samplehcp + geom_tiplab2(size=fontsize, hjust=hjust)
	}else{
		samplehcp <- samplehcp + geom_tiplab(size=fontsize, hjust=hjust)	
	}
	samplehcp <- samplehcp + labs(title=paste0("Hierarchical Cluster of Samples ", "(", obj@distmethod, ")"))
	if (settheme){
		samplehcp <- samplehcp + 
				theme(plot.title = element_text(face="bold", lineheight=25,	hjust=0.5))
	}
	return(samplehcp)
}




#' @title getpca
#' @param obj the phyloseq object
#' @param ... additional parameters
#' @export
getpca <- function(obj,...){
	UseMethod("getpca")
}

#' @title getpca
#' @param data data.frame, nrow sample * ncol feature
#' @param sampleda data.frame, nrow sample * ncol factors.
#' @param method character, the standardization methods for 
#' community ecologists. see \code{\link[vegan]{decostand}}
#' @importFrom vegan decostand
#' @export
getpca.default <- function(data,
				   sampleda=NULL,
				   method="log",
				   ...){
	if (!is.null(method)){
		data <- decostand(data, method=method)
	}
	pca <- prcomp(data, ...)
	#varcontrib <- getvarct(pca) 
	pca <- new("pcasample",
			   pca=pca,
			   #varcontrib=varcontrib,
			   sampleda=sampleda)
	return(pca)
}

#' @title getpca
#' @param obj phyloseq class
#' @param ... additional parameters
#' @export
getpca.phyloseq <- function(obj,...){
	otuda <- checkotu(obj)
	sampleda <- checksample(obj)
	pca <- getpca.default(otuda, sampleda=sampleda,...)
	return(pca)
}

#' @export
ggordpoint <- function(obj, ...){
	UseMethod("ggordpoint")
}

#' @title ggordpoint
#' @param obj prcomp or pcasample class
#' @param pc integer vector, the component index.
#' @param mapping set of aesthetic mapping of ggplot2, default is NULL.
#' @param sampleda data.frame, nrow sample * ncol factors, default is NULL.
#' @param factorNames vector, the names of factors contained sampleda.
#' @param factorLevels list, the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this.
#' @param poinsize numeric, the size of point, default is 2.
#' @param linesize numeric, the line size of segment, default is 0.3.
#' @param arrowsize numeric, the size of arrow, default is 1.5.
#' @param arrowlinecolour character, the color of segment, default is grey.
#' @param biplot logical, whether plot the species, default is TRUE.
#' @param topn integer, the number species have top important contribution, 
#' default is 5.
#' @param settheme logical, whether set the theme for the plot, default is TRUE.
#' @param speciesannot logical, whether plot the species, default is TRUE.
#' @param textsize numeric, the size of text, default is 2.5.
#' @param fontface character, the font face, default is "blod.italic".
#' @param fontfamily character, the font family, default is "sans".
#' @param textlinesize numeric, the segment size in \code{\link[ggrepel]{geom_text_repel}}.
#' @param ... additional parameters, see \code{\link[ggrepel]{geom_text_repel}}.
#' @importFrom ggplot2 ggplot geom_point geom_segment aes_ 
#' @importFrom ggrepel geom_text_repel
#' @export
ggordpoint.default <-  function(obj,
						   pc=c(1,2),
						   mapping=NULL, 
						   sampleda=NULL, 
						   factorNames=NULL,
						   factorLevels=NULL,
						   poinsize=2,
						   linesize=0.3,
						   arrowsize=1.5,
						   arrowlinecolour="grey",
						   biplot=TRUE,
						   topn=5,
						   settheme=TRUE,
						   speciesannot=FALSE,
						   textsize=2.5,
						   fontface="bold.italic",
						   fontfamily="sans",
						   textlinesize=0.02,
						   ...){
	plotcoord <- getcoord(obj)[,pc]
	ev <- obj$sdev^2
	vp <- ev*100/sum(ev)
	xlab_text <- paste("PC",pc[1], "(", round(vp[pc[1]], 2), "%)")
	ylab_text <- paste("PC",pc[2], "(", round(vp[pc[2]], 2), "%)")
	if(is.null(mapping)){
		mapping <- aes_(x=as.formula(paste0("~",colnames(plotcoord)[1])), 
						y=as.formula(paste0("~", colnames(plotcoord)[2])))
	}
	if(!is.null(sampleda)){
		plotcoord <- merge(plotcoord, sampleda, by=0)
		if (!is.null(factorNames)){
			tmpfactormap <- getfactormap(factorNames)
		} else{
			tmpfactormap <- getfactormap(colnames(sampleda))
		}
		mapping <- modifyList(mapping, tmpfactormap)
		if (!is.null(factorLevels)){
			plotcoord <- setfactorlevels(plotcoord, factorLevels)
		}
	}
	p <- ggplot() + 
		 geom_point(data=plotcoord, 
					mapping=mapping, 
					size=poinsize) + 
		 xlab(xlab_text) + 
		 ylab(ylab_text)
	if (biplot){
		varcontrib <- getvarct(obj)
		varcontr <- varcontrib$VarContribution[,pc]
		varlist <- names(sort(rowSums(varcontr), decreasing=T))[1:topn]
		biplotcoord <- varcontrib$VarCoordinates[match(varlist, 
												rownames(varcontrib$VarCoordinates)),
												pc, drop=FALSE]
		biplotmapping <- aes_(x=0, 
							  y=0, 
							  xend=as.formula(paste0("~",colnames(biplotcoord)[1])), 
							  yend=as.formula(paste0("~",colnames(biplotcoord)[2])))
		p <- p + geom_segment(data=biplotcoord, 
							  mapping=biplotmapping,
							  arrow=arrow(length=unit(arrowsize, "mm")),
							  colour=arrowlinecolour,
							  size=linesize)
		if(speciesannot){
			biplotcoord$tax <- rownames(biplotcoord)
			textmapping <- aes_(x=as.formula(paste0("~", colnames(biplotcoord)[1])),
								y=as.formula(paste0("~", colnames(biplotcoord)[2])),
								label=~tax)
			p <- p + geom_text_repel(data=biplotcoord,
									 mapping=textmapping,
									 fontface=fontface,
									 family=fontfamily,
									 size=textsize,
									 segment.size=textlinesize,...)
		}
	}
   	if (settheme){
		p <- p + 
			 geom_vline(xintercept = 0,linetype='dashed',size=0.3,alpha=0.7)+
			 geom_hline(yintercept = 0,linetype='dashed',size=0.3,alpha=0.7)+
			 theme_bw() +
			 theme(panel.grid=element_blank())
	}
	return(p)	
}

#' @export
ggordpoint.pcasample <- function(obj,...){
	pcaobj <- obj@pca
	sampleda <- obj@sampleda
	p <- ggordpoint.default(pcaobj, 
					   sampleda=sampleda,...)
	return(p)
}


#' @export
getcoord <- function(obj){
	UseMethod("getcoord")
}

#' @export
getcoord.prcomp <- function(obj){
	coord <- obj$x
	return(coord)
}

#' @importFrom ggplot2 aes_
#' @keywords internal
getfactormap <- function(namelist){
	if (length(namelist)==1){
		tmpfactormap <- aes_(color=as.formula(paste0("~", namelist[1])),
							 shape=as.formula(paste0("~", namelist[1])))
	}else{
		tmpfactormap <- aes_(color=as.formula(paste0("~", namelist[2])),
							 shape=as.formula(paste0("~", namelist[2])))
	}
	return(tmpfactormap)
}

#' @export
getvarct <- function(obj,...){
	UseMethod("getvarct")
}

#' @export
getvarct.prcomp <- function(obj){
	varcorr <- data.frame(t(apply(obj$rotation,1, varcorf, obj$sdev)),
						  check.names=FALSE)
	varcorr2 <- varcorr^2
	componentcos <- apply(varcorr2, 2, sum)
	#print(componentcos)
	varcontrib <- data.frame(t(apply(varcorr2, 
									 1, 
									 contribution, 
									 componentcos)), 
							 check.names=FALSE)
	res <- list(VarContribution=varcontrib, 
				VarCoordinates=varcorr)
	attr(res, "class") <- "VarContrib"
	return(res)
}

#' @export
getvarct.pcasample <- function(obj){
	pcaobj <- obj@pca
	res <- getvarct(pcaobj)
	return(res)
}

#' @keywords internal
varcorf <- function(varp, sdev){
	varcor <- varp * sdev
	return(varcor)
}

#' @keywords internal
contribution <- function(x, y ){
	contrib <- x*100/y
	return(contrib)
}


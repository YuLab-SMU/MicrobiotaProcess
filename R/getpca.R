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

#' @title
#' @param 
#' @param
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
						   arrowsize=1,
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

#' @title
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


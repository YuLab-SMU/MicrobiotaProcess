#' @title Performs a principal components analysis
#' @param obj phyloseq, phyloseq class or data.frame
#' shape of data.frame is nrow sample * ncol feature.
#' @param sampleda data.frame, nrow sample * ncol factors.
#' @param method character, the standardization methods for
#' community ecologists. see \code{\link[vegan]{decostand}}.
#' @param ... additional parameters, see\code{\link[stats]{prcomp}}.
#' @return pcasample class, contained prcomp class and sample information.
#' @export
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' subGlobal <- subset_samples(GlobalPatterns, 
#'          SampleType %in% c("Feces", "Mock", "Ocean", "Skin"))
#' pcares <- getpca(subGlobal, method="hellinger")
#' #pcaplot <- ggordpoint(pcares, biplot=TRUE, 
#' #                      speciesannot=TRUE,
#' #                      factorNames=c("SampleType"), ellipse=TRUE)
getpca <- function(obj,...){
    UseMethod("getpca")
}

#' @method getpca default
#' @importFrom vegan decostand
#' @importFrom stats prcomp
#' @rdname getpca
#' @export
getpca.default <- function(obj,
    sampleda=NULL,
    method="hellinger",
    ...){
    if (!is.null(method)){
    	obj <- decostand(obj, method=method)
    }
    pca <- prcomp(obj, ...)
    #varcontrib <- getvarct(pca) 
    pca <- new("pcasample",
    		   pca=pca,
    		   #varcontrib=varcontrib,
    		   sampleda=sampleda)
    return(pca)
}

#' @method getpca phyloseq
#' @rdname getpca
#' @export
getpca.phyloseq <- function(obj, method="hellinger", ...){
    otuda <- checkotu(obj)
    sampleda <- checksample(obj)
    pca <- getpca.default(otuda, sampleda=sampleda, method=method, ...)
    return(pca)
}

#' @title ordination plotter based on ggplot2.
#' @param obj prcomp class or pcasample class,
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
#' @param ellipse logical, whether add confidence ellipse to ordinary plot, default is FALSE.
#' @param ellipse_pro numeric, confidence value for the ellipse, default is 0.9.
#' @param ellipse_alpha numeric, the alpha of ellipse, default is 0.2.
#' @param biplot logical, whether plot the species, default is FALSE.
#' @param topn integer or vector, the number species have top important contribution, default is 5.
#' @param settheme logical, whether set the theme for the plot, default is TRUE.
#' @param speciesannot logical, whether plot the species, default is FALSE. 
#' @param textsize numeric, the size of text, default is 2.5.
#' @param fontface character, the font face, default is "blod.italic".
#' @param fontfamily character, the font family, default is "sans".
#' @param textlinesize numeric, the segment size in \code{\link[ggrepel]{geom_text_repel}}.
#' @param ... additional parameters, see \code{\link[ggrepel]{geom_text_repel}}. 
#' @return point figures of PCA or PCoA.
#' @author ShuangbinXu
#' @export
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' subGlobal <- subset_samples(GlobalPatterns,
#'          SampleType %in% c("Feces", "Mock", "Ocean", "Skin"))
#' pcares <- getpca(subGlobal, method="hellinger")
#' pcaplot <- ggordpoint(pcares, biplot=TRUE,
#'                     speciesannot=TRUE,
#'                      factorNames=c("SampleType"), ellipse=TRUE)
ggordpoint <- function(obj, ...){
	UseMethod("ggordpoint")
}

#' @method ggordpoint default
#' @importFrom ggplot2 ggplot geom_point geom_segment aes_ labs arrow unit geom_vline geom_hline theme_bw element_blank element_text  
#' @importFrom ggrepel geom_text_repel
#' @rdname ggordpoint
#' @export
ggordpoint.default <-  function(obj, pc=c(1,2), mapping=NULL, sampleda=NULL, factorNames=NULL, factorLevels=NULL,
    poinsize=2, linesize=0.3, arrowsize=1.5, arrowlinecolour="grey", ellipse=FALSE, ellipse_pro=0.9, 
    ellipse_alpha=0.2, biplot=FALSE, topn=5, settheme=TRUE, speciesannot=FALSE, textsize=2.5,
    fontface="bold.italic", fontfamily="sans", textlinesize=0.02, ...){
    plotcoordclass <- getcoord(obj,pc)
    plotcoord <- plotcoordclass@coord
    xlab_text <- plotcoordclass@xlab	
    ylab_text <- plotcoordclass@ylab
    title_text <- plotcoordclass@title
    if(is.null(mapping)){mapping <- aes_(x=as.formula(paste0("~",colnames(plotcoord)[1])), y=as.formula(paste0("~", colnames(plotcoord)[2])))}
    if(!is.null(sampleda)){
    	plotcoord <- merge(plotcoord, sampleda, by=0)
    	if (!is.null(factorNames)){
    		tmpfactormap <- getfactormap(factorNames)
    		ellipsemapping <- getellipsemap(factorNames)
    	} else{
    		tmpfactormap <- getfactormap(colnames(sampleda))
    		ellipsemapping <- getellipsemap(colnames(sampleda))
    	}
    	mapping <- modifyList(mapping, tmpfactormap)
    	ellipsemapping <- modifyList(mapping, ellipsemapping)
    	if (!is.null(factorLevels)){plotcoord <- setfactorlevels(plotcoord, factorLevels)}
    }
    p <- ggplot() + geom_point(data=plotcoord, mapping=mapping, size=poinsize) + labs(x=xlab_text, y=ylab_text, title=title_text)
    if (ellipse){p <- p + geom_ord_ellipse(data=plotcoord, mapping=ellipsemapping, ellipse_pro=ellipse_pro, alpha=ellipse_alpha, show.legend=FALSE)}
    if (biplot){
    	varcontrib <- getvarct(obj)
    	varcontr <- varcontrib$VarContribution[,pc]
    	tmpvars <- names(sort(rowSums(varcontr), decreasing=TRUE))
    	varlist <- getvarlist(namevector=tmpvars, n=topn)
    	biplotcoord <- varcontrib$VarCoordinates[match(varlist, rownames(varcontrib$VarCoordinates)),pc, drop=FALSE]
    	biplotcoord <- data.frame(biplotcoord, check.names=FALSE)
    	biplotmapping <- aes_(x=0, y=0, xend=as.formula(paste0("~",colnames(biplotcoord)[1])), yend=as.formula(paste0("~",colnames(biplotcoord)[2])))
    	p <- p + geom_segment(data=biplotcoord, mapping=biplotmapping, arrow=arrow(length=unit(arrowsize, "mm")), colour=arrowlinecolour, size=linesize)
    	if(speciesannot){
    		biplotcoord$tax <- rownames(biplotcoord)
    		textmapping <- aes_(x=as.formula(paste0("~", colnames(biplotcoord)[1])), y=as.formula(paste0("~", colnames(biplotcoord)[2])),label=~tax)
    		p<-p+geom_text_repel(data=biplotcoord, mapping=textmapping, fontface=fontface, family=fontfamily, size=textsize,segment.size=textlinesize,...)
    		}
    	}
       	if (settheme){
    		p <- p + geom_vline(xintercept = 0,linetype='dashed',size=0.3,alpha=0.7)+ geom_hline(yintercept = 0,linetype='dashed',size=0.3,alpha=0.7)+
    			 theme_bw() + theme(panel.grid=element_blank(), plot.title = element_text(face="bold",lineheight=25,hjust=0.5))
    	}
    return(p)	
}

#' @keywords internal
getvarlist <- function(namevector, n){
    if (inherits(n, "integer") || inherits(n, "numeric")){
    	if (length(n)==2){
    		varnames <- namevector[n[1]:n[2]]
    	}
    	if (length(n)> 2){
    		varnames <- namevector[n]
    	}
    	if (length(n)==1){
    		varnames <- namevector[seq_len(n)]
    	}
    }
    if (inherits(n, "character")){
    	if (length(n)>1){
    		varnames <- n[match(n, namevector)]
    	}else{
    		stop("the n should be a character vector, integer  or integer vector")
    	}
    }
    return(varnames)
}

#' @method ggordpoint pcasample
#' @rdname ggordpoint
#' @export
ggordpoint.pcasample <- function(obj,...){
    pcaobj <- obj@pca
    sampleda <- obj@sampleda
    p <- ggordpoint.default(pcaobj, 
    				   sampleda=sampleda,...)
    return(p)
}

#' @title get ordination coordinates.
#' @param obj object,prcomp class or pcoa class
#' @param pc integer vector, the component index.
#' @return ordplotClass object.
#' @export
#' @examples
#' require(graphics)
#' data(USArrests)
#' pcares <- prcomp(USArrests, scale = TRUE)
#' coordtab <- getcoord(pcares,pc=c(1, 2))
#' coordtab2 <- getcoord(pcares, pc=c(2, 3))
getcoord <- function(obj, pc){
    UseMethod("getcoord")
}

#' @method getcoord prcomp
#' @rdname getcoord
#' @export
getcoord.prcomp <- function(obj, pc){
    coord <- obj$x[,pc]
    ev <- obj$sdev^2
    vp <- ev*100/sum(ev)
    tmpvp1 <- round(vp[pc[1]],2)
    tmpvp2 <- round(vp[pc[2]],2)
    xlab_text <- paste0("PC",pc[1], "(", tmpvp1, "%)")
    ylab_text <- paste0("PC",pc[2], "(", tmpvp2, "%)")
    title_text <- paste0("PCA - ","PC",pc[1], " VS PC",pc[2])
    ordplotclass <- new("ordplotClass",
    					coord=coord,
    					xlab=xlab_text,
    					ylab=ylab_text,
    					title=title_text)
    return(ordplotclass)
}


#' @importFrom ggplot2 aes_
#' @keywords internal
getfactormap <- function(namelist){
    if (length(namelist)==1){
    	tmpfactormap <- aes_(color=as.formula(paste0("~", namelist[1])))
    }else{
    	tmpfactormap <- aes_(color=as.formula(paste0("~", namelist[1])),
    						 shape=as.formula(paste0("~", namelist[2])))
    }
    return(tmpfactormap)
}

#' @importFrom ggplot2 aes_
#' @importFrom stats as.formula
#' @keywords internal
getellipsemap <- function(namelist){
    tmpellipsemap <- aes_(fill=as.formula(paste0("~", namelist[1])),
    					  group=as.formula(paste0("~", namelist[1])))
    return(tmpellipsemap)
} 


#' @title get the contribution of variables
#' @param obj prcomp class or pcasample class
#' @param ... additional parameters.
#' @return the VarContrib class, contained the 
#' contribution and coordinate of features.
#' @export
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' subGlobal <- subset_samples(GlobalPatterns,
#'          SampleType %in% c("Feces", "Mock", "Ocean", "Skin"))
#' pcares <- getpca(subGlobal, method="hellinger") 
#' varres <- getvarct(pcares)
getvarct <- function(obj,...){
    UseMethod("getvarct")
}

#' @method getvarct prcomp
#' @rdname getvarct
#' @export
getvarct.prcomp <- function(obj,...){
    varcorr <- data.frame(t(apply(obj$rotation,1, varcorf, obj$sdev)),
    					  check.names=FALSE)
    varcorr2 <- varcorr^2
    componentcos <- apply(varcorr2, 2, sum)
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

#' @method getvarct pcasample
#' @rdname getvarct
#' @export 
getvarct.pcasample <- function(obj,...){
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


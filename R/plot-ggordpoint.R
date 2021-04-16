#' @title ordination plotter based on ggplot2.
#' @param obj prcomp class or pcasample class,
#' @param pc integer vector, the component index. 
#' @param mapping set of aesthetic mapping of ggplot2, default is NULL
#' when your want to set it by yourself, only alpha can be setted, and
#' the first element of factorNames has been setted to map 'fill', and
#' the second element of factorNames has been setted to map 'starshape',
#' you can add 'scale_starshape_manual' of 'ggstar' to set the shapes.
#' @param sampleda data.frame, nrow sample * ncol factors, default is NULL. 
#' @param factorNames vector, the names of factors contained sampleda.
#' @param factorLevels list, the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this.
#' @param poinsize numeric, the size of point, default is 2.
#' @param linesize numeric, the line size of segment, default is 0.3. 
#' @param arrowsize numeric, the size of arrow, default is 1.5.
#' @param arrowlinecolour character, the color of segment, default is grey.
#' @param ellipse logical, whether add confidence ellipse to ordinary plot, default is FALSE.
#' @param showsample logical, whether show the labels of sample, default is FALSE.
#' @param labelfactor character, the factor want to be show in label, default is NULL.
#' @param ellipse_pro numeric, confidence value for the ellipse, default is 0.9.
#' @param ellipse_alpha numeric, the alpha of ellipse, default is 0.2.
#' @param biplot logical, whether plot the species, default is FALSE.
#' @param topn integer or vector, the number species have top important contribution, default is 5.
#' @param settheme logical, whether set the theme for the plot, default is TRUE.
#' @param speciesannot logical, whether plot the species, default is FALSE. 
#' @param stroke numeric, the line size of points, default is 0.1.
#' @param fontsize numeric, the size of text, default is 2.5.
#' @param fontface character, the font face, default is "blod.italic".
#' @param fontfamily character, the font family, default is "sans".
#' @param textlinesize numeric, the segment size in \code{\link[ggrepel]{geom_text_repel}}.
#' @param ... additional parameters, see \code{\link[ggrepel]{geom_text_repel}}. 
#' @return point figures of PCA or PCoA.
#' @author Shuangbin Xu
#' @export
#' @examples
#' #don't run in examples
#' #library(phyloseq)
#' #data(GlobalPatterns)
#' #subGlobal <- subset_samples(GlobalPatterns,
#' #         SampleType %in% c("Feces", "Mock", "Ocean", "Skin"))
#' #pcares <- get_pca(subGlobal, method="hellinger")
#' #pcaplot <- ggordpoint(pcares, biplot=TRUE,
#' #                    speciesannot=TRUE,
#' #                     factorNames=c("SampleType"), ellipse=TRUE)
ggordpoint <- function(obj, ...){
	UseMethod("ggordpoint")
}

#' @method ggordpoint default
#' @importFrom ggplot2 ggplot geom_segment aes_string labs arrow unit geom_vline geom_hline theme_bw element_blank element_text  
#' @importFrom ggstar geom_star scale_starshape_manual 
#' @importFrom ggrepel geom_text_repel
#' @rdname ggordpoint
#' @export
ggordpoint.default <-  function(obj, pc=c(1,2), mapping=NULL, sampleda=NULL, factorNames=NULL, factorLevels=NULL,
    poinsize=2, linesize=0.3, arrowsize=1.5, arrowlinecolour="grey", ellipse=FALSE, showsample=FALSE, ellipse_pro=0.9, 
    ellipse_alpha=0.2, biplot=FALSE, topn=5, settheme=TRUE, speciesannot=FALSE, fontsize=2.5, labelfactor=NULL, stroke=0.1,
    fontface="bold.italic", fontfamily="sans", textlinesize=0.02, ...){
    plotcoordclass <- get_coord(obj, pc)
    plotcoord <- plotcoordclass@coord
    xlab_text <- plotcoordclass@xlab	
    ylab_text <- plotcoordclass@ylab
    title_text <- plotcoordclass@title
    defaultmapping <- aes_string(x=colnames(plotcoord)[1], y=colnames(plotcoord)[2])
    if(is.null(mapping)){
        mapping <- defaultmapping
    }else{
        mapping <- modifyList(defaultmapping, mapping)
    }
    if(!is.null(sampleda)){
        plotcoord <- merge(plotcoord, sampleda, by=0)
        if (!is.null(factorNames)){
            tmpfactormap <- get_factormap(factorNames)
            ellipsemapping <- get_ellipsemap(factorNames)
        }else{
            tmpfactormap <- get_factormap(colnames(sampleda))
            ellipsemapping <- get_ellipsemap(colnames(sampleda))
        }
        mapping <- modifyList(mapping, tmpfactormap)
        ellipsemapping <- modifyList(mapping, ellipsemapping)
        ellipsemapping <- modifyList(ellipsemapping,aes_string(starshape=NULL,size=NULL))
        if (!is.null(factorLevels)){plotcoord <- setfactorlevels(plotcoord, factorLevels)}
    }
    p <- ggplot() + geom_star(data=plotcoord, mapping=mapping, size=poinsize, starstroke=stroke) + labs(x=xlab_text, y=ylab_text, title=title_text)
    if ("starshape" %in% names(mapping)){
        shapes <- c(13, 15, 12, 6, 1, 2, 9, 29, 27, 5, 14, 22, 11, 23)[seq_len(length(unique(as.vector(plotcoord[[rlang::as_name(mapping$starshape)]]))))]
        p <- p + scale_starshape_manual(values=shapes)
    }
    if (ellipse){p <- p + geom_ord_ellipse(data=plotcoord,mapping=ellipsemapping,ellipse_pro=ellipse_pro, alpha=ellipse_alpha, fill=NA, show.legend=FALSE, lty=3)}
    if (showsample){
        labelmapping <- modifyList(defaultmapping, aes_string(label="Row.names"))
        if (!is.null(labelfactor)){
            labelmapping <- modifyList(labelmapping, aes_string(label=labelfactor))
        }
        p <- p + geom_text_repel(data=plotcoord, mapping=labelmapping, size=fontsize, segment.size=textlinesize, ...)
    }
    if (biplot){
        varcontrib <- get_varct(obj)
        varcontr <- varcontrib$VarContribution[,pc]
        tmpvars <- names(sort(rowSums(varcontr), decreasing=TRUE))
        varlist <- get_varlist(namevector=tmpvars, n=topn)
        biplotcoord <- varcontrib$VarCoordinates[match(varlist, rownames(varcontrib$VarCoordinates)),pc, drop=FALSE]
        biplotcoord <- data.frame(biplotcoord, check.names=FALSE)
        biplotmapping <- aes_string(x="0", y="0", xend=colnames(biplotcoord)[1], yend=colnames(biplotcoord)[2])
        p <- p + geom_segment(data=biplotcoord, mapping=biplotmapping, arrow=arrow(length=unit(arrowsize, "mm")), colour=arrowlinecolour, size=linesize)
   	    if(speciesannot){
            biplotcoord$tax <- rownames(biplotcoord)
            textmapping <- aes_string(x=colnames(biplotcoord)[1], y=colnames(biplotcoord)[2],label="tax")
            p<-p+geom_text_repel(data=biplotcoord, mapping=textmapping, fontface=fontface, family=fontfamily, size=fontsize,segment.size=textlinesize,...)
        }
   	}
    if (settheme){
        p <- p + geom_vline(xintercept = 0,linetype='dashed',size=0.3,alpha=0.7)+ geom_hline(yintercept = 0,linetype='dashed',size=0.3,alpha=0.7)+
                 theme_bw() + theme(panel.grid=element_blank(), plot.title = element_text(face="bold",lineheight=25,hjust=0.5))
    }
    return(p)
}

#' @keywords internal
get_varlist <- function(namevector, n){
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
#' coordtab <- get_coord(pcares,pc=c(1, 2))
#' coordtab2 <- get_coord(pcares, pc=c(2, 3))
get_coord <- function(obj, pc){
    UseMethod("get_coord")
}

#' @method get_coord prcomp
#' @rdname get_coord
#' @export
get_coord.prcomp <- function(obj, pc){
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


#' @importFrom ggplot2 aes_string
#' @keywords internal
get_factormap <- function(namelist){
    if (length(namelist)==1){
    	tmpfactormap <- aes_string(fill=namelist[1])
    }else{
    	tmpfactormap <- aes_string(fill=namelist[1],
                                   starshape=namelist[2])
    }
    return(tmpfactormap)
}

#' @importFrom ggplot2 aes_string
#' @importFrom stats as.formula
#' @keywords internal
get_ellipsemap <- function(namelist){
    tmpellipsemap <- aes_string(color=namelist[1],
                                group=namelist[1])
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
#' #pcares <- get_pca(subGlobal, method="hellinger") 
#' #varres <- get_varct(pcares)
get_varct <- function(obj,...){
    UseMethod("get_varct")
}

#' @method get_varct prcomp
#' @rdname get_varct
#' @export
get_varct.prcomp <- function(obj,...){
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

#' @method get_varct pcasample
#' @rdname get_varct
#' @export 
get_varct.pcasample <- function(obj,...){
    pcaobj <- obj@pca
    res <- get_varct(pcaobj)
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


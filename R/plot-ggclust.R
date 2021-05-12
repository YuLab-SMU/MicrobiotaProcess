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
#' #don't run in examples
#' #library(phyloseq)
#' #library(ggtree)
#' #library(ggplot2)
#' #data(GlobalPatterns)
#' #subGlobal <- subset_samples(GlobalPatterns,
#' #         SampleType %in% c("Feces", "Mock", "Ocean", "Skin"))
#' #hcsample <- get_clust(subGlobal, distmethod="jaccard",
#' #                  method="hellinger", hclustmethod="average")
#' #hc_p <- ggclust(hcsample, layout = "rectangular",
#' #                pointsize=1, fontsize=0,
#' #                factorNames=c("SampleType")) +
#' #        theme_tree2(legend.position="right",
#' #                    plot.title = element_text(face="bold", lineheight=25,hjust=0.5))
ggclust <- function(obj,...){
    UseMethod("ggclust")
}

#' @method ggclust clustplotClass
#' @rdname ggclust
#' @importFrom ggtree ggtree %<+% geom_tippoint geom_tiplab
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
        if (!"sample" %in% colnames(sampleda)){
            sampleda <- data.frame(sample=rownames(sampleda),sampleda, check.names=FALSE)
        }
        rownames(sample) <- NULL
        if(!is.null(factorNames)){
            tmpfactormap <- set_factormap(factorNames)	
    	}else{
            tmpfactormap <- set_factormap(colnames(sampleda)[-1])
        }
        if(!is.null(factorLevels)){
            sampleda <- setfactorlevels(sampleda, factorLevels)
        }
        samplehcp <- samplehcp %<+% sampleda + 
    			geom_tippoint(tmpfactormap, size=pointsize, ...)
    }
    samplehcp <- samplehcp + 
                 geom_tiplab(size=fontsize, hjust=hjust) + 
                 labs(title=paste0("Hierarchical Cluster of Samples ", "(", obj@distmethod, ")"))
    if (settheme){
        samplehcp <- samplehcp + 
                     theme(plot.title = element_text(face="bold", lineheight=25, hjust=0.5))
    }
    return(samplehcp)
}

#' @importFrom ggplot2 aes_string
#' @keywords internal
set_factormap <- function(namelist){
    if (length(namelist)==1){
        tmpfactormap <- aes_string(color=namelist[1])
    }else{
        tmpfactormap <- aes_string(color=namelist[1],
                                   shape=namelist[2])
    }
    return(tmpfactormap)
}

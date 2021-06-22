#' @title Hierarchical cluster analysis for the samples
#' @param obj phyloseq, phyloseq class or dist class, or 
#' data.frame, data.frame, default is nrow samples * ncol features.
#' @param taxa_are_rows logical, if the features of data.frame(obj) 
#' is in column, it should set FALSE.
#' @param distmethod character, the method of dist, when the 
#' obj is data.frame or phyloseq default is "euclidean". see also 
#' \code{\link[MicrobiotaProcess]{get_dist}}.
#' @param sampleda data.frame, nrow sample * ncol factor. default is NULL.
#' @param method character, the standardization methods for community 
#' ecologists, see also \code{\link[vegan]{decostand}}
#' @param hclustmethod character, the method of hierarchical cluster, 
#' default is average.
#' @param tree phylo, the phylo class, see also \code{\link[ape]{as.phylo}}.
#' @param ..., additional parameters.
#' @return treedata object.
#' @author Shuangbin Xu
#' @export
#' @examples
#' \dontrun{
#' library(phyloseq)
#' data(GlobalPatterns)
#' subGlobal <- subset_samples(GlobalPatterns, 
#'          SampleType %in% c("Feces", "Mock", "Ocean", "Skin"))
#' hcsample <- get_clust(subGlobal, distmethod="jaccard",
#'                   method="hellinger", hclustmethod="average")
#' }
get_clust <- function(obj,...){
    UseMethod("get_clust")
}

#' @method get_clust dist
#' @rdname get_clust
#' @importFrom ape as.phylo
#' @importFrom stats hclust
#' @export
get_clust.dist <- function(obj,
                           distmethod,
                           sampleda=NULL,
                           hclustmethod="average",
                           ...){
    if (missing(distmethod) && is.null(attr(obj, "distmethod"))){
        #stop("the method of distance should be provided!")
        distmethod <- NULL
    }
    if (!is.null(attr(obj, "distmethod"))){
    	distmethod <- attr(obj, "distmethod")
    }
    if (is.null(attr(obj, "distmethod")) && !missing(distmethod)){
    	distmethod <- distmethod
    }
    hclustobj <- hclust(obj, 
                        method=hclustmethod, 
                        ...)
    phyloobj <- as.phylo(hclustobj)
    if (!is.null(sampleda)){
        sampleda %<>% rownames_to_column(var="label")
        clustplot <- phyloobj %>% full_join(sampleda, by="label") 
    }else{
        clustplot <- phyloobj %>% as.treedata
    }
    attr(clustplot, "distmethod") <- distmethod
    return(clustplot)
}

#' @method get_clust data.frame
#' @rdname get_clust
#' @export
get_clust.data.frame <- function(obj, 
                              distmethod="euclidean",
                              taxa_are_rows=FALSE,
                              sampleda=NULL,
                              tree=NULL,
                              method="hellinger",
                              hclustmethod="average",
                              ...){
    distobj <- get_dist(obj, 
                        distmethod=distmethod,
                        taxa_are_rows=taxa_are_rows,
                        sampleda=sampleda,
                        tree=tree,
                        method=method, ...)
    phyloobj <- get_clust.dist(distobj, 
                               distmethod=distmethod,
                               sampleda=sampleda,
                               hclustmethod=hclustmethod)
    return(phyloobj)
}

#' @method get_clust phyloseq
#' @rdname get_clust
#' @export
get_clust.phyloseq <- function(obj, 
                               distmethod="euclidean", 
                               method="hellinger",
                               hclustmethod="average",
                               ...){
    distobj <- get_dist(obj,
                        distmethod=distmethod,
                        method=method,
                        ...)
    sampleda <- checksample(obj)
    phyloobj <- get_clust.dist(distobj, 
                               sampleda=sampleda, 
                               hclustmethod=hclustmethod)
    return(phyloobj)
}

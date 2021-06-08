#' @title NRT (Nearest Relative Index) and NTI (Nearest Taxon Index)
#' 
#' @description
#' calculate NRT and NTI of sample. It is a wrapper method of picante::ses.mpd
#' and picante::ses.mntd
#' @param obj object, data.frame of (nrow sample * ncol taxonomy(feature)) 
#' or phyloseq.
#' @param mindepth numeric, Subsample size for rarefying community.
#' @param sampleda data.frame, sample information, row sample * column factors.
#' @param tree tree object, it can be phylo object or treedata object.
#' @param abundance.weighted logical, whether calculate mean nearest taxon distances for each species 
#' weighted by species abundance, default is TRUE. 
#' @param ... additional arguments see also "ses.mpd" and "ses.mntd" of "picante".
#' @return alphasample object contained NRT and NTI.
#' @rdname get_NRI_NTI-methods
#' @author Shuangbin Xu
#' @export
setGeneric("get_NRI_NTI", function(obj, ...){standardGeneric("get_NRI_NTI")})

#' @aliases get_NRI_NTI,matrix
#' @rdname get_NRI_NTI-methods
#' @export
setMethod("get_NRI_NTI", "matrix", function(obj, mindepth, sampleda, tree, abundance.weighted=TRUE,...){
    if (missing(mindepth) || is.null(mindepth)){
           mindepth <- min(rowSums(obj))
    }
    obj <- rrarefy(obj, mindepth)
    treedist <- cal_treedist(tree=tree)
    resnri <- picante::ses.mpd(samp=obj, dis=treedist, abundance.weighted=abundance.weighted, ...)
    resnti <- picante::ses.mntd(samp=obj, dis=treedist, abundance.weighted=abundance.weighted, ...)
    resda <- data.frame(NRI = -resnri$mpd.obs.z,
                      NTI = -resnti$mntd.obs.z
                      )
    if (missing(sampleda)){
        sampleda <- NULL
    }
    res <- new("alphasample",
               alpha=resda,
               sampleda=sampleda)
    return(res)
})     

#' @aliases get_NRI_NTI,data.frame
#' @rdname get_NRI_NTI-methods
#' @export
setMethod("get_NRI_NTI", "data.frame", function(obj, mindepth, sampleda, tree, abundance.weighted=TRUE, ...){
    obj <- obj[,colSums(obj)>0,drop=FALSE]
    obj <- as.matrix(obj)
    res <- get_NRI_NTI(obj=obj, mindepth=mindepth, sampleda=sampleda, tree=tree, abundance.weighted=abundance.weighted, ...)
    return(res)
})

#' @aliases get_NRI_NTI,phyloseq
#' @rdname get_NRI_NTI-methods
#' @export
setMethod("get_NRI_NTI", "phyloseq", function(obj, mindepth, abundance.weighted=TRUE, ...){
    otuda <- checkotu(obj)
    sampleda <- checksample(obj)
    if (is.null(obj@phy_tree)){
        stop ("The tree should be provided, but the phyloseq does not have phy_tree slot!")
    }
    res <- get_NRI_NTI(obj=otuda, 
                       mindepth=mindepth, 
                       sampleda=sampleda, 
                       tree=obj@phy_tree, 
                       abundance.weighted=abundance.weighted,
                       ...)
    return(res)
})

cal_treedist <- function(tree){
    if (inherits(tree, "phylo")){
        treedist <- ape::cophenetic.phylo(tree)
    }else if (inherits(tree, "treedata")){
        treedist <- ape::cophenetic.phylo(tree@phylo)
    }else{
        stop("the tree should be phylo object or treedata object of tidytree")
    }
    return (treedist)
}

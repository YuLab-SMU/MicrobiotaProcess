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
        clustplot <- phyloobj %>% treeio::full_join(sampleda, by="label") 
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

#' Hierarchical cluster analysis for the samples with MPSE or tbl_mpse object
#' @rdname mp_cal_clust-methods
#' @param .data the MPSE or tbl_mpse object
#' @param .abundance the name of abundance to be calculated.
#' @param distmethod the method of distance.
#' @param hclustmethod the method of hierarchical cluster
#' @param action a character "add" will return a MPSE object with the cluster 
#' result as a attributes, and it can be extracted with 'object %>% mp_extract_cluster()', 
#' "only" or "get" will return 'treedata' object, default is 'get'.
#' @param ... additional parameters
#' @return update object with the action argument, the treedata object 
#' contained hierarchical cluster analysis of sample, it can be visualized 
#' with 'ggtree' directly.
#' @export
#' @author Shuangbin Xu
#' @examples
#' library(ggtree)
#' library(ggplot2)
#' data(mouse.time.mpse)
#' res <- mouse.time.mpse %>%
#'  mp_decostand(.abundance=Abundance) %>% 
#'  mp_cal_clust(.abundance=hellinger, distmethod="bray")
#' res
#' res %>%
#'  ggtree() + 
#'  geom_tippoint(aes(color=time))
setGeneric("mp_cal_clust", function(.data, .abundance, distmethod="bray", hclustmethod="average", action="get", ...)standardGeneric("mp_cal_clust"))

.internal_cal_clust <-  function(.data, .abundance, distmethod="bray", hclustmethod="average", action="get", ...){
    action %<>% match.arg(c("get", "only", "add"))

    .abundance <- rlang::enquo(.abundance)

    if (inherits(.data, "MPSE")){
       flag <- !distmethod %in% colnames(.data@colData)
    }else{
       flag <- !distmethod %in% colnames(.data)
    }

    if (flag){
        if (rlang::quo_is_missing(.abundance)){
            rlang::abort("The .abundance must be required, when the distmethod is not present in the object.")
        }else{
            .data %<>% mp_cal_dist(.abundance=!!.abundance, distmethod=distmethod, action="add")
        }
    }
    
    distobj <- .data %>% mp_extract_dist(distmethod=distmethod)

    res <- distobj %>%
           hclust(method=hclustmethod, ...) %>%
           ape::as.phylo() %>%
           tidytree::left_join(
               y = .data %>% 
                    mp_extract_sample(),
               by=c("label"="Sample")
           )
    if (action %in% c("get", "only")){
        return(res)
    }else if (action=="add"){
        message("The result was added to the internal attributes of the object")
        message("It can be extracted via object %>% mp_extract_internal_attr(name='SampleClust') !")
        .data %<>%
             add_internal_attr(object=res, name="SampleClust")
        return(.data)
    }
}

#' @rdname mp_cal_clust-methods
#' @aliases mp_cal_clust,MPSE
#' @exportMethod mp_cal_clust
setMethod("mp_cal_clust", signature(.data="MPSE"), .internal_cal_clust)

#' @rdname mp_cal_clust-methods
#' @aliases mp_cal_clust,tbl_mpse
#' @exportMethod mp_cal_clust
setMethod("mp_cal_clust", signature(.data="tbl_mpse"), .internal_cal_clust)

#' @rdname mp_cal_clust-methods
#' @aliases mp_cal_clust,grouped_df_mpse
#' @exportMethod mp_cal_clust
setMethod("mp_cal_clust", signature(.data="grouped_df_mpse"), .internal_cal_clust)

#' @title NRI (Nearest Relative Index) and NTI (Nearest Taxon Index)
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
#' @param force logical whether calculate the index even the count of otu is
#' not rarefied, default is FALSE. If it is TRUE, meaning the rarefaction is not be
#' performed automatically.
#' @param seed integer a random seed to make the result reproducible, default is 123.
#' @param ... additional arguments see also "ses.mpd" and "ses.mntd" of "picante".
#' @return alphasample object contained NRT and NTI.
#' @rdname get_NRI_NTI-methods
#' @author Shuangbin Xu
#' @export
setGeneric("get_NRI_NTI", function(obj, ...){standardGeneric("get_NRI_NTI")})

#' @aliases get_NRI_NTI,matrix
#' @rdname get_NRI_NTI-methods
#' @export
setMethod("get_NRI_NTI", "matrix", function(obj, mindepth, sampleda, tree, abundance.weighted=TRUE, force=FALSE, seed=123, ...){
    if (missing(mindepth) || is.null(mindepth)){
           mindepth <- min(rowSums(obj))
    }
    if (obj %>% rowSums %>% var != 0 && !force){
        obj <- vegan::rrarefy(obj, mindepth)
    }
    treedist <- cal_treedist(tree=tree)
    resnri <- withr::with_seed(seed, picante::ses.mpd(samp=obj, dis=treedist, abundance.weighted=abundance.weighted, ...))
    resnti <- withr::with_seed(seed, picante::ses.mntd(samp=obj, dis=treedist, abundance.weighted=abundance.weighted, ...))
    resda <- data.frame(NRI = abs(resnri$mpd.obs.z),
                      NTI = abs(resnti$mntd.obs.z)
                      ) %>% magrittr::set_rownames(rownames(resnri))
    
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

#' Calculating NRI (Nearest Relative Index) and NTI (Nearest Taxon Index) with MPSE or tbl_mpse object
#' @param .data object, MPSE or tbl_mpse object
#' @param .abundance The column name of OTU abundance column to be calculate.
#' @param action character it has three options, "add" joins the new information
#' to the input tbl (default), "only" return a non-redundant tibble with the just
#' new information, ang 'get' return a 'alphasample' object.
#' @param abundance.weighted logical, whether calculate mean nearest taxon distances for each species
#' weighted by species abundance, default is TRUE.
#' @param force logical whether calculate the alpha index even the '.abundance' is
#' not rarefied, default is FALSE.
#' @param seed integer a random seed to make the result reproducible, default is 123.
#' @param ... additional arguments see also "ses.mpd" and "ses.mntd" of "picante".
#' @return update object.
#' @rdname mp_cal_NRI_NTI-methods
#' @author Shuangbin Xu
#' @export
setGeneric("mp_cal_NRI_NTI", function(.data, .abundance, action="add", abundance.weighted=TRUE, force=FALSE, seed=123, ...)standardGeneric("mp_cal_NRI_NTI"))

#' @rdname mp_cal_NRI_NTI-methods
#' @aliases mp_cal_NRI_NTI,MPSE
#' @exportMethod mp_cal_NRI_NTI
setMethod("mp_cal_NRI_NTI", signature(.data="MPSE"), function(.data, .abundance, action="add", abundance.weighted=TRUE, force=FALSE, seed=123, ...){
    action %<>% match.arg(c("add", "only", "get"))

    res <- .internal_cal_NRI_NTI(.data=.data,
                                 .abundance=!!rlang::enquo(.abundance),
                                 abundance.weighted=abundance.weighted,
                                 force=force,
                                 seed = seed,
                                 ...)
    if (action=="get"){
        res@sampleda <- .data %>%
                         mp_extract_sample() %>%
                         tibble::column_to_rownames(var="Sample")
        return(res)
    }

    da <- res@alpha %>% as_tibble(rownames="Sample")

    da <- .data %>%
          mp_extract_sample() %>%
          left_join(da, by="Sample") #%>%
          #column_to_rownames(var="Sample")
    if (action=="add"){
        .data@colData <- da %>% 
                         column_to_rownames(var="Sample") %>% 
                         S4Vectors::DataFrame(check.names=FALSE)
        return(.data)
    }else if (action=="only"){
        return(da)
    }

})
          
.internal_cal_NRI_NTI <- function(.data, .abundance, abundance.weighted=TRUE, force=FALSE, seed=123, ...){

    .abundance <- rlang::enquo(.abundance)
    if (rlang::quo_is_null(.abundance)){
        .abundance <- as.symbol("RareAbundance")
    }

    otutree <- .data %>%
               mp_extract_tree(type="otutree")

    if (is.null(otutree)){
        rlang::abort("The otutree is required to calculate the NRI and NTI !")
    }

    if (!valid_rare(.data, .abundance=.abundance) && !force){
        glue::glue("The rarefied abundance of species might not be provided. Rarefaction of all
                    observations is performed automatically. If you still want to calculate the
                    alpha index with the '.abundance', you can set 'force=TRUE'. ")
        .data <- mp_rrarefy(.data=.data, ...)
        .abundance <- as.symbol("RareAbundance")
    }

    indexda <- .data %>% 
               mp_extract_assays(.abundance=!!.abundance, byRow=FALSE) %>%
               as.matrix() %>%
               get_NRI_NTI(tree=otutree@phylo, force=TRUE, abundance.weighted=abundance.weighted, seed=seed, ...)

    return(indexda)
}

.internal_cal_NRI_NTI_ <- function(.data, .abundance, action="add", abundance.weighted=TRUE, force=FALSE, seed=123, ...){
    action %<>% match.arg(c("add", "only", "get"))
    res <- .internal_cal_NRI_NTI(.data=.data,
                                 .abundance=!!rlang::enquo(.abundance),
                                 force=force,
                                 seed = seed,
                                 ...)
    if (action=="get"){
        res@sampleda <- .data %>%
                         mp_extract_sample() %>%
                         tibble::column_to_rownames(var="Sample")
        return(res)
    }

    da <- res@alpha %>% as_tibble(rownames="Sample")

    da <- .data %>%
          mp_extract_sample() %>%
          left_join(da, by="Sample") #%>%
          #column_to_rownames(var="Sample")
    if (action=="add"){
        samplevar <- .data %>% attr("samplevar")
        assaysvar <- .data %>% attr("assaysvar")
        othernm <- colnames(.data)[!colnames(.data) %in% c("OTU", "Sample", assaysvar, samplevar)]
        .data %<>% left_join(da, by="Sample") %>%
                   select(c("OTU", "Sample", assaysvar, samplevar, colnames(da), othernm))
        return(.data)
    }else if (action=="only"){
        return(da)
    }
}

#' @rdname mp_cal_NRI_NTI-methods
#' @aliases mp_cal_NRI_NTI,tbl_mpse
#' @exportMethod mp_cal_NRI_NTI
setMethod("mp_cal_NRI_NTI", signature(.data="tbl_mpse"), .internal_cal_NRI_NTI_)

#' @rdname mp_cal_NRI_NTI-methods
#' @aliases mp_cal_NRI_NTI,grouped_df_mpse
#' @exportMethod mp_cal_NRI_NTI
setMethod("mp_cal_NRI_NTI", signature(.data="grouped_df_mpse"), .internal_cal_NRI_NTI_)

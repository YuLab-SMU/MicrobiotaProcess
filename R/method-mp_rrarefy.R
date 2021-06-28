##' @docType methods
##' @name mp_rrarefy 
##' @rdname mp_rrarefy-methods
##' @title mp_rrarefy method
##' @param obj phyloseq or tbl_mpse object
##' @param raresize integer Subsample size for rarefying community.
##' @param trimOTUs logical Whether to remove the otus that are no 
##' longer present in any sample after rarefaction
##' @param seed a random seed to make the rrarefy reproducible,
##' default is 123.
##' @return update object
##' @export
setGeneric("mp_rrarefy", function(obj, raresize, trimOTUs=TRUE, seed=123){standardGeneric("mp_rrarefy")}) 

##' @rdname mp_rrarefy-methods
##' @aliases mp_rrarefy,matrix
##' @exportMethod mp_rrarefy
##' @source 
##' mp_rrarefy for matrix object is a wrapper method of vegan::rrarefy from the vegan 
##' package. 
##' @seealso
##' \link[vegan]{rrarefy}
setMethod("mp_rrarefy", signature(obj="matrix"), function(obj, raresize, trimOTUs=TRUE, seed=123){
    if (missing(raresize)||is.null(raresize)){
        raresize <- min(rowSums(obj))
    }
    if (is.na(seed) || is.null(seed)){
        message("To reproduce and it was not be set, the random seed will set to 123 automatically.")
        seed <- 123
    }
    res <- withr::with_seed(seed, vegan::rrarefy(x=obj, sample=raresize))
    if (trimOTUs){
        removeOTU <- colSums(res)==0
        if (sum(removeOTU) > 0){
            message(sum(removeOTU), " OTUs were removed because they are no longer present in any sample after rarefaction!")
        }
        res <- res[, !removeOTU]
    }
    return(res)
})


##' @rdname mp_rrarefy-methods
##' @aliases mp_rrarefy,data.frame
##' @exportMethod mp_rrarefy
setMethod("mp_rrarefy", signature(obj="data.frame"), function(obj, raresize, trimOTUs=TRUE, seed=123){
    obj <- as.matrix(obj)
    res <- mp_rrarefy(obj=obj, raresize=raresize, trimOTUs=trimOTUs, seed=seed)
    res <- data.frame(res) 
    return(res)
})

##' @rdname mp_rrarefy-methods
##' @aliases mp_rrarefy,phyloseq
##' @exportMethod mp_rrarefy
setMethod("mp_rrarefy", signature(obj="phyloseq"), function(obj, raresize, trimOTUs=TRUE, seed=123){
    otuda <- get_otudata(obj)
    res <- mp_rrarefy(obj=otuda, raresize=raresize, trimOTUs=trimOTUs, seed=seed)
    obj@otu_table <- otu_table(res, taxa_are_rows=FALSE)
    if (!is.null(obj@sam_data)){
        obj@sam_data <- obj@sam_data[rownames(obj@sam_data) %in% rownames(res),,drop=FALSE]
    }
    if (!is.null(obj@tax_table)){
        obj@tax_table <- obj@tax_table[rownames(obj@tax_table) %in% colnames(res),,drop=FALSE]
    }
    if (!is.null(obj@phy_tree)){
        obj@phy_tree <- ape::keep.tip(obj@phy_tree, tip=colnames(res))
    }
    if (!is.null(obj@refseq)){
        obj@refseq <- obj@refseq[colnames(res)]
    }
    return (obj)
})

##' @rdname mp_rrarefy-methods
##' @aliases mp_rrarefy,tbl_mpse
##' @exportMethod mp_rrarefy
setMethod("mp_rrarefy", signature(obj="tbl_mpse"), function(obj, raresize, trimOTUs=TRUE, seed=123){
    obj <- obj %>% as.phyloseq()
    res <- mp_rrarefy(obj=obj, raresize=raresize, trimOTUs=trimOTUs, seed=seed)
    res %<>% as_tibble()
    return (res)
})

##' @rdname mp_rrarefy-methods
##' @aliases mp_rrarefy,grouped_df_mpse
##' @exportMethod mp_rrarefy
setMethod("mp_rrarefy", signature(obj="grouped_df_mpse"), function(obj, raresize, trimOTUs=TRUE, seed=123){
    obj %<>% ungroup
    res <- mp_rrarefy(obj=obj, raresize=raresize, trimOTUs=trimOTUs, seed=seed)
    res %<>% as_tibble()
    return (res)
})

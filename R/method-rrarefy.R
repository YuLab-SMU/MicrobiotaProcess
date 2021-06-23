##' @docType methods
##' @name rrarefy 
##' @rdname rrarefy-methods
##' @title rrarefy method
##' @param obj phyloseq or tbl_ps object
##' @param raresize integer Subsample size for rarefying community.
##' @param trimOTUs logical Whether to remove the otus that are no 
##' longer present in any sample after rarefaction
##' @param seed a random seed to make the rrarefy reproducible,
##' default is 123.
##' @return update object
##' @export
setGeneric("rrarefy", function(obj, raresize, trimOTUs=TRUE, seed=123){standardGeneric("rrarefy")}) 

##' @rdname rrarefy-methods
##' @aliases rrarefy,matrix
##' @exportMethod rrarefy
##' @source 
##' rrarefy for matrix object is a wrapper method of vegan::rrarefy from the vegan 
##' package. 
##' @seealso
##' \link[vegan]{rrarefy}
setMethod("rrarefy", signature(obj="matrix"), function(obj, raresize, trimOTUs=TRUE, seed=123){
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


##' @rdname rrarefy-methods
##' @aliases rrarefy,data.frame
##' @exportMethod rrarefy
setMethod("rrarefy", signature(obj="data.frame"), function(obj, raresize, trimOTUs=TRUE, seed=123){
    obj <- as.matrix(obj)
    res <- rrarefy(obj=obj, raresize=raresize, trimOTUs=trimOTUs, seed=seed)
    res <- data.frame(res) 
    return(res)
})

##' @rdname rrarefy-methods
##' @aliases rrarefy,phyloseq
##' @exportMethod rrarefy
setMethod("rrarefy", signature(obj="phyloseq"), function(obj, raresize, trimOTUs=TRUE, seed=123){
    otuda <- get_otudata(obj)
    res <- rrarefy(obj=otuda, raresize=raresize, trimOTUs=trimOTUs, seed=seed)
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

##' @rdname rrarefy-methods
##' @aliases rrarefy,tbl_ps
##' @exportMethod rrarefy
setMethod("rrarefy", signature(obj="tbl_ps"), function(obj, raresize, trimOTUs=TRUE, seed=123){
    obj <- obj %>% as.phyloseq()
    res <- rrarefy(obj=obj, raresize=raresize, trimOTUs=trimOTUs, seed=seed)
    res %<>% as_tibble()
    return (res)
})

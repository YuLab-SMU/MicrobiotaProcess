##' @docType methods
##' @name mp_rrarefy 
##' @rdname mp_rrarefy-methods
##' @title mp_rrarefy method
##' @param .data MPSE or tbl_mpse object
##' @param raresize integer Subsample size for rarefying community.
##' @param trimOTU logical Whether to remove the otus that are no 
##' longer present in any sample after rarefaction
##' @param seed a random seed to make the rrarefy reproducible,
##' default is 123. 
##' @return update object
##' @export
##' @author Shuangbin Xu
##' @examples
##' data(mouse.time.mpse)
##' mouse.time.mpse %>% mp_rrarefy()
setGeneric("mp_rrarefy", function(.data, raresize, trimOTU=FALSE, seed=123){standardGeneric("mp_rrarefy")}) 

##' @rdname mp_rrarefy-methods
##' @aliases mp_rrarefy,matrix
##' @exportMethod mp_rrarefy
##' @source 
##' mp_rrarefy for matrix object is a wrapper method of vegan::rrarefy from the vegan 
##' package. 
##' @seealso
##' \link[vegan]{rrarefy}
setMethod("mp_rrarefy", signature(.data="matrix"), function(.data, raresize, trimOTU=FALSE, seed=123){
    if (missing(raresize)||is.null(raresize)){
        raresize <- min(rowSums(.data))
    }
    if (is.na(seed) || is.null(seed)){
        message("To reproduce and it was not be set, the random seed will set to 123 automatically.")
        seed <- 123
    }
    res <- withr::with_seed(seed, vegan::rrarefy(x=.data, sample=raresize))
    if (trimOTU){
        removeOTU <- colSums(res)==0
        if (sum(removeOTU) > 0){
            message(sum(removeOTU), " OTUs were removed because they are no longer present in any sample after ",
                    "rarefaction, if you want to keep them you can set 'trimOTU = FALSE' !")
        }
        res <- res[, !removeOTU]
    }
    return(res)
})


##' @rdname mp_rrarefy-methods
##' @aliases mp_rrarefy,data.frame
##' @exportMethod mp_rrarefy
setMethod("mp_rrarefy", signature(.data="data.frame"), function(.data, raresize, trimOTU=FALSE, seed=123){
    .data <- as.matrix(.data)
    res <- mp_rrarefy(.data=.data, raresize=raresize, trimOTU=trimOTU, seed=seed)
    res <- data.frame(res) 
    return(res)
})

# ##' @rdname mp_rrarefy-methods
# ##' @aliases mp_rrarefy,phyloseq
# ##' @exportMethod mp_rrarefy
# setMethod("mp_rrarefy", signature(obj="phyloseq"), function(obj, raresize, trimOTU=TRUE, seed=123){
#     otuda <- get_otudata(obj)
#     res <- mp_rrarefy(obj=otuda, raresize=raresize, trimOTU=trimOTU, seed=seed)
#     obj@otu_table <- otu_table(res, taxa_are_rows=FALSE)
#     if (!is.null(obj@sam_data)){
#         obj@sam_data <- obj@sam_data[rownames(obj@sam_data) %in% rownames(res),,drop=FALSE]
#     }
#     if (!is.null(obj@tax_table)){
#         obj@tax_table <- obj@tax_table[rownames(obj@tax_table) %in% colnames(res),,drop=FALSE]
#     }
#     if (!is.null(obj@phy_tree)){
#         obj@phy_tree <- ape::keep.tip(obj@phy_tree, tip=colnames(res))
#     }
#     if (!is.null(obj@refseq)){
#         obj@refseq <- obj@refseq[colnames(res)]
#     }
#     return (obj)
# })


##' @rdname mp_rrarefy-methods
##' @aliases mp_rrarefy,MPSE
##' @exportMethod mp_rrarefy
setMethod("mp_rrarefy", signature(.data="MPSE"), function(.data, raresize, trimOTU=FALSE, seed=123){
    allassays <- SummarizedExperiment::assays(.data) %>% as.list()
    if ("RareAbundance" %in% names(allassays)){
        message("The RareAbundance was in the MPSE object, please check whether it has been rarefied !")
        return(.data)
    }
    rare <- mp_rrarefy(.data=t(allassays[["Abundance"]]), raresize=raresize, trimOTU=FALSE, seed=123) %>% t()
    SummarizedExperiment::assays(.data)@listData <- c(allassays, list(RareAbundance=rare))
    if (trimOTU){
        removeOTU <- rowSums(rare)==0
        if (sum(removeOTU) > 0){
            message(sum(removeOTU), " OTUs were removed because they are no longer present in any sample after ",
                    "rarefaction, if you want to keep them you can set 'trimOTU = FALSE' !")
            .data <- .data[!rownames(.data) %in% rownames(rare[removeOTU,]), ]
        }
    }
    return(.data)
})

##' @rdname mp_rrarefy-methods
##' @aliases mp_rrarefy,tbl_mpse
##' @exportMethod mp_rrarefy
setMethod("mp_rrarefy", signature(.data="tbl_mpse"), function(.data, raresize, trimOTU=FALSE, seed=123){
    assaysvar <- attr(.data, "assaysvar")
    otutree <- attr(.data, "otutree")
    taxatree <- attr(.data, "taxatree")
    if ("RareAbundance" %in% assaysvar){
        message("The RareAbundance was in the MPSE object, please check whether it has been rarefied !")
        return(.data)
    }
    tmpotu <- .data %>% 
              tibble::as_tibble() %>% 
              select(c("OTU", "Sample", "Abundance")) %>%
              tidyr::pivot_wider(names_from="OTU", values_from="Abundance") %>% 
              tibble::column_to_rownames(var="Sample")
    rare <- mp_rrarefy(.data=tmpotu, raresize=raresize, trimOTU=FALSE, seed=seed) %>%
            t() %>% 
            tibble::as_tibble(rownames="OTU") %>% 
            tidyr::pivot_longer(!"OTU", names_to="Sample", values_to="RareAbundance")
    othernms <- colnames(.data)[!colnames(.data) %in% c("OTU", "Sample", assaysvar)]
    res <- .data %>% left_join(rare, by=c("OTU", "Sample")) %>% 
           select(c("OTU", "Sample", assaysvar, "RareAbundance", othernms))
    if (trimOTU){
        removeOTU <- res %>% 
                     dplyr::group_by(.data$OTU) %>% 
                     dplyr::summarise(Total=sum(.data$RareAbundance)) %>% 
                     filter(.data$Total==0) %>% 
                     select("OTU") %>%
                     unlist(use.names=FALSE)
        if (length(removeOTU) > 0){
            message(length(removeOTU), " OTUs were removed because they are no longer present in any sample after ",
                    "rarefaction, if you want to keep them you can set 'trimOTU = FALSE' !")
            res %<>% dplyr::filter(!.data$OTU %in% removeOTU)
        }
    }
    res <- add_attr.tbl_mpse(x1 = res, x2 = .data)
    attr(res, "assaysvar") <- c(assaysvar, "RareAbundance")
    return (res)
})

##' @rdname mp_rrarefy-methods
##' @aliases mp_rrarefy,grouped_df_mpse
##' @exportMethod mp_rrarefy
setMethod("mp_rrarefy", signature(.data="grouped_df_mpse"), function(.data, raresize, trimOTU=FALSE, seed=123){
    tmpgroups <- attr(.data, "groups")
    groupvars <- names(tmpgroups)[names(tmpgroups) != ".rows"]
    groupvars <- lapply(groupvars, function(i) as.symbol(i))
    .data %<>% ungroup
    res <- mp_rrarefy(.data=.data, raresize=raresize, trimOTU=trimOTU, seed=seed)
    res <- do.call(group_by, c(list(res), groupvars))
    return (res)
})

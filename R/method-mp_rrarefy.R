##' @docType methods
##' @name mp_rrarefy 
##' @rdname mp_rrarefy-methods
##' @title mp_rrarefy method
##' @param .data MPSE or tbl_mpse object
##' @param .abundance the name of OTU(feature) abundance column,
##' default is Abundance.
##' @param raresize integer Subsample size for rarefying community.
##' @param trimOTU logical Whether to remove the otus that are no 
##' longer present in any sample after rarefaction
##' @param trimSample logical whether to remove the samples that 
##' do not have enough abundance (raresize), default is FALSE.
##' @param seed a random seed to make the rrarefy reproducible,
##' default is 123. 
##' @param ... additional parameters, meaningless now.
##' @return update object
##' @seealso [mp_extract_assays()] and [mp_decostand()]
##' @export
##' @author Shuangbin Xu
##' @examples
##' data(mouse.time.mpse)
##' mouse.time.mpse %>% mp_rrarefy()
setGeneric("mp_rrarefy", function(.data, .abundance=NULL, raresize, trimOTU=FALSE, trimSample=FALSE, seed=123, ...){standardGeneric("mp_rrarefy")}) 

.internal_mp_rrarefy <- function(.data, raresize, trimOTU=FALSE, trimSample=FALSE, seed=123, ...){
    if (missing(raresize) || is.null(raresize)){
        raresize <- min(rowSums(.data))
    }
    if (is.na(seed) || is.null(seed)){
        message("To reproduce and it was not be set, the random seed will set to 123 automatically.")
        seed <- 123
    }
    if (trimSample){
        .data <- .data[rowSums(.data) >= raresize, ,drop=FALSE]
    }
    res <- withr::with_seed(seed, vegan::rrarefy(x=.data, sample=raresize)) %>% suppressWarnings()
    if (trimOTU){
        removeOTU <- colSums(res)==0
        if (sum(removeOTU) > 0){
            message(sum(removeOTU), " OTUs were removed because they are no longer present in any sample after ",
                    "rarefaction, if you want to keep them you can set 'trimOTU = FALSE' !")
        }
        res <- res[, !removeOTU]
    }
    return(res)
}

##' @rdname mp_rrarefy-methods
##' @aliases mp_rrarefy,MPSE
##' @exportMethod mp_rrarefy
setMethod("mp_rrarefy", signature(.data="MPSE"), function(.data, .abundance=NULL, raresize, trimOTU=FALSE, trimSample=FALSE, seed=123, ...){
    .abundance <- rlang::enquo(.abundance)
    if (rlang::quo_is_null(.abundance)){
        .abundance <- as.symbol("Abundance")
    }
    allassays <- SummarizedExperiment::assays(.data) %>% as.list()
    if ("RareAbundance" %in% names(allassays)){
        message("The RareAbundance was in the MPSE object, please check whether it has been rarefied !")
        return(.data)
    }
    rare <- .internal_mp_rrarefy(.data=t(allassays[[rlang::as_name(.abundance)]]), raresize=raresize, trimOTU=FALSE, trimSample=FALSE, seed=123) %>% t()
    SummarizedExperiment::assays(.data)@listData <- c(allassays, list(RareAbundance=rare))
    if (trimOTU){
        removeOTU <- rowSums(rare)==0
        if (sum(removeOTU) > 0){
            message(sum(removeOTU), " OTUs were removed because they are no longer present in any sample after ",
                    "rarefaction, if you want to keep them you can set 'trimOTU = FALSE' !")
            .data <- .data[!rownames(.data) %in% rownames(rare[removeOTU, , drop=FALSE]), ]
        }
    }

    if (trimSample){
        removeSample <- colSums(rare) < max(colSums(rare))
        if (sum(removeSample) > 0){
            message(paste0(sum(removeSample), ifelse(sum(removeSample)>1, " samples were", " samples was"), 
                           " removed because they do not have enough abundance when raresize = ", raresize,
                           ", if you want to keep them you can set 'trimSample = FALSE' !"))
            .data <- .data[, !colnames(.data) %in% colnames(rare[, removeSample, drop=FALSE])]
        }
    }
    return(.data)
})

##' @rdname mp_rrarefy-methods
##' @aliases mp_rrarefy,tbl_mpse
##' @exportMethod mp_rrarefy
setMethod("mp_rrarefy", signature(.data="tbl_mpse"), function(.data, .abundance=NULL, raresize, trimOTU=FALSE, trimSample=FALSE, seed=123, ...){
    .abundance <- rlang::enquo(.abundance)
    if (rlang::quo_is_null(.abundance)){
        .abundance <- as.symbol("Abundance")
    }
    assaysvar <- attr(.data, "assaysvar")
    otutree <- attr(.data, "otutree")
    taxatree <- attr(.data, "taxatree")
    if ("RareAbundance" %in% assaysvar){
        message("The RareAbundance was in the MPSE object, please check whether it has been rarefied !")
        return(.data)
    }
    tmpotu <- .data %>% mp_extract_assays(.abundance=!!.abundance, byRow=FALSE)
    rare <- .internal_mp_rrarefy(.data=tmpotu, raresize=raresize, trimOTU=FALSE, trimSample=FALSE, seed=seed, ...) %>%
            t() %>% 
            tibble::as_tibble(rownames = "OTU") %>% 
            tidyr::pivot_longer(!"OTU", names_to="Sample", values_to="RareAbundance")
    othernms <- colnames(.data)[!colnames(.data) %in% c("OTU", "Sample", assaysvar)]
    res <- .data %>% left_join(rare, by=c("OTU", "Sample"), suffix=c("", ".y")) %>% 
           select(c("OTU", "Sample", assaysvar, "RareAbundance", othernms))
    if (trimOTU){
        removeOTU <- res %>% 
                     dplyr::group_by(.data$OTU) %>% 
                     dplyr::summarise(Total=sum(.data$RareAbundance)) %>% 
                     filter(.data$Total==0) %>% 
                     pull("OTU") %>%
                     unique()
        if (length(removeOTU) > 0){
            message(length(removeOTU), " OTUs were removed because they are no longer present in any sample after ",
                    "rarefaction, if you want to keep them you can set 'trimOTU = FALSE' !")
            res %<>% dplyr::filter(!.data$OTU %in% removeOTU)
        }
    }
    if (trimSample){
        sampleTotal <- res %>% 
                       dplyr::group_by(.data$Sample) %>%
                       dplyr::summarise(Total=sum(.data$RareAbundance)) %>%
                       dplyr::pull(var=.data$Total, name=.data$Sample)

        removeSample <- sampleTotal < max(sampleTotal)
        if (sum(removeSample) > 0){
            message(paste0(sum(removeSample), ifelse(sum(removeSample)>1, " samples were", " samples was"), 
                           " removed because they do not have enough abundance when raresize = ", raresize,
                           ", if you want to keep them you can set 'trimSample = FALSE' !"))
            res %<>% dplyr::filter(!.data$Sample %in% names(removeSample[removeSample]))
        }
    }
    res <- add_attr.tbl_mpse(x1 = res, x2 = .data)
    attr(res, "assaysvar") <- c(assaysvar, "RareAbundance")
    return (res)
})

##' @rdname mp_rrarefy-methods
##' @aliases mp_rrarefy,grouped_df_mpse
##' @exportMethod mp_rrarefy
setMethod("mp_rrarefy", signature(.data="grouped_df_mpse"), function(.data, .abundance=NULL, raresize, trimOTU=FALSE, trimSample=FALSE, seed=123, ...){
    .abundance <- rlang::enquo(.abundance)
    if (rlang::quo_is_null(.abundance)){
        .abundance <- as.symbol("Abundance")
    }
    tmpgroups <- attr(.data, "groups")
    groupvars <- names(tmpgroups)[names(tmpgroups) != ".rows"]
    groupvars <- lapply(groupvars, function(i) as.symbol(i))
    .data %<>% ungroup
    res <- mp_rrarefy(.data=.data, .abundance=!!.abundance, raresize=raresize, trimOTU=trimOTU, trimSample=trimSample, seed=seed, ...)
    res <- do.call(group_by, c(list(res), groupvars))
    return (res)
})

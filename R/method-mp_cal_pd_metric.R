#' @title calculating related phylogenetic alpha metric 
#' @param obj object, data.frame of (nrow sample * ncol taxonomy(feature)) 
#' or phyloseq.
#' @param mindepth numeric, Subsample size for rarefying community.
#' @param sampleda data.frame, sample information, row sample * column factors.
#' @param tree tree object, it can be phylo object or treedata object.
#' @param metric the related phylogenetic metric, options is 'NRI', 'NTI', 'PD', 'PAE', 'HAED', 'EAED', 'IAC', 'all',
#' default is 'PAE', meaning all the metrics ('NRI', 'NTI', 'PD', 'PAE', 'HAED', 'EAED', 'IAC').
#' @param abundance.weighted logical, whether calculate mean nearest taxon distances for each species 
#' weighted by species abundance, default is FALSE.
#' @param force logical whether calculate the index even the count of otu is
#' not rarefied, default is FALSE. If it is TRUE, meaning the rarefaction is not be
#' performed automatically.
#' @param seed integer a random seed to make the result reproducible, default is 123.
#' @param ... additional arguments, meaningless now.
#' @return alphasample object contained NRT and NTI.
#' @rdname get_NRI_NTI-methods
#' @author Shuangbin Xu
#' @export
setGeneric("get_NRI_NTI", function(obj, ...){standardGeneric("get_NRI_NTI")})

#' @aliases get_NRI_NTI,matrix
#' @rdname get_NRI_NTI-methods
#' @export
setMethod("get_NRI_NTI", "matrix", function(obj, mindepth, sampleda, tree, metric=c('PAE', 'NRI', 'NTI', 'PD', 'HAED', 'EAED', 'IAC', 'all'),
                                            abundance.weighted=FALSE, force=FALSE, seed=123, ...){
    if (missing(mindepth) || is.null(mindepth)){
           mindepth <- min(rowSums(obj))
    }
    if (obj %>% rowSums %>% var != 0 && !force){
        obj <- vegan::rrarefy(obj, mindepth)
    }
    metric %<>% match.arg(c('NRI', 'NTI', 'PD', 'PAE', 'HAED', 'EAED', 'IAC', 'all'))
    if (metric == 'all'){
        res <- .internal_cal_all_pd_metric(obj, tree, weighted.abund = abundance.weighted, seed = seed, ...)
    }

    res <- switch(metric,
        NRI = withr::with_seed(seed, .internal_cal_nri(obj, cal_treedist(tree), weighted.abund = abundance.weighted, ...)),
        NTI = withr::with_seed(seed, .internal_cal_nti(obj, cal_treedist(tree), weighted.abund = abundance.weighted, ...)),
        PD  = .internal_cal_pd(obj, tree),
        PAE = .internal_cal_pae(obj, tree),
        HAED = .internal_cal_haed(obj, tree, ...),
        EAED = .internal_cal_eaed(obj, tree, ...),
        IAC = .internal_cal_iac(obj, tree),
        all = res
    )  
    if (metric != 'all'){
        res <- data.frame(res) %>% 
               magrittr::set_colnames(metric)
    }
    if (missing(sampleda)){
        sampleda <- NULL
    }
    res <- new("alphasample",
               alpha=res,
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

#' @title Calculating related phylogenetic alpha metric with MPSE or tbl_mpse object
#'
#' @param .data object, MPSE or tbl_mpse object
#' @param .abundance The column name of OTU abundance column to be calculate.
#' @param action character it has three options, "add" joins the new information
#' to the input tbl (default), "only" return a non-redundant tibble with the just
#' new information, ang 'get' return a 'alphasample' object.
#' @param metric the related phylogenetic metric, options is 'NRI', 'NTI', 'PD', 'PAE', 'HAED', 'EAED', 'IAC', 'all',
#' default is 'PAE', 'all' meaning all the metrics ('NRI', 'NTI', 'PD', 'PAE', 'HAED', 'EAED', 'IAC').
#' @param abundance.weighted logical, whether calculate mean nearest taxon distances for each species
#' weighted by species abundance, default is TRUE.
#' @param force logical whether calculate the alpha index even the '.abundance' is
#' not rarefied, default is FALSE.
#' @param seed integer a random seed to make the result reproducible, default is 123.
#' @param ... additional arguments see also "ses.mpd" and "ses.mntd" of "picante".
#' @return update object.
#' @rdname mp_cal_pd_metric-methods
#' @author Shuangbin Xu
#' @export
#'
#' @references Cadotte, M.W., Jonathan Davies, T., Regetz, J., Kembel, S.W., Cleland,
#' E. and Oakley, T.H. (2010), Phylogenetic diversity metrics for ecological communities:
#' integrating species richness, abundance and evolutionary history. Ecology Letters,
#' 13: 96-105. https://doi.org/10.1111/j.1461-0248.2009.01405.x.
#' 
#' Webb, C. O. (2000). Exploring the phylogenetic structure of ecological communities: 
#' an example for rain forest trees. The American Naturalist, 156(2), 145-155.
#' https://doi.org/10.1086/303378.
#' @examples
#' \dontrun{
#'   suppressPackageStartupMessages(library(curatedMetagenomicData))
#'   xx <- curatedMetagenomicData('ZellerG_2014.relative_abundance', dryrun=F)
#'   xx[[1]] %>% as.mpse -> mpse
#'   mpse %<>% 
#'     mp_cal_pd_metric(
#'       .abundance = Abundance, 
#'       force = TRUE,
#'       metric = 'PAE'
#'     )
#'   mpse %>% 
#'     mp_plot_alpha(
#'       .alpha = PAE,
#'       .group = disease
#'   )
#' }
setGeneric("mp_cal_pd_metric", function(
    .data, 
    .abundance, 
    action = "add", 
    metric = c('PAE', 'NRI', 'NTI', 'PD', 'HAED', 'EAED', 'all'),
    abundance.weighted = FALSE, 
    force = FALSE, 
    seed = 123, 
    ...)
standardGeneric("mp_cal_pd_metric")
)

#' @rdname mp_cal_pd_metric-methods
#' @aliases mp_cal_pd_metric,MPSE
#' @exportMethod mp_cal_pd_metric
setMethod("mp_cal_pd_metric", signature(.data = "MPSE"), function(
    .data, 
    .abundance, 
    action = "add", 
    metric = c('PAE', 'NRI', 'NTI', 'PD', 'HAED', 'EAED', 'all'),
    abundance.weighted = FALSE, 
    force = FALSE, 
    seed = 123, 
    ...){
    action %<>% match.arg(c("add", "only", "get"))
    metric <- rlang::enquo(metric)
    metric <- quo.vect_to_str.vect(metric)    
    metric %<>% match.arg(c('NRI', 'NTI', 'PD', 'PAE', 'HAED', 'EAED', 'IAC', 'all'))
    res <- .internal_cal_pd_metric(.data=.data,
                                 .abundance=!!rlang::enquo(.abundance),
                                 abundance.weighted=abundance.weighted,
                                 metric = metric,
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
          left_join(da, by="Sample", suffix=c("", ".y")) #%>%
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
          
.internal_cal_pd_metric <- function(.data, .abundance, abundance.weighted = FALSE, metric = 'NRI', force=FALSE, seed=123, ...){

    .abundance <- rlang::enquo(.abundance)
    if (rlang::quo_is_null(.abundance)){
        .abundance <- as.symbol("RareAbundance")
    }

    otutree <- .data %>%
               mp_extract_tree(type="otutree") %>%
               suppressMessages()

    if (is.null(otutree)){
        rlang::abort("The otutree is required to calculate the NRI and NTI !")
    }

    if (!valid_rare(.data, .abundance=.abundance) && !force){
        trash <- try(silent = TRUE,
                     expr = {
                         .data <- mp_rrarefy(.data = .data, ...)
                     }
                 )
        if (inherits(trash, "try-error")){
            stop_wrap("The 'Abundance' column cannot be rarefied, please check whether it is integer (count).
                       Or you can set 'force=TRUE' to calculate the result of 'NRI' and 'NTI' without rarefaction.
                      ")
        }

        message_wrap("The rarefied abundance of species might not be provided. Rarefaction of all
                      observations is performed automatically using 'Abundance' column. If you still
                      want to calculate the result of 'NRI' and 'NTI' with the specified '.abundance',
                      you can set 'force=TRUE'. ")
        .abundance <- as.symbol("RareAbundance")           
    }

    indexda <- .data %>% 
               mp_extract_assays(.abundance=!!.abundance, byRow=FALSE) %>%
               as.matrix() %>%
               get_NRI_NTI(tree=otutree@phylo, force=TRUE, metric = metric, abundance.weighted=abundance.weighted, seed=seed, ...)

    return(indexda)
}

.internal_cal_pd_metric_ <- function(
    .data, 
    .abundance, 
    action = "add", 
    metric = c('PAE', 'NRI', 'NTI', 'PD', 'HAED', 'EAED', 'all'), 
    abundance.weighted = TRUE, 
    force = FALSE, 
    seed = 123, 
    ...){
    action %<>% match.arg(c("add", "only", "get"))
    metric <- rlang::enquo(metric)
    metric <- quo.vect_to_str.vect(metric) 
    metric %<>% match.arg(c('NRI', 'NTI', 'PD', 'PAE', 'HAED', 'EAED', 'IAC', 'all'))
    res <- .internal_cal_pd_metric(.data = .data,
                                 .abundance = !!rlang::enquo(.abundance),
                                 abundance.weighted = abundance.weighted, 
                                 metric = metric,
                                 force = force,
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
          left_join(da, by="Sample", suffix=c("", ".y")) #%>%
          #column_to_rownames(var="Sample")
    if (action=="add"){
        samplevar <- .data %>% attr("samplevar")
        assaysvar <- .data %>% attr("assaysvar")
        othernm <- colnames(.data)[!colnames(.data) %in% c("OTU", "Sample", assaysvar, samplevar)]
        .data %<>% left_join(da, by="Sample", suffix=c("", ".y")) %>%
                   select(c("OTU", "Sample", assaysvar, samplevar, colnames(da), othernm))
        return(.data)
    }else if (action=="only"){
        return(da)
    }
}

#' @rdname mp_cal_pd_metric-methods
#' @aliases mp_cal_pd_metric,tbl_mpse
#' @exportMethod mp_cal_pd_metric
setMethod("mp_cal_pd_metric", signature(.data="tbl_mpse"), .internal_cal_pd_metric_)

#' @rdname mp_cal_pd_metric-methods
#' @aliases mp_cal_pd_metric,grouped_df_mpse
#' @exportMethod mp_cal_pd_metric
setMethod("mp_cal_pd_metric", signature(.data="grouped_df_mpse"), .internal_cal_pd_metric_)

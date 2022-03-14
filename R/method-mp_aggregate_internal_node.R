#' calculate the mean/median (relative) abundance of internal nodes according to their children tips.
#' @rdname mp_aggregate_internal_node-methods
#' @param .data MPSE object which must contain otutree slot, required
#' @param .abundance the column names of abundance.
#' @param force logical whether calculate the (relative) abundance forcibly when the abundance
#' is not be rarefied, default is FALSE.
#' @param relative logical whether calculate the relative abundance.
#' @param aggregate_fun function the method to calculate the (relative) abundance of internal nodes
#' according to their children tips, default is mean.
#' @param action character, "add" joins the new information to the otutree slot if it exists (default). 
#' In addition, "only" return a non-redundant tibble with the just new information. "get" return 'otutree'
#' slot, which is a treedata object.
#' @param ... additional parameters, meaningless now.
#' @return a object according to 'action' argument.
#' @export
#' @examples
#' \dontrun{
#'   library(curatedMetagenomicData) 
#'   xx <- curatedMetagenomicData('ZellerG_2014.relative_abundance', dryrun=F)
#'   xx[[1]] %>% as.mpse -> mpse
#'   otu.tree <- mpse %>% 
#'     mp_aggregate_internal_node(
#'       .abundance = Abundance, 
#'       force = TRUE, 
#'       relative = FALSE,
#'       action = 'get' # other option is 'add' or 'only'.
#'     )
#'   otu.tree
#' }
setGeneric("mp_aggregate_internal_node", 
           function(.data, .abundance = NULL, force = FALSE, relative = TRUE, aggregate_fun = mean, action='get', ...)
               standardGeneric("mp_aggregate_internal_node")
)

.aggregate_internal_node <- function(.data, .abundance = NULL, force = FALSE, relative = TRUE, aggregate_fun = mean, action = 'get', ...){
    otu.tree <- .data %>% mp_extract_otutree()
    if (is.null(otu.tree)){
        stop_wrap("The object did not contain otutree slot.")
    }
    action %<>% match.arg(c("add", "get", "only"))

    .abundance <- rlang::enquo(.abundance)

    if (rlang::quo_is_null(.abundance)){
        .abundance <- as.symbol("RareAbundance")
    }

    if (!valid_rare(.data, .abundance=.abundance) && !force){
        trash <- try(silent = TRUE,
                     expr = {
                         .data <- mp_rrarefy(.data = .data, ...)
                     }
                 )
        if (inherits(trash, "try-error")){
            stop_wrap("The 'Abundance' column cannot be rarefied, please check whether it is integer (count).
                       Or you can set 'force=TRUE' to aggregate the abundance without rarefaction.
                      ")
        }

        message_wrap("The rarefied abundance of species might not be provided. Rarefaction of all
                      observations is performed automatically using 'Abundance' column. If you still
                      want to aggregate with the specified '.abundance', you can set 'force=TRUE'. ")
        .abundance <- as.symbol("RareAbundance")
    }
    da <- .data %>% mp_extract_assays(.abundance = !!.abundance)
    if (relative){
        da <- apply(da, 2, function(x)x/sum(x) * 100)
        .abundance <- as.symbol(paste0("Rel", rlang::as_name(.abundance)))
    }
    da %<>% 
        tibble::as_tibble(rownames='OTU') %>% 
        tidyr::pivot_longer(
          cols = !.data$OTU, 
          names_to = 'Sample', 
          values_to = rlang::as_name(.abundance)
        ) %>%
        dplyr::left_join(
          otu.tree %>% 
            select('node', 'label'),
          by=c('OTU'='label')
        ) %>%
        select(-'OTU')
    if (is.null(otu.tree@phylo$node.label)){
        otu.tree@phylo <- ape::makeNodeLabel(otu.tree@phylo)
    }
    
    sample.da <- .data %>% mp_extract_sample() %>% remove_MP_internal_res()
    index.name <- paste0(rlang::as_name(.abundance), 'BySample')
    inode.da <- treeio::offspring(
                   otu.tree, 
                   .node = otu.tree %>% 
                           dplyr::filter(!.data$isTip, keep.td = FALSE) %>% 
                           pull(.data$node), 
                   type='tips'
                ) %>%
                lapply(
                   .internal_cal_inode, 
                   x = da, 
                   tree = otu.tree, 
                   fun = aggregate_fun,
                   abundance = rlang::as_name(.abundance)
                ) %>% 
                dplyr::bind_rows(.id = 'node') %>%
                dplyr::left_join(sample.da, by = "Sample") %>%
                tidyr::nest(!!rlang::sym(index.name) := - .data$node) %>%
                dplyr::mutate_at("node", as.integer)
    da %<>% 
        dplyr::left_join(sample.da, by = 'Sample') %>%
        tidyr::nest(!!rlang::sym(index.name) := - .data$node) %>%
        dplyr::bind_rows(inode.da)
    if (action == 'only'){
        da <- otu.tree %>% select('node', 'label') %>% dplyr::left_join(da, by='node')
        return (da)
    }
    if (index.name %in% treeio::get.fields(otu.tree)){
        otu.tree %<>% dplyr::select(- index.name, keep.td = TRUE)
    }
    otu.tree %<>% left_join(da, by='node')

    if (action == 'get'){
        return(otu.tree)
    }else if (action == 'add'){
        if (inherits(.data, 'MPSE')){
            otutree(.data) <- otu.tree
        }else{
            attr(.data, 'otutree') <- otu.tree 
        }
        return(.data)
    }
}

.internal_cal_inode <- function(node2tips, x, tree, fun, abundance){
    x %<>% dplyr::filter(.data$node %in% node2tips) %>% 
          dplyr::group_by(.data$Sample) %>%
          dplyr::summarize(dplyr::across(abundance, fun, .name=abundance))
    return(x)
}

#' @rdname mp_aggregate_internal_node-methods
#' @aliases mp_aggregate_internal_node,MPSE
#' @exportMethod mp_aggregate_internal_node
setMethod("mp_aggregate_internal_node", signature(.data = "MPSE"), .aggregate_internal_node)

#' @rdname mp_aggregate_internal_node-methods
#' @aliases mp_aggregate_internal_node,tbl_mpse
#' @exportMethod mp_aggregate_internal_node
setMethod("mp_aggregate_internal_node", signature(.data = 'tbl_mpse'), .aggregate_internal_node)

#' @rdname mp_aggregate_internal_node-methods
#' @aliases mp_aggregate_internal_node,grouped_df_mpse
#' @exportMethod mp_aggregate_internal_node
setMethod("mp_aggregate_internal_node", signature(.data = 'grouped_df_mpse'), .aggregate_internal_node)

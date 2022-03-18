#' Calculating the balance score of internal nodes (clade) according to the mean/median abundance of their binary children tips.
#'
#' @rdname mp_balance_clade-methods
#' @param .data MPSE object which must contain otutree slot, required
#' @param .abundance the column names of abundance.
#' @param force logical whether calculate the (relative) abundance forcibly when the abundance
#' is not be rarefied, default is FALSE.
#' @param relative logical whether calculate the relative abundance.
#' @param aggregate_fun function the method to calculate the (relative) abundance of internal nodes
#' according to their children tips, default is 'mean', other options are 'median' and 'geometric.mean'.
#' @param pseudonum numeric add a pseudo numeric to avoid the error of division in calculation, default 
#' is 0.001 .
#' @param action character, "add" joins the new information to the otutree slot if it exists (default).
#' In addition, "only" return a non-redundant tibble with the just new information. "get" return 'otutree'
#' slot, which is a treedata object.
#' @param ... additional parameters, meaningless now.
#' @return a object according to 'action' argument.
#' @export
setGeneric("mp_balance_clade",
           function(.data, 
                    .abundance = NULL, 
                    force = FALSE, 
                    relative = TRUE, 
                    aggregate_fun = c('mean', 'median', 'geometric.mean'), 
                    pseudonum = .001, 
                    action='get', ...)
               standardGeneric("mp_balance_clade")
)

.balance_clade <- function(.data,
                           .abundance = NULL,
                           force = FALSE, 
                           relative = TRUE, 
                           aggregate_fun = c('mean', 'median', 'geometric.mean'), 
                           pseudonum = .001, action = 'get', ...){
    aggregate_fun %<>% match.arg(c('mean', 'median', 'geometric.mean'))
    otu.tree <- .data %>% mp_extract_otutree() %>% suppressMessages()
    if (is.null(otu.tree)){
        warning_wrap("The object did not contain otutree slot.")
        return (NULL)
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
                       Or you can set 'force=TRUE' to calculate the balance of clade without rarefaction.
                      ")
        }

        message_wrap("The rarefied abundance of species might not be provided. Rarefaction of all
                      observations is performed automatically using 'Abundance' column. If you still
                      want to calculate the balance of clade with the specified '.abundance', you can 
                      set 'force=TRUE'. ")
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
    #if (is.character(aggregate_fun)){
    #    aggregate_fun <- rlang::as_function(aggregate_fun)
    #}
    sample.da <- .data %>% mp_extract_sample() %>% remove_MP_internal_res()
    index.name <- paste0('BalanceBy', rlang::as_name(.abundance), 'BySample')
    inodes <- otu.tree %>% .extract_nodes()
     
    inodes2binary <- extract_binary_offspring(otu.tree, inodes)
    if (!is_binary_tree(inodes2binary)){
        stop_wrap("The otutree is not a binary tree")
    }

    da <- lapply(
         inodes2binary,
         .internal_balance_clade,
         x = da,
         fun = aggregate_fun,
         abundance = rlang::as_name(.abundance),
         pseudonum = pseudonum
      ) %>%
      dplyr::bind_rows(.id = 'node') %>%
      dplyr::left_join(sample.da, by = "Sample") %>%
      tidyr::nest(!!rlang::sym(index.name) := - .data$node) %>%
      dplyr::mutate_at("node", as.integer)

    if (action == 'only'){
        da <- otu.tree %>% 
              dplyr::filter(!.data$isTip, keep.td = FALSE) %>%
              select('node', 'label') %>%
              dplyr::left_join(da, by='node')
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

#' @rdname mp_balance_clade-methods
#' @aliases mp_balance_clade,MPSE
#' @exportMethod mp_balance_clade
setMethod("mp_balance_clade", signature(.data = "MPSE"), .balance_clade)

#' @rdname mp_balance_clade-methods
#' @aliases mp_balance_clade,tbl_mpse
#' @exportMethod mp_balance_clade
setMethod("mp_balance_clade", signature(.data = 'tbl_mpse'), .balance_clade)

#' @rdname mp_balance_clade-methods
#' @aliases mp_balance_clade,grouped_df_mpse
#' @exportMethod mp_balance_clade
setMethod("mp_balance_clade", signature(.data = 'grouped_df_mpse'), .balance_clade)

##' extract the binary offspring of the specified internal nodes
##' @param .data phylo or treedata object
##' @param .node the internal nodes
##' @param type the type of binary offspring, options are 'tips' (default),
##' 'all', 'internal'.
##' @param ... additional parameter, meaningless now.
##' @export
extract_binary_offspring <- function(.data, .node, type = 'tips', ...){
    UseMethod('extract_binary_offspring')
}

##' @method extract_binary_offspring phylo
##' @export
extract_binary_offspring.phylo <- function(.data, .node, type = 'tips', ...){
    res <- lapply(.node,
                  function(i).internl_extract_binary(.data = .data, .node = i, type = type)
           ) %>%
           stats::setNames(.node)
    return(res)
}

##' @method extract_binary_offspring treedata
##' @export
extract_binary_offspring.treedata <- function(.data, .node, type = 'all', ...) {
    extract_binary_offspring(
         .data = as.phylo(.data),
         .node,
         type = type,
         ...)
}

.internl_extract_binary <- function(.data, .node, type = 'all'){
    child.nodes <- treeio::child(.data = .data, .node = .node)
    res <- lapply(child.nodes, function(i)treeio::offspring(.data, .node = i, tiponly = TRUE)) %>%
        suppressMessages() %>%
        stats::setNames(child.nodes)

    index <- lapply(res, length)==0
    res[index] <- as.numeric(names(res[index]))
    return(res)
}

.internal_balance_clade <- function(node2binary, x, fun, abundance, pseudonum=.001){
    newnm <- paste0('BalanceBy', abundance)
    res <- lapply(node2binary, function(i){
            x %>% 
            dplyr::filter(.data$node %in% i) %>%
            dplyr::group_by(.data$Sample) %>%
            dplyr::summarize(Abun = .internal.aggregate(
                   x = .data[[abundance]],
                   fun = fun, 
                   pseudonum = pseudonum 
              )
            )
         }
    )
    
    res2 <- data.frame(Sample=res[[1]][,1,drop=TRUE]) 
    res2[[newnm]] <- log((res[[1]][,2,drop=TRUE] + pseudonum) / (res[[2]][,2,drop=TRUE] + pseudonum))
    res2 %<>% select('Sample', newnm)
    return(res2)
}

is_binary_tree <- function(x){
    all(lapply(x, length) == 2)
}

geometric.mean <- function(x, pseudonum, na.rm = TRUE){
    y <- exp(mean(log(x + pseudonum), na.rm = na.rm)) - pseudonum
    return(y)
}

.internal.aggregate <- function(x, fun, pseudonum = 0.01){
    xx <- switch(fun,
        mean = mean(x, na.rm=TRUE),
        median = median(x, na.rm=TRUE),
        geometric.mean = geometric.mean(x, pseudonum=pseudonum)
    )
    return(xx)
}

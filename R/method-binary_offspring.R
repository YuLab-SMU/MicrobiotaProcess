#' Calculating the balance score of internal nodes (clade) according to the mean/median abundance of their binary children tips.
#'
#' @rdname mp_balance_clade-methods
#' @param .data MPSE object which must contain otutree slot, required
#' @param .abundance the column names of abundance.
#' @param force logical whether calculate the (relative) abundance forcibly when the abundance
#' is not be rarefied, default is FALSE.
#' @param relative logical whether calculate the relative abundance.
#' @param balance_fun function the method to calculate the (relative) abundance of internal nodes
#' according to their children tips, default is 'geometric.mean', other options are 'mean' and 'median'.
#' @param pseudonum numeric add a pseudo numeric to avoid the error of division in calculation, default 
#' is 0.001 .
#' @param action character, "add" joins the new information to the otutree slot if it exists (default).
#' In addition, "only" return a non-redundant tibble with the just new information. "get" return a new 'MPSE' 
#' object, and the 'OTU' column is the internal nodes and 'Abundance' column is the balance scores.
#' @param ... additional parameters, meaningless now.
#' @return a object according to 'action' argument.
#' @export
#' @references Morton JT, Sanders J, Quinn RA, McDonald D, Gonzalez A, Vázquez-Baeza Y, 
#' Navas-Molina JA, Song SJ, Metcalf JL, Hyde ER, Lladser M, Dorrestein PC, 
#' Knight R. 2017. Balance trees reveal microbial niche differentiation. 
#' mSystems 2:e00162-16. https://doi.org/10.1128/mSystems.00162-16.
#'
#' Justin D Silverman, Alex D Washburne, Sayan Mukherjee, Lawrence A David. 
#' A phylogenetic transform enhances analysis of compositional microbiota data. 
#' eLife 2017;6:e21887. https://doi.org/10.7554/eLife.21887.001.
#'
#' @examples
#' \dontrun{
#'   suppressPackageStartupMessages(library(curatedMetagenomicData))
#'   xx <- curatedMetagenomicData('ZellerG_2014.relative_abundance', dryrun=F)
#'   xx[[1]] %>% as.mpse -> mpse
#'   mpse.balance.clade <- mpse %>%
#'     mp_balance_clade(
#'       .abundance = Abundance,
#'       force = TRUE,
#'       relative = FALSE,
#'       action = 'get',
#'       pseudonum = .01
#'     )
#'   mpse.balance.clade 
#'
#'   # Performing the Euclidean distance or PCA.
#'
#'   mpse.balance.clade %>%
#'     mp_cal_dist(.abundance = Abundance, distmethod = 'euclidean') %>%
#'     mp_plot_dist(.distmethod = 'euclidean', .group = disease, group.test = T)
#'
#'   mpse.balance.clade %>%
#'     mp_adonis(.abundance = Abundance, distmethod = 'euclidean', permutation = 9999)
#' 
#'   mpse.balance.clade %>%
#'     mp_cal_pca(.abundance = Abundance) %>% 
#'     mp_plot_ord(.group = disease)
#'
#'   # Detecting the signal balance nodes.
#'   mpse.balance.clade %>% mp_diff_analysis(
#'       .abundance = Abundance,
#'       force = TRUE,
#'       relative = FALSE,
#'       .group = disease,
#'       fc.method = 'compare_mean'
#'   )
#' }
setGeneric("mp_balance_clade",
           function(.data, 
                    .abundance = NULL, 
                    force = FALSE, 
                    relative = TRUE, 
                    balance_fun = c('geometric.mean', 'mean', 'median'), 
                    pseudonum = .001, 
                    action='get', ...)
               standardGeneric("mp_balance_clade")
)

.balance_clade <- function(.data,
                           .abundance = NULL,
                           force = FALSE, 
                           relative = TRUE, 
                           balance_fun = c('geometric.mean', 'mean', 'median'), 
                           pseudonum = .001, action = 'get', ...){
    balance_fun %<>% match.arg(c('geometric.mean', 'mean', 'median'))
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
        .abundance <- as.symbol(paste0("Rel", rlang::as_name(.abundance), 'BySample'))
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
    #if (is.character(balance_fun)){
    #    balance_fun <- rlang::as_function(balance_fun)
    #}
    sample.da <- .data %>% mp_extract_sample() %>% remove_MP_internal_res()
    index.name <- paste0('BalanceBy', rlang::as_name(.abundance))
    inodes <- otu.tree %>% .extract_nodes()
     
    inodes2binary <- extract_binary_offspring(otu.tree, inodes)
    if (!is_binary_tree(inodes2binary)){
        stop_wrap("The otutree is not a binary tree")
    }

    da <- lapply(
         inodes2binary,
         .internal_balance_clade,
         x = da,
         fun = balance_fun,
         abundance = rlang::as_name(.abundance),
         pseudonum = pseudonum
      ) %>%
      dplyr::bind_rows(.id = 'node')

    node2offspring <- .balance_offspring_node(inodes2binary, otu.tree)

    if (action == 'get'){
        node2label <- otu.tree %>% dplyr::select(.data$node, .data$label, keep.td=FALSE)
        da %<>% tidyr::pivot_wider(id_cols='node', names_from='Sample', values_from=index.name) %>%
            dplyr::mutate_at("node", as.integer) %>%
            dplyr::left_join(y = node2label, by = 'node') %>%
            dplyr::select(-.data$node) %>% 
            tibble::column_to_rownames(var='label')
        
        mpse <- MPSE(assays=list(Abundance=da), colData=sample.da %>% tibble::column_to_rownames(var='Sample'))
        mpse %<>% dplyr::left_join(node2label %>% 
                                   left_join(node2offspring, by='node'), 
                                   by=c('OTU'='label'))
        message_wrap('The new MPSE object (internal nodes as features) will be created.')
        return(mpse)
    }

    da %<>%
      dplyr::left_join(sample.da, by = "Sample") %>%
      tidyr::nest(!!rlang::sym(index.name) := - .data$node) %>%
      dplyr::mutate_at("node", as.integer)

    if (action == 'only'){
        da <- otu.tree %>% 
              dplyr::filter(!.data$isTip, keep.td = FALSE) %>%
              select('node', 'label') %>%
              dplyr::left_join(da, by='node') %>%
              dplyr::left_join(node2offspring, by='node')
        return (da)
    }

    if (index.name %in% treeio::get.fields(otu.tree)){
        otu.tree %<>% dplyr::select(- index.name, keep.td = TRUE)
    }
    otu.tree %<>% left_join(da, by='node') %>% left_join(node2offspring, by='node')

    if (action == 'add'){
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
extract_binary_offspring.treedata <- function(.data, .node, type = 'tips', ...) {
    extract_binary_offspring(
         .data = as.phylo(.data),
         .node,
         type = type,
         ...)
}

.internl_extract_binary <- function(.data, .node, type = 'tips'){
    child.nodes <- treeio::child(.data = .data, .node = .node)
    res <- lapply(child.nodes, function(i)treeio::offspring(.data, .node = i, tiponly = TRUE, type=type)) %>%
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
    weight.s <- .weight.subtree.count(node2binary)
    res2 <- data.frame(Sample=res[[1]][,1,drop=TRUE]) 
    res2[[newnm]] <- weight.s * log((res[[1]][,2,drop=TRUE] + pseudonum) / (res[[2]][,2,drop=TRUE] + pseudonum))
    res2 %<>% select('Sample', newnm)
    return(res2)
}

is_binary_tree <- function(x){
    all(lapply(x, length) == 2)
}

geometric.mean <- function(x, pseudonum = 0, na.rm = TRUE){
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

.weight.subtree.count <- function(x){
    x <- lapply(x, length)
    x <- do.call(`*`, x) / do.call(`+`, x)
    return(sqrt(x))
}

.balance_offspring_node <- function(x, otu.tree){
    otu.tree %<>% select(.data$node, .data$label)
    x <- lapply(x,function(i)lapply(i, function(j)otu.tree %>% dplyr::filter(.data$node %in% j) %>% pull(.data$label)))
    #x %<>% purrr::map_depth(2,function(j)otu.tree %>% dplyr::filter(.data$node %in% j) %>% pull(.data$label))
    res <- x %>% lapply(., function(x) do.call('cbind', list(x)) %>% 
                 as.data.frame() %>% 
                 setNames('offspringTiplabel') %>% 
                 tibble::rownames_to_column(var='Children') %>%
                 dplyr::mutate_at('Children', as.integer) %>% 
                 dplyr::mutate(Clade=c('down', 'up'))) %>% 
           dplyr::bind_rows(.id='node') %>% 
           tibble::as_tibble() %>%
           dplyr::mutate(pseudolabel=unlist(lapply(.data$offspringTiplabel, function(x)paste0(unlist(x), collapse=';')))) %>%
           dplyr::group_by(.data$node) %>%
           dplyr::mutate(pseudolabel=paste0(.data$node, ":", paste0(.data$pseudolabel, collapse='/'))) %>%
           dplyr::ungroup() %>% 
           dplyr::mutate_at('node', as.integer) %>%
           tidyr::nest(Balance_offspring=-.data$node)
    return(res)
}

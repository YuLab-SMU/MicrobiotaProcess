#' Calculating the balance score of internal nodes (clade) according to the geometric.mean/mean/median abundance of their binary children tips.
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
#'     mp_adonis(.abundance = Abundance, .formula=~disease, distmethod = 'euclidean', permutation = 9999)
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
    if (!ape::is.binary(otu.tree@phylo)){
        warning_wrap("The otutree is not a binary tree, it would be converted to 
                     binary tree using multi2di of ape automatically!")
        otu.tree@phylo <- ape::multi2di(otu.tree@phylo)
        if (!ape::is.binary(otu.tree@phylo)){
            stop_wrap('The multi2di can not convert the otutree to a binary tree automatically, 
                      please convert it to binary tree manually !')
        }
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
    if (any(da<0)){
        stop_wrap('The ', rlang::as_name(.abundance), ' can not contain negative value.')
    }
    da <- da + pseudonum
    
    if (is.null(otu.tree@phylo$node.label)){
        otu.tree@phylo <- ape::makeNodeLabel(otu.tree@phylo)
    }
    sample.da <- .data %>% mp_extract_sample() %>% remove_MP_internal_res()
    index.name <- paste0('BalanceBy', rlang::as_name(.abundance))
    inodes <- .nodeId(otu.tree, type = 'internal')

    tip2ances <- .convert_balance_tip2ancestors(otu.tree@phylo, inodes)
    inode.balance <- .internal_balance(x = da, tip2ances = tip2ances, fun = balance_fun)
    
    node2offspring <- .balance_offspring_node(tip2ances, inodes, otu.tree)

    if (action == 'get'){
        node2label <- otu.tree %>% dplyr::select(.data$node, .data$label, keep.td=FALSE)
        node2name <- treeio::nodeid(otu.tree, otu.tree@phylo$node.label)
        names(node2name) <- otu.tree@phylo$node.label
        node2name <- node2name[match(colnames(inode.balance), node2name)]
        colnames(inode.balance) <- names(node2name)
        mpse <- MPSE(assays=list(Abundance=t(inode.balance)), colData=sample.da %>% tibble::column_to_rownames(var='Sample'))
        mpse %<>% dplyr::left_join(node2label %>% 
                                   left_join(node2offspring, by='node'), 
                                   by=c('OTU'='label'))
        message_wrap('The new MPSE object (internal nodes as features) will be created.')
        return(mpse)
    }

    inode.balance %<>%
        as_tibble(rownames = 'Sample') %>%
        tidyr::pivot_longer(
          cols = !.data$Sample,
          names_to = 'node',
          values_to = index.name
        ) %>%
        dplyr::left_join(sample.da, by = 'Sample') %>%
        tidyr::nest(!!rlang::sym(index.name) := - .data$node) %>%
        dplyr::mutate_at("node", as.integer)    

    if (action == 'only'){
        da <- otu.tree %>% 
              dplyr::filter(!.data$isTip, keep.td = FALSE) %>%
              select('node', 'label') %>%
              dplyr::left_join(inode.balance, by='node') %>%
              dplyr::left_join(node2offspring, by='node')
        return (da)
    }

    if (index.name %in% treeio::get.fields(otu.tree)){
        otu.tree %<>% dplyr::select(- index.name, keep.td = TRUE)
    }
    otu.tree %<>% left_join(inode.balance, by='node') %>% left_join(node2offspring, by='node')

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
                  function(i).internal_extract_binary(.data = .data, .node = i, type = type)
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

.internal_extract_binary <- function(.data, .node, type = 'tips'){
    child.nodes <- treeio::child(.data = .data, .node = .node)
    res <- lapply(child.nodes, function(i)treeio::offspring(.data, .node = i, tiponly = TRUE, type=type)) %>%
        suppressMessages() %>%
        stats::setNames(child.nodes)

    index <- lapply(res, length)==0
    res[index] <- as.numeric(names(res[index]))
    return(res)
}

#.internal_balance_clade <- function(node2binary, x, fun, abundance, pseudonum=.001){
#    newnm <- paste0('BalanceBy', abundance)
#    res <- lapply(node2binary, function(i){
#            x %>% 
#            dplyr::filter(.data$node %in% i) %>%
#            dplyr::group_by(.data$Sample) %>%
#            dplyr::summarize(Abun = .internal.aggregate(
#                   x = .data[[abundance]],
#                   fun = fun, 
#                   pseudonum = pseudonum 
#              )
#            )
#         }
#    )
#    weight.s <- .weight.subtree.count(node2binary)
#    res2 <- data.frame(Sample=res[[1]][,1,drop=TRUE]) 
#    res2[[newnm]] <- weight.s * log((res[[1]][,2,drop=TRUE] + pseudonum) / (res[[2]][,2,drop=TRUE] + pseudonum))
#    res2 %<>% select('Sample', newnm)
#    return(res2)
#}
# 
# is_binary_tree <- function(x){
#     all(lapply(x, length) == 2)
# }

geometric.mean <- function(x, pseudonum = 0, na.rm = TRUE, nozero=FALSE){
    if (all(x == 0)){
        return(0)
    }
    if (nozero){
        y <- exp(mean(log(x[x>0]), na.rm = na.rm))
    }else{
        y <- exp(mean(log(x + pseudonum), na.rm = na.rm)) - pseudonum
    }
    return(y)
}

#.internal.aggregate <- function(x, fun, pseudonum = 0.01){
#    xx <- switch(fun,
#        mean = mean(x, na.rm=TRUE),
#        median = median(x, na.rm=TRUE),
#        geometric.mean = geometric.mean(x, pseudonum=pseudonum)
#    )
#    return(xx)
#}
#
#.weight.subtree.count <- function(x){
#    x <- lapply(x, length)
#    x <- do.call(`*`, x) / do.call(`+`, x)
#    return(sqrt(x))
#}

.balance_offspring_node <- function(x, inode, tree){
    children <- lapply(inode, function(i)treeio::child(tree, i)) %>%
                stats::setNames(inode) %>%
                dplyr::bind_cols() %>%
                dplyr::mutate(Clade=c('down', 'up')) %>%
                tidyr::pivot_longer(cols = !.data$Clade, names_to='node', values_to='Children')
    offspringtips <- x %>% tibble::as_tibble(rownames='offspringTiplabel') %>%
        tidyr::pivot_longer(cols = !.data$offspringTiplabel, names_to='node', values_to='Clade') %>%
        dplyr::filter(!is.na(.data$Clade)) %>%
        dplyr::mutate(Clade=ifelse(.data$Clade==1, 'down', 'up')) %>%
        dplyr::group_by(.data$node, .data$Clade) %>%
        dplyr::mutate(offspringTiplabel = list(.data$offspringTiplabel)) %>%
        dplyr::distinct() %>%
        dplyr::ungroup(.data$Clade) %>%
        dplyr::mutate(pseudolabel=paste0(.data$node, ':', paste0(lapply(.data$offspringTiplabel, paste0, collapse=';'), collapse='/')))
    res <- dplyr::left_join(children, offspringtips, by=c('node'='node', 'Clade'='Clade')) %>%
           dplyr::mutate_at('node', as.integer) %>%
           tidyr::nest(Balance_offspring = -.data$node)
    return(res)    
}

.convert_balance_tip2ancestors <- function(tree, inode){
    tmp <- treeio::nodeid(tree, tree$tip.label)
    names(tmp) <- tree$tip.label
    res <- lapply(inode, function(x){
        children <- treeio::child(tree, x)
        unlist(mapply(.build_up_and_down,
               children, c(1, -1),
               MoreArgs=list(tree = tree),
               SIMPLIFY = FALSE
        ))
    }) %>%
    stats::setNames(inode)
    res <- as.matrix(dplyr::bind_rows(res))
    rownames(res) <- inode
    colnames(res) <- names(tmp[match(colnames(res), tmp)])
    return (t(res))
}

.build_up_and_down <- function(x, y, tree){
    newx <- treeio::offspring(tree, x, tiponly = TRUE)
    if (length(newx)==0){
        newx <- x
    }
    res <- rep(y, length(newx))
    names(res) <- newx
    return(res)
}

.internal_balance <- function(x, tip2ances, fun){
    x <- x[match(rownames(tip2ances), rownames(x)), ,drop=FALSE]
    res <- lapply(seq_len(ncol(x)), function(i){
       res <- x[, i] * tip2ances
       res <- switch(fun,
            mean = apply(res, 2, .balance_mean, na.rm = TRUE),
            median = apply(res, 2, .balance_median, na.rm = TRUE),
            geometric.mean = apply(res, 2, .balance_geometric.mean, na.rm = TRUE)
          )
       }
    ) %>%
    do.call('rbind', .)
    rownames(res) <- colnames(x)
    return(res)
}

.balance_mean <- function(x, pseudonum = 0, na.rm = TRUE){
    if (na.rm) x <- x[!is.na(x)]
    weights <- .weight.subtree(x)
    weights * log(mean(x[x > 0]) / mean(abs(x[x<0])))
}

.balance_median <- function(x, pseudonum = 0, na.rm = TRUE){
    if (na.rm) x <- x[!is.na(x)]
    weights <- .weight.subtree(x)
    weights * log(median(x[x > 0]) / median(abs(x[x<0])))
}

.balance_geometric.mean <- function(x, pseudonum = 0, na.rm = TRUE){
    if (na.rm) x <- x[!is.na(x)]
    weights <- .weight.subtree(x)
    weights * log(geometric.mean(x[x > 0])  / geometric.mean(abs(x[x<0])))
}

.weight.subtree <- function(x){
    i <- sum(x>0)
    j <- sum(x<0)
    sqrt((i * j) / (i + j))
}

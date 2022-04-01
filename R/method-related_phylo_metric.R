.internal_cal_all_pd_metric <- function(x, tree, weighted.abund = FALSE, seed = 123, ...){
    tip.dis <- cal_treedist(tree)
    params <- list(...)
    if ("method" %in% names(params)){
        HAED = .internal_cal_haed(x, tree, method = params$method)
    }else{
        HAED = .internal_cal_haed(x, tree)
    }
    res <- list(
        NRI = withr::with_seed(seed, .internal_cal_nri(x, tip.dis, weighted.abund = weighted.abund, ...)),
        NTI = withr::with_seed(seed, .internal_cal_nti(x, tip.dis, weighted.abund = weighted.abund, ...)),
        PD = .internal_cal_pd(x, tree),
        PAE = .internal_cal_pae(x, tree),
        HAED = HAED,
        EAED = HAED / log(rowSums(x)) 
    ) 
    res <- do.call(cbind, res)
    return (res)
}


.internal_cal_pae <- function(x, tree, flag.zero = TRUE){
    edge <- tree$edge
    istip <- ! edge[,2] %in% edge[,1]
    tips.edge.length <- tree$edge.length[istip]
    x <- x[, match(tree$tip.label, colnames(x))]
    xx <- x == 0
    PD <- .internal_cal_pd(xx, tree)
    TL.x <- .generate_tip.len(xx, tree, tips.edge.length, flag.zero = flag.zero)
    sumabun.tip <- rowSums(TL.x * (x - 1))
    meanabun.tip <- (rowSums(x)/colSums(apply(x, 1, function(i)i>0)) - 1) * rowSums(TL.x)
    pae <- (PD + sumabun.tip)/(PD + meanabun.tip)
    #pae <- as.data.frame(pae) %>% stats::setNames(nm='PAE')
    return (pae)
}

.generate_tip.len <- function(x, tree, tips.edge.length, flag.zero = TRUE){
    if (!all(is.logical(x))){
        x <- x == 0
    }
    if (flag.zero && any(x)){
        res <- lapply(seq_len(nrow(x)), function(i)
                   ape::drop.tip(tree, names(x[i,][x[i,]])) %>% 
                   tidytree::as_tibble() %>% 
                   dplyr::mutate(isTip=!.data$node %in% .data$parent) %>% 
                   dplyr::filter(.data$isTip) %>% 
                   dplyr::pull(.data$branch.length, name=.data$label)) %>% 
        dplyr::bind_rows() 
        res[is.na(res)] <- 0
        tmp <- setdiff(colnames(x), colnames(res))
        if (length(tmp) > 0){
            tmp.x <- data.frame(matrix(0, ncol=length(tmp), nrow=nrow(x)))
            colnames(tmp.x) <- tmp
            res <- cbind(res, tmp.x)
        }
    }else{
        res <- t((!t(x)) * tips.edge.length)
    }
    return(res[,match(tree$tip.label, colnames(res))])
}

.internal_cal_nri <- function(x, tree, weighted.abund = FALSE, bootnum = 999){
    if (inherits(tree, 'phylo') || inherits(tree, 'treedata')){
        tip.dis <- cal_treedist(tree)
    }else{
        tip.dis <- tree
    }
    mpd.obs <- .internal_mpd(x, tip.dis, weighted.abund = weighted.abund)
    mpd.random <- replicate(bootnum, .internal_mpd(x, .make_random_label(tip.dis), weighted.abund = weighted.abund))
    mpd.mean <- apply(mpd.random, 1, mean, na.rm = TRUE)
    mpd.sd <- apply(mpd.random, 1, sd, na.rm = TRUE)
    nri <- - (mpd.obs - mpd.mean) / mpd.sd
    return(nri)

}

.internal_mpd <- function(x, tree, weighted.abund = FALSE){
    if (inherits(tree, 'phylo') || inherits(tree, 'treedata')){
        tip.dis <- cal_treedist(tree)
    }else{
        tip.dis <- tree
    }
    res <- lapply(seq_len(nrow(x)), function(i){
        sp <- names(x[i,][x[i,] != 0])
        if (length(sp) > 0){
            tmp.dis <- tip.dis[sp, sp]
            if (weighted.abund){
                abund <- x[i, sp, drop = FALSE]
                weighted.abund <- t(as.matrix(abund)) %*% abund
                res <- stats::weighted.mean(tmp.dis, weighted.abund) 
            }else{
                res <- mean(tmp.dis[lower.tri(tmp.dis)]) 
            }
        }else{
            res <- NA
        }
        return(res)
    }) %>% unlist() 
    names(res) <- rownames(x)
    return(res)
}

.make_random_label <- function(x){
    rownames(x) <- colnames(x) <- sample(rownames(x))
    return(x)
}

.internal_cal_nti <- function(x, tree, weighted.abund = FALSE, bootnum = 999){
    if (inherits(tree, 'phylo') || inherits(tree, 'treedata')){
        tip.dis <- cal_treedist(tree)
    }else{
        tip.dis <- tree
    }
    mntd.obs <- .internal_mntd(x, tip.dis, weighted.abund = weighted.abund)
    mntd.random <- replicate(bootnum, .internal_mntd(x, .make_random_label(tip.dis), weighted.abund = weighted.abund))
    mntd.mean <- apply(mntd.random, 1, mean, na.rm = TRUE)
    mntd.sd <- apply(mntd.random, 1, sd, na.rm = TRUE)
    nti <- - (mntd.obs - mntd.mean) / mntd.sd
    return(nti)
}

.internal_mntd <- function(x, tree, weighted.abund = FALSE){
    if (inherits(tree, 'phylo') || inherits(tree, 'treedata')){
        tip.dis <- cal_treedist(tree)
    }else{
        tip.dis <- tree
    }
    res <- lapply(seq_len(nrow(x)), function(i){
        sp <- names(x[i,][x[i,] != 0])
        if (length(sp) > 0){
            tmp.dis <- tip.dis[sp, sp]
            diag(tmp.dis) <- NA
            tmp.mntd <- apply(tmp.dis, 1, min, na.rm = TRUE)
            if (weighted.abund){
                res <-  stats::weighted.mean(tmp.mntd, x[i, sp])
            }else{
                res <- mean(tmp.mntd)
            }
        }else{
            res <- NA
        }
        return(res)
    }) %>% unlist()    
    names(res) <- rownames(x)
    return(res)
}

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

.internal_cal_haed <- function(x, tree, method='multi.abund', ...){
    #AED <- .internal_cal_aed(x = x, tree = tree, action = 'get')
    AED <- .internal_cal_aed(x = x, tree = tree)
    PD <- .internal_cal_pd(x, tree)
    if (method == 'multi.abund' || !all.equal(round(x), x)){
        tmp <- as.matrix(AED) * x / PD
        #tmp <- apply(tmp,1,function(i)i[i>0])
        tmp2 <- log(tmp)
        tmp2[is.infinite(as.matrix(tmp2))] <- 0
        res <- - rowSums(tmp * tmp2)
        #tmp <- lapply(seq_len(nrow(x)), function(i){
        #           unname(as.matrix(AED[i,])) * unname(as.vector(x[i, names(AED[i,])]))/PD[[i]]
        #         }
        #       )
    }else{
        tmp <- lapply(seq_len(nrow(x)), function(i){
                   rep(unname(as.matrix(AED[i,])), unname(as.vector(x[i, names(AED[i,])])))/PD[[i]]
                 }
               )
        res <- lapply(tmp, function(i)-sum(i * log(i))) %>% unlist() %>%
            stats::setNames(rownames(x))
    } 
    return(res)
}

.internal_cal_eaed <- function(x, tree){
    Haed <- .internal_cal_haed(x, tree)
    Eaed <- Haed/log(rowSums(x))
}

.internal_cal_aed <- function(x, tree, action='only', ...){
    #x <- x + pseudonum
    tree <- .internal_drop_zero_tip(x, tree)
    tips2ancestor <- lapply(tree, .convert_tips2ancestors_sbp)
    edge.len <- lapply(tree, .extract_edge)
    tips2ancestor <- mapply(function(x, y){x[,match(names(y), colnames(x))]}, 
                            tips2ancestor, edge.len, SIMPLIFY=FALSE)
    if (length(tree) != nrow(x)){
        tips2ancestor <- rep(tips2ancestor, nrow(x))
        edge.len <- rep(edge.len, nrow(x))
    }
    res <- lapply(seq_len(nrow(x)), function(i){
        sp <- x[i, match(rownames(tips2ancestor[[i]]), colnames(x))]
        xx <- sp * tips2ancestor[[i]]
        AEDi <- tcrossprod(x=edge.len[[i]], y=prop.table(xx, margin=2))
        AEDi <- AEDi / sp
        apply(AEDi, 2, function(i)i)
      }
    )
    if (action == 'get'){
        names(res) <- rownames(x)
        return(res)
    }
    
    res %<>% 
    dplyr::bind_rows() %>%
    as.data.frame() 
    res[is.na(res)] <- 0
    rownames(res) <- rownames(x)
    tmp <- setdiff(colnames(x), colnames(res))
    if (length(tmp)>0){
        tmp.x <- data.frame(matrix(0, ncol=length(tmp), nrow=nrow(x)))
        colnames(tmp.x) <- tmp
        res <- cbind(res, tmp.x)
    }
    res <- res[,match(colnames(x), colnames(res))]
    return(res)
}

.convert_tips2ancestors_sbp <- function(tree, include.root = FALSE){
    all.nodes <- .nodeId(tree)
    if (!include.root){
        all.nodes <- setdiff(all.nodes, treeio::rootnode(tree))
    }   
    tip.nodes <- .nodeId(tree, type = 'tips')
    if(inherits(tree, "treedata")){
        tiplabels <- tree@phylo$tip.label
    }else{
        tiplabels <- tree$tip.label
    }
    lapply(tip.nodes, treeio::ancestor, .data=tree) %>% 
        mapply(append, tip.nodes, ., SIMPLIFY=FALSE) %>% 
        lapply(.,function(i) all.nodes %in% i) %>%
        stats::setNames(tip.nodes) %>%
        do.call(rbind, .) %>%
        magrittr::set_colnames(all.nodes) %>%
        magrittr::set_rownames(tiplabels)
}

.extract_edge <- function(tree, type='all', include.root = FALSE){
    if (inherits(tree, 'treedata')){
        tree <- tree@phylo
    }
    node2edge <- tree$edge.length %>% stats::setNames(tree$edge[, 2])
    if (include.root){
        rootid <- treeio::rootnode(tree)
        node2edge <- c(node2edge, tree$root.edge %>% stats::setNames(rootid))
    }
    node2edge <- node2edge[order(names(node2edge))]
    if (type == 'all'){
        return(node2edge)
    }
    nodes <- .nodeId(tree, type = type)
    return(node2edge[match(nodes, names(node2edge))])
}

.nodeId <- function(tree, type='all'){
    type <- match.arg(type, c('all', 'tips', 'internal'))
    if (inherits(tree, 'treedata')){
        tree <- tree@phylo
    }
    nodes <- unique(as.vector(tree$edge))
    if (type == 'all'){
        return(nodes)
    }
    edge <- tree$edge
    
    tips <- edge[!edge[,2] %in% edge[,1], 2]
    if (type == 'tips'){
        return(tips)
    }else if(type == 'internal'){
        return(setdiff(nodes, tips))
    }
}

.internal_cal_pd <- function(x, tree){
    if (inherits(tree, 'treedata')){
        tree <- tree@phylo    
    }
    tree <- .internal_drop_zero_tip(x, tree)
    pd <- unlist(lapply(tree, function(i)sum(i$edge.length)))
    if (length(pd)!=nrow(x)){
        pd <- rep(pd, nrow(x))
    }
    names(pd) <- rownames(x)

    return(pd)
}

.internal_drop_zero_tip <- function(x, tree){
    if (!all(is.logical(x))){
        x <- x == 0 
    }
    if (any(x)){
        tree <- lapply(seq_len(nrow(x)), function(i)
                   ape::drop.tip(tree, names(x[i,][x[i,]])))
    }else{
        tree <- list(tree)
    }
    return(tree)
}

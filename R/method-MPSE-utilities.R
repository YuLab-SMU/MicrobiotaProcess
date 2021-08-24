convert_to_treedata2 <- function(x){
    x <- x %>%
        tibble::add_column(Root="r__root", .before=1) %>%
        dplyr::mutate(OTU=rownames(.))
    datalist <- list()
    clnm <- colnames(x)
    for (i in seq_len(ncol(x)-1)){
        tmpdat <- x[,c(i, i+1)]
        colnames(tmpdat) <- c("parent", "child")
        tmpdat %<>% dplyr::mutate(nodeClass=clnm[i+1], nodeDepth=i) %>%
                    dplyr::distinct()
        datalist[[i]] <- tmpdat
    }
    datalist <- do.call(rbind, datalist)
    isTip <- !as.vector(datalist$child) %in% as.vector(datalist$parent)
    index <- rep(NA, length(isTip))
    index[isTip] <- seq(1,sum(isTip))
    index[!isTip] <- seq(sum(isTip) + 2, length(isTip) + 1)
    mapping <- data.frame(node=index, labelnames=as.vector(datalist$child), isTip)
    indxx <- match(mapping$labelnames, datalist$child)
    mapping$nodeClass <- datalist[indxx, "nodeClass"]
    mapping$nodeDepth <- datalist[indxx, "nodeDepth"]
    parentnode <- mapping[match(datalist$parent, mapping$labelnames),"node"]
    childnode <- mapping[match(datalist$child, mapping$labelnames),"node"]
    edges <- cbind(parentnode, childnode) %>% as.matrix() 
    colnames(edges) <- NULL
    edges[is.na(edges)] <- sum(isTip) + 1
    root <- data.frame(node=sum(isTip)+1, labelnames="r__root",
                       isTip=FALSE, nodeClass="Root", nodeDepth=0)
    mapping <- rbind(mapping, root)
    mapping <- mapping[order(mapping$node),]
    node.label <- as.vector(mapping$labelnames)[!mapping$isTip]
    tip.label <- as.vector(mapping$labelnames)[mapping$isTip]
    mapping <- mapping[,colnames(mapping) %in% c("node", "nodeClass", "nodeDepth")]
    taxphylo <- structure(list(edge=edges, node.label=node.label,
                               tip.label=tip.label,
                               Nnode = length(node.label)),
                          class="phylo")
    res <- treeio::treedata(phylo=taxphylo, data=as_tibble(mapping))
    res
}

taxatree_to_tb <- function(x){
    x %<>% as_tibble(x)
    extrada <- x %>% select(-c("parent", "node", "nodeDepth")) %>%
               dplyr::filter(.data$nodeClass=="OTU") %>% select(-c("nodeClass"))
    clnm <- x %>% dplyr::select("nodeClass", "nodeDepth") %>%
            dplyr::distinct() %>% dplyr::arrange(.data$nodeDepth) %>%
            dplyr::select("nodeClass") %>% unlist(use.names=FALSE)

    x$parent <- x[match(x$parent, x$node),]$label
    d <- x %>%
         dplyr::filter(.data$nodeDepth!=0) %>%
         dplyr::select("parent", "label", "nodeDepth") %>%
         dplyr::group_split(.data$nodeDepth) %>%
         as.list() %>%
         purrr::map(select, -"nodeDepth")
    for (i in seq_len(length(d))){
        d[[i]] %<>% setNames(c(clnm[i], clnm[i+1]))
    }
    d %<>% purrr::reduce(dplyr::right_join) %>%
        suppressMessages() %>%
        select(-c("Root")) 
    if (ncol(extrada)>1){
        d %<>% dplyr::left_join(extrada, by=c("OTU"="label"), suffix=c("", ".y"))
    }    
    d %<>% column_to_rownames(var="OTU")
    return (d)   
}

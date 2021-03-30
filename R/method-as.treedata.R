#' @title convert dataframe contained hierarchical relationship or other classes to treedata class
#' @param data data.frame, such like the tax_table of phyloseq.
#' @param type character, the type of datasets, default is "species", if the dataset is not about species,                                                                                                         #' such as dataset of kegg function, you should set it to "others".
#' @param ..., additional parameters.
#' @return treedata class.
#' @author Shuangbin Xu
#' @importFrom tibble as_tibble
#' @importFrom tidytree treedata
#' @export
#' @examples
#' data(hmp_aerobiosis_small)
#' head(taxda)
#' treedat <- convert_to_treedata(taxda)
convert_to_treedata <- function(data, type="species", ...){
    if (!"fillNA" %in% names(attributes(data))){
        data <- fillNAtax(data, type=type)
    }
    data <- data.frame(root=rep("r__root", nrow(data)), data)
    datalist <- list()
    for (i in seq_len(ncol(data)-1)){
    	tmpdat <- data[,c(i, i+1)]
    	colnames(tmpdat) <- c("parent", "child")
    	datalist[[i]] <- tmpdat
    }
    datalist <- do.call("rbind", datalist)
    datalist <- datalist[!duplicated(datalist),]
    isTip <- !as.vector(datalist$child) %in% as.vector(datalist$parent)
    index <- c()
    index[isTip] <- seq(1,sum(isTip))
    index[!isTip] <- seq(sum(isTip)+2,length(isTip)+1)
    mapping <- data.frame(node=index, labelnames=as.vector(datalist$child), isTip)
    mapping$nodeClass <- unlist(lapply(as.vector(mapping$labelnames),
    		   		       function(x)(unlist(strsplit(x,"__"))[1])))
    mapping$nodeSize <- 1
    parentnode <- mapping[match(as.vector(datalist$parent), as.vector(mapping$labelnames)),]$node 
    childnode <- mapping[match(as.vector(datalist$child), as.vector(mapping$labelnames)),]$node
    edges <- cbind(parentnode, childnode) 
    colnames(edges) <- NULL
    edges[is.na(edges)] <- sum(isTip) + 1
    root <- data.frame(node=sum(isTip)+1,labelnames="r__root",
    		       isTip=FALSE, nodeClass="r", nodeSize=1)
    mapping <- rbind(root, mapping)
    mapping <- mapping[order(mapping$node),]
    node.label <- as.vector(mapping$labelnames)[!mapping$isTip]
    tip.label <- as.vector(mapping$labelnames)[mapping$isTip]
    mapping <- mapping[,colnames(mapping) %in% c("node", "nodeClass", "nodeSize")]
    taxphylo <- structure(list(edge=edges, node.label=node.label,
                               tip.label=tip.label, edge.length=rep(0.5, nrow(edges)),
                               Nnode = length(node.label)), class="phylo")
    res <- treedata(phylo=taxphylo, data=as_tibble(mapping))
}

#' convert taxonomyTable to treedata
#'
#' @title as.treedata
#' @param tree object, This is for taxonomyTable class, 
#' so it should be a taxonomyTable.
#' @param ... additional parameters. 
#' @method as.treedata taxonomyTable
#' @rdname as.treedata
#' @export
#' @examples
#' data(test_otu_data)
#' tree <- as.treedata(phyloseq::tax_table(test_otu_data))
as.treedata.taxonomyTable <- function(tree, ...){
    convert_to_treedata(data.frame(tree, check.names=FALSE))
}


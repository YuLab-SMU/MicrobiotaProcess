#' @importFrom tidytree as.treedata
#' @export
tidytree::as.treedata

#' @title convert dataframe contained hierarchical relationship or other classes to treedata class
#' @method as.treedata data.frame
#' @rdname as.treedata
#' @author Shuangbin Xu
#' @importFrom tibble as_tibble
#' @export
as.treedata.data.frame <- function(data,...){
	data <- fillNAtax(data)
	data <- data.frame(root=rep("r__root", nrow(data)), data)
	datalist <- list()
	for (i in 1:(ncol(data)-1)){
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
	mapping <- data.frame(node=index, 
						  labelnames=as.vector(datalist$child), 
						  isTip)
	mapping$nodeClass <- unlist(lapply(as.vector(mapping$labelnames),
									   function(x)(unlist(strsplit(x,"__"))[1])))
	#if (is.null(nodeSize)){
	mapping$nodeSize <- 1
	#}
	parentnode <- mapping[match(as.vector(datalist$parent), as.vector(mapping$labelnames)),]$node 
	childnode <- mapping[match(as.vector(datalist$child), as.vector(mapping$labelnames)),]$node
	edges <- cbind(parentnode, childnode) 
	colnames(edges) <- NULL
	edges[is.na(edges)] <- sum(isTip) + 1 
	root <- data.frame(node=sum(isTip)+1,
					   labelnames="r__root",
					   isTip=FALSE,
					   nodeClass="r",
					   nodeSize=1)
	mapping <- rbind(root, mapping)
	mapping <- mapping[order(mapping$node),]
	node.label <- as.vector(mapping$labelnames)[!mapping$isTip]
	tip.label <- as.vector(mapping$labelnames)[mapping$isTip]
	taxphylo <- structure(list(edge=edges,
							   node.label=node.label,
							   tip.label=tip.label,
							   edge.length=rep(0.5, nrow(edges)),
							   Nnode = length(node.label)), 
						  class="phylo") 
	res <- new("treedata",
			   phylo=taxphylo, 
			   data=as_tibble(mapping))

}

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
    #mapping$nodeSize <- 1
    parentnode <- mapping[match(as.vector(datalist$parent), as.vector(mapping$labelnames)),]$node 
    childnode <- mapping[match(as.vector(datalist$child), as.vector(mapping$labelnames)),]$node
    edges <- cbind(parentnode, childnode) 
    colnames(edges) <- NULL
    edges[is.na(edges)] <- sum(isTip) + 1
    root <- data.frame(node=sum(isTip)+1,labelnames="r__root",
    		       isTip=FALSE, nodeClass="r")#, nodeSize=1)
    mapping <- rbind(root, mapping)
    mapping <- mapping[order(mapping$node),]
    node.label <- as.vector(mapping$labelnames)[!mapping$isTip]
    tip.label <- as.vector(mapping$labelnames)[mapping$isTip]
    mapping <- mapping[,colnames(mapping) %in% c("node", "nodeClass")]#, "nodeSize")]
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


#' @method as.treedata tbl_ps
#' @export
as.treedata.tbl_ps <- function(tree, use_taxatree=TRUE, tiplevel="OTU", ...){
    tr <- attr(tree, "tree")
    taxavar <- attr(tree, "taxavar")
    if (!is.null(tr) && !use_taxatree){
        treeda <- tr %>% as_tibble()
        extrada <- tree %>% nest()
        tiplevel <- "OTU"
    }else{
        if (use_taxatree && !is.null(taxavar)){
            if (!is.null(taxavar)){
                if (tiplevel != "OTU"){tiplevel <- taxavar[length(taxavar)]}
                taxavar <- taxavar[taxavar!=tiplevel]
                taxavar2 <- c(taxavar, tiplevel)
                if (tiplevel != "OTU"){
                    taxada <- tree %>% select(taxavar2) %>% distinct() %>% 
                               column_to_rownames(var="OTU") 
                    treeda1 <- taxada %>% data.frame() %>% fillNAtax()
                    taxavar <- taxavar[taxavar != "OTU"]
                }else{
                    treeda1 <- tree %>% select(taxavar2) %>% distinct() %>% 
                               data.frame() %>% fillNAtax()
                }
                treeda <- convert_to_treedata(data=treeda1) %>% as_tibble()
                extrada <- tree %>% select(-taxavar) %>% nest()
                if (tiplevel == "OTU"){
                    isTip <- !treeda$node %in% treeda$parent
                    treeda$label[isTip] <- restore_name(treeda$label[isTip])
                }else{
                    extrada[[tiplevel]] <- as.vector(treeda1[match(extrada$OTU, rownames(treeda1)),tiplevel])
                    extrada <- extrada %>% select(-c("OTU")) %>% nest() 
                }
            }else{
                stop("The tax table is empty in the object!")
            }
        }else{
            stop("The tree slot is empty, you can use the taxa tree via set use_taxatree=TRUE")
        }
    }
    treeda %<>% full_join(extrada, by=c("label"=tiplevel)) %>% as.treedata()
    return(treeda)
}

restore_name <- function(x){
    x <- strsplit(x, "__")
    x <- lapply(x, function(x){
          if (length(x)==1){
              res <- x
          }
          if (length(x)==2){
              res <- x[2]
          }
          if (length(x)>2){
              res <- paste(x[-1], collapse="__")
          }
          return (res)
     }
    ) %>% unlist()
    return (x)
}

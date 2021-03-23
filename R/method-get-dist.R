#' @title calculate distance
#' @param obj phyloseq, phyloseq class or data.frame
#' nrow sample * ncol feature. 
#' @param method character, default is hellinger, 
#' see alse \code{\link[vegan]{decostand}} 
#' @param distmethod character, default is "euclidean", 
#' see also \code{\link[phyloseq]{distanceMethodList}}
#' @param taxa_are_rows logical, default is FALSE.
#' @param sampleda data.frame, nrow sample * ncol factors.
#' @param tree object, the phylo class, see also \code{\link[ape]{as.phylo}}.
#' @param ..., additional parameters.
#' @return distance class contianed distmethod and originalD attr
#' @export
#' @examples
#' data(test_otu_data)
#' distclass <- get_dist(test_otu_data)
#' hcsample <- get_clust(distclass)
get_dist <- function(obj,...){
    UseMethod("get_dist")
}

#' @method get_dist data.frame
#' @rdname get_dist
#' @importFrom vegan decostand
#' @importFrom phyloseq otu_table 
#' @export
get_dist.data.frame <- function(obj, 
    distmethod="euclidean",
    taxa_are_rows=FALSE,	
    sampleda=NULL,
    tree=NULL,
    method="hellinger",
    ...){
    objphyloseq <- new("phyloseq",
                       otu_table=otu_table(obj, 
                       taxa_are_rows=taxa_are_rows),
                       sam_data=phyloseq::sample_data(sampleda),
                       phy_tree=tree)
    return(get_dist.phyloseq(objphyloseq, 
                             distmethod=distmethod, 
                             method=method,
                             ...))
    
}

#' @method get_dist phyloseq
#' @importFrom phyloseq distance taxa_are_rows phy_tree
#' @seealso \code{\link[phyloseq]{distance}}
#' @rdname get_dist
#' @export
get_dist.phyloseq <- function(obj, distmethod="euclidean", method="hellinger",...){
    tmpmethod <- gsub("^(u.*)*unifrac$", "unifrac", distmethod, ignore.case = TRUE)
    tmpmethod <- gsub("^w.*unifrac$", "wunifrac", distmethod, ignore.case = TRUE) 
    tree <- obj@phy_tree
    if (tmpmethod=="unifrac" || tmpmethod=="wunifrac"){
    	if(is.null(tree)){
    		stop("The tree should be required when the distmethod is `WeightUniFrac` or `UnWeightUniFrac`")
    	}
    }
    if (!is.null(method)){
    	if (taxa_are_rows(obj@otu_table)){
    		tmpotu <- t(obj@otu_table)
    	}else{
    		tmpotu <- data.frame(obj@otu_table)
    	}
    	obj@otu_table <- otu_table(decostand(tmpotu, method=method), 
                                   taxa_are_rows=FALSE)
    }
    disres <- distance(obj, method=distmethod, type="sample", ...)
    attr(disres, "distmethod") <- distmethod
    attr(disres, "originalD") <- data.frame(obj@otu_table, check.names=FALSE)
    return(disres)
}

#' @title calculate distance
#' @export
getdist <- function(obj,...){
	UseMethod("getdist")
}

#' @method getdist default
#' @param data data.frame, nrow sample * ncol feature.
#' @param method character, default is hellinger, see alse \code{\link[vegan]{decostand}}
#' @param distmethod character, default is "euclidean", see also \code{\link[phyloseq]{distanceMethodList}}
#' @param taxa_are_rows logical, default is FALSE.
#' @param sampleda data.frame, nrow sample * ncol factors.
#' @param tree phylo, see also \code{\link[ape]{phylo}}
#' @rdname getdist
#' @importFrom vegan decostand
#' @importFrom phyloseq otu_table 
#' @export
getdist.default <- function(data, 
							distmethod="euclidean",
						    taxa_are_rows=FALSE,	
							sampleda=NULL,
							tree=NULL,
							type="sample",
							method="hellinger",
							...){
	tmpmethod <- gsub("^(u.*)*unifrac$", "unifrac", distmethod, ignore.case = TRUE)
	tmpmethod <- gsub("^w.*unifrac$", "wunifrac", distmethod, ignore.case = TRUE)
	if (method=="unifrac" || method=="wunifrac"){
		if(is.null(tree)){
			stop("The tree should be required when the method is `WeightUniFrac` or `UnWeightUniFrac`")
		}
	}
	objphyloseq <- new("phyloseq",
					   otu_table=otu_table(data, 
									  taxa_are_rows=taxa_are_rows),
					   sam_data=sampleda,
					   phy_tree=tree)

	return(getdist.phyloseq(objphyloseq, 
							distmethod=tmpmethod, 
							method=method,
							type=type, ...))
	
}

#' @method getdist phyloseq
#' @importFrom phyloseq distance taxa_are_rows
#' @seealso \code{\link[phyloseq]{distance}}
#' @rdname getdist
#' @export
getdist.phyloseq <- function(obj, distmethod="euclidean", type="sample", method="hellinger",...){
	if (!is.null(method)){
		if (taxa_are_rows(obj@otu_table)){
			tmpotu <- t(obj@otu_table)
		}else{
			tmpotu <- data.frame(obj@otu_table)
		}
		obj@otu_table <- otu_table(decostand(tmpotu, method=method), 
								   taxa_are_rows=FALSE)
	}
	disres <- distance(obj, method=distmethod, type=type, ...)
	attr(disres, "distmethod") <- distmethod
	return(disres)
}

#' @title performs principal coordinate analysis (PCoA)
#' @rdname getpcoa
#' @export
getpcoa <- function(obj, ...){
	UseMethod("getpcoa")
}

#' @method getpca default
#' @rdname getpcoa
#' @export
getpcoa.default <- function(data, 
							distmethod="euclidean", 
							taxa_are_rows=FALSE,
							sampleda=NULL, 
							tree=NULL,
							type="sample",
							...){
	tmpdist <- getdist.default(data=data, 
							   distmethod=distmethod,
							   taxa_are_rows=taxa_are_rows,
							   sampleda=sampleda,
							   tree=tree,
							   type=type,...)

	pcoares <- getpcoa.dist(tmpdist, distmethod=distmethod, sampleda=sampleda)
	return(pcoares)
}

#' @method getpcoa dist
#' @importFrom ape pcoa
#' @rdname getpcoa
#' @export
getpcoa.dist <- function(obj, distmethod, sampleda=NULL, ...){
	if (missing(distmethod)){
		if (!is.null(attr(obj, "distmethod"))){
			distmethod <- attr(obj, "distmethod")
		}else{
			distmethod <- NULL
		}
	}
	pcoares <- pcoa(obj, ...)
	attr(pcoares, "distmethod") <- distmethod
	res <- new("pcasample", 
			   pca=pcoares,
			   sampleda=sampleda)
	return(res)
}

#' @method getpcoa phyloseq
#' @rdname getpcoa
#' @export
getpcoa.phyloseq <- function(obj,distmethod="euclidean",...){
	sampleda <- checksample(obj)
	tmpdist <- getdist.phyloseq(obj, distmethod=distmethod,...)
	pcoares <- getpcoa.dist(tmpdist, distmethod=distmethod, sampleda=sampleda)
	return(pcoares)
}

#' @method getcoord pcoa
#' @rdname getcoord
#' @export
getcoord.pcoa <- function(obj, pc){
	coord <- obj$vector[,pc]
	
}



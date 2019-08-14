#' @title calculate distance
#' @param obj phyloseq, phyloseq class.
#' @param data data.frame, nrow sample * ncol feature. 
#' @param method character, default is hellinger, see alse \code{\link[vegan]{decostand}} 
#' @param distmethod character, default is "euclidean", see also \code{\link[phyloseq]{distanceMethodList}}
#' @param taxa_are_rows logical, default is FALSE.
#' @param sampleda data.frame, nrow sample * ncol factors.
#' @param tree phylo, see also \code{\link[ape]{phylo}}.
#' @param ..., additional parameters.
#' @return distance class contianed distmethod and originalD attr
#' @export
#' @examples
#' data(test_otu_data)
#' distclass <- getdist(test_otu_data)
#' hcsample <- getclust(distclass)
getdist <- function(obj,...){
	UseMethod("getdist")
}

#' @method getdist default
#' @rdname getdist
#' @importFrom vegan decostand
#' @importFrom phyloseq otu_table 
#' @export
getdist.default <- function(data, 
							distmethod="euclidean",
						    taxa_are_rows=FALSE,	
							sampleda=NULL,
							tree=NULL,
							method="hellinger",
							...){
	objphyloseq <- new("phyloseq",
					   otu_table=otu_table(data, 
									  taxa_are_rows=taxa_are_rows),
					   sam_data=sampleda,
					   phy_tree=tree)
	return(getdist.phyloseq(objphyloseq, 
							distmethod=distmethod, 
							method=method,
							type="sample", ...))
	
}

#' @method getdist phyloseq
#' @importFrom phyloseq distance taxa_are_rows phy_tree
#' @seealso \code{\link[phyloseq]{distance}}
#' @rdname getdist
#' @export
getdist.phyloseq <- function(obj, distmethod="euclidean", type="sample", method="hellinger",...){
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
	disres <- distance(obj, method=distmethod, type=type, ...)
	attr(disres, "distmethod") <- distmethod
	attr(disres, "originalD") <- data.frame(obj@otu_table, check.names=FALSE)
	return(disres)
}

#' @title performs principal coordinate analysis (PCoA)
#' @param data data.frame, numeric data.frame nrow sample * ncol features.
#' @param obj phyloseq, the phyloseq class or dist class.
#' @param sampleda data.frame, nrow sample * ncol factor, default is NULL.
#' @param distmethod character, the method of distance, see also \code{\link[phyloseq]{distance}}
#' @param taxa_are_rows logical, if feature of data is column, it should be set FALSE.
#' @param tree phylo, the phylo class, default is NULL, when use unifrac method, it should be
#' required.
#' @param method character, the standardization method for community ecologists, default is hellinger,
#' if the data has be normlized, it shound be set NULL.
#' @param ..., additional parameter, see \code{\link[MicrobiotaProcess]{getdist}}.
#' @return pcasample object, contained prcomp or pcoa and sampleda (data.frame).
#' @author Shuangbin Xu
#' @export
#' @examples
#' library(phyloseq)
#' data(GlobalPatterns)
#' subGlobal <- subset_samples(GlobalPatterns, SampleType %in% c("Feces", "Mock", "Ocean", "Skin"))
#' pcoares <- getpcoa(subGlobal, distmethod="euclidean", method="hellinger")
#' # pcoaplot <- ggordpoint(pcoares, biplot=FALSE,
#' #                        speciesannot=FALSE,
#' #                        factorNames=c("SampleType"), 
#' #                        ellipse=FALSE)
getpcoa <- function(obj, ...){
	UseMethod("getpcoa")
}

#' @method getpcoa default
#' @rdname getpcoa
#' @export
getpcoa.default <- function(data, 
							distmethod="euclidean", 
							taxa_are_rows=FALSE,
							sampleda=NULL, 
							tree=NULL,
							#type="sample",
							method="hellinger",
							...){
	tmpdist <- getdist.default(data=data, 
							   distmethod=distmethod,
							   taxa_are_rows=taxa_are_rows,
							   sampleda=sampleda,
							   tree=tree,
							   #type=type,
							   method=method,
							   ...)
	data <- attr(tmpdist, "originalD")
	pcoares <- getpcoa.dist(tmpdist, 
							distmethod=distmethod, 
							data=data,
							sampleda=sampleda)
	return(pcoares)
}

#' @method getpcoa dist
#' @importFrom ape pcoa
#' @importFrom stats cov
#' @rdname getpcoa
#' @export
getpcoa.dist <- function(obj, distmethod, data=NULL, sampleda=NULL, method="hellinger", ...){
	if (missing(distmethod)){
		if (!is.null(attr(obj, "distmethod"))){
			distmethod <- attr(obj, "distmethod")
		}else{
			distmethod <- NULL
		}
	}
	pcoares <- pcoa(obj, ...)
	attr(pcoares, "distmethod") <- distmethod
	if (!is.null(data)){
		n <- nrow(data)
		points.stand <- scale(pcoares$vectors)
		S <- cov(data, points.stand)
		tmpEig <- pcoares$values$Eigenvalues
		tmpEig <- tmpEig[seq_len(dim(S)[2])]
		diagtmp <- diag((tmpEig/(n-1))^(-0.5))
		U <- S %*% diagtmp
		colnames(U) <- colnames(pcoares$vectors)
		attr(pcoares, "varcorr") <- U
	}
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
	otuda <- attr(tmpdist, "originalD")
	pcoares <- getpcoa.dist(tmpdist, 
							distmethod=distmethod, 
							data=otuda, 
							sampleda=sampleda)
	return(pcoares)
}

#' @method getcoord pcoa
#' @rdname getcoord
#' @export
getcoord.pcoa <- function(obj, pc){
	coord <- obj$vector[,pc]
	vp <- obj$values$Relative_eig
	vp <- vp*100/sum(vp)
	tmpvp1 <- round(vp[pc[1]], 2)
	tmpvp2 <- round(vp[pc[2]], 2)
	xlab_text <- paste0("PCoA", pc[1], "(", tmpvp1, "%)")
	ylab_text <- paste0("PCoA", pc[2], "(", tmpvp2, "%)")
	title_text <- paste0("PCoA - PCoA",pc[1], " VS PCoA",pc[2], " (", attr(obj, "distmethod"), ")")
	ordplotclass <- new("ordplotClass",
						coord=coord,
						xlab=xlab_text,
						ylab=ylab_text,
						title=title_text
						)
	return(ordplotclass)
}

#' @method getvarct pcoa
#' @rdname getvarct
#' @export
getvarct.pcoa <- function(obj){
	if (is.null(attr(obj,"varcorr"))){
		stop("The pcoa class have not `varcorr` attr")
	}
	else{
		varcorr <- attr(obj, "varcorr")
		varcorr2 <- varcorr^2
		componentcos <- apply(varcorr2, 2, sum)
		varcontrib <- data.frame(t(apply(varcorr2,
										 1, 
										 contribution,
										 componentcos)),
								 check.names=FALSE)

		res <- list(VarContribution=varcontrib,
					 VarCoordinates=varcorr)
		attr(res, "class") <- "VarContrib"
		return(res)
	}
}

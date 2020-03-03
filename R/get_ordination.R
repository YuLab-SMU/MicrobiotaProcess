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

#' @method get_dist default
#' @rdname get_dist
#' @importFrom vegan decostand
#' @importFrom phyloseq otu_table 
#' @export
get_dist.default <- function(obj, 
    distmethod="euclidean",
    taxa_are_rows=FALSE,	
    sampleda=NULL,
    tree=NULL,
    method="hellinger",
    ...){
    objphyloseq <- new("phyloseq",
                       otu_table=otu_table(obj, 
                       taxa_are_rows=taxa_are_rows),
                       sam_data=sampleda,
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

#' @title performs principal coordinate analysis (PCoA)
#' @param data data.frame, numeric data.frame nrow sample * ncol features.
#' @param obj phyloseq, the phyloseq class or dist class.
#' @param sampleda data.frame, nrow sample * ncol factor, default is NULL.
#' @param distmethod character, the method of distance, 
#' see also \code{\link[phyloseq]{distance}}
#' @param taxa_are_rows logical, if feature of data is column, 
#' it should be set FALSE.
#' @param tree phylo, the phylo class, default is NULL, 
#' when use unifrac method, it should be required.
#' @param method character, the standardization method for 
#' community ecologists, default is hellinger, if the data 
#' has be normlized, it shound be set NULL.
#' @param ..., additional parameter, see also
#' \code{\link[MicrobiotaProcess]{get_dist}}.
#' @return pcasample object, contained prcomp or 
#' pcoa and sampleda (data.frame).
#' @author Shuangbin Xu
#' @export
#' @examples
#' #don't run in examples
#' #library(phyloseq)
#' #data(GlobalPatterns)
#' #subGlobal <- subset_samples(GlobalPatterns, 
#' #              SampleType %in% c("Feces", "Mock", "Ocean", "Skin"))
#' #pcoares <- get_pcoa(subGlobal, 
#' #                   distmethod="euclidean",
#' #                   method="hellinger")
#' # pcoaplot <- ggordpoint(pcoares, biplot=FALSE,
#' #                        speciesannot=FALSE,
#' #                        factorNames=c("SampleType"), 
#' #                        ellipse=FALSE)
get_pcoa <- function(obj, ...){
    UseMethod("get_pcoa")
}

#' @method get_pcoa default
#' @rdname get_pcoa
#' @export
get_pcoa.default <- function(obj, 
    distmethod="euclidean", 
    taxa_are_rows=FALSE,
    sampleda=NULL, 
    tree=NULL,
    #type="sample",
    method="hellinger",
    ...){
    tmpdist <- get_dist.default(obj, 
                                distmethod=distmethod,
                                taxa_are_rows=taxa_are_rows,
                                sampleda=sampleda,
                                tree=tree,
                                #type=type,
                                method=method,
                                ...)
    data <- attr(tmpdist, "originalD")
    pcoares <- get_pcoa.dist(tmpdist, 
                             distmethod=distmethod, 
                             data=data,
                             sampleda=sampleda)
    return(pcoares)
}

#' @method get_pcoa dist
#' @importFrom ape pcoa
#' @importFrom stats cov
#' @rdname get_pcoa
#' @export
get_pcoa.dist <- function(obj, distmethod, data=NULL, sampleda=NULL, method="hellinger", ...){
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

#' @method get_pcoa phyloseq
#' @rdname get_pcoa
#' @export
get_pcoa.phyloseq <- function(obj,distmethod="euclidean",...){
    sampleda <- checksample(obj)
    tmpdist <- get_dist.phyloseq(obj, distmethod=distmethod,...)
    otuda <- attr(tmpdist, "originalD")
    pcoares <- get_pcoa.dist(tmpdist, 
                             distmethod=distmethod, 
                             data=otuda, 
                             sampleda=sampleda)
    return(pcoares)
}

#' @method get_coord pcoa
#' @rdname get_coord
#' @export
get_coord.pcoa <- function(obj, pc){
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

#' @method get_varct pcoa
#' @rdname get_varct
#' @export
get_varct.pcoa <- function(obj,...){
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

#' @title Performs a principal components analysis
#' @param obj phyloseq, phyloseq class or data.frame
#' shape of data.frame is nrow sample * ncol feature.
#' @param sampleda data.frame, nrow sample * ncol factors.
#' @param method character, the standardization methods for
#' community ecologists. see \code{\link[vegan]{decostand}}.
#' @param ... additional parameters, see\code{\link[stats]{prcomp}}.
#' @return pcasample class, contained prcomp class and sample information.
#' @export
#' @examples
#' # don't run in examples
#' #library(phyloseq)
#' #data(GlobalPatterns)
#' #subGlobal <- subset_samples(GlobalPatterns, 
#' #         SampleType %in% c("Feces", "Mock", "Ocean", "Skin"))
#' #pcares <- get_pca(subGlobal, method="hellinger")
#' #pcaplot <- ggordpoint(pcares, biplot=TRUE, 
#' #                      speciesannot=TRUE,
#' #                      factorNames=c("SampleType"), ellipse=TRUE)
get_pca <- function(obj,...){
    UseMethod("get_pca")
}

#' @method get_pca data.frame
#' @importFrom vegan decostand
#' @importFrom stats prcomp
#' @rdname get_pca
#' @export
get_pca.data.frame <- function(obj,
    sampleda=NULL,
    method="hellinger",
    ...){
    if (!is.null(method)){
    	obj <- decostand(obj, method=method)
    }
    pca <- prcomp(obj, ...)
    #varcontrib <- get_varct(pca) 
    pca <- new("pcasample",
    		   pca=pca,
    		   #varcontrib=varcontrib,
    		   sampleda=sampleda)
    return(pca)
}

#' @method get_pca phyloseq
#' @rdname get_pca
#' @export
get_pca.phyloseq <- function(obj, method="hellinger", ...){
    otuda <- checkotu(obj)
    sampleda <- checksample(obj)
    pca <- get_pca.data.frame(otuda, sampleda=sampleda, method=method, ...)
    return(pca)
}

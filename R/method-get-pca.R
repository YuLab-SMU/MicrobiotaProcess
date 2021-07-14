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

#' Principal Components Analysis with MPSE or tbl_mpse object
#' @rdname mp_cal_pca-methods
#' @param .data MPSE or tbl_mpse object
#' @param .abundance the name of abundance to be calculated.
#' @param .dim integer The number of dimensions to be returned, default is 3.
#' @param action character "add" joins the pca result to the object, "only" return
#' a non-redundant tibble with the pca result. "get" return 'pcasample' object can
#' be visualized with 'ggordpoint'. 
#' @param ... additional parameters see also 'prcomp'
#' @return update object or tbl according to the action.
#' @export
setGeneric("mp_cal_pca", function(.data, .abundance, .dim=3, action="add", ...)standardGeneric("mp_cal_pca"))

#' @rdname mp_cal_pca-methods
#' @aliases mp_cal_pca,MPSE
#' @exportMethod mp_cal_pca
setMethod("mp_cal_pca", signature(.data="MPSE"), function(.data, .abundance, .dim=3, action="add", ...){

    action %<>% match.arg(c("add", "only", "get"))
    
    .abundance <- rlang::enquo(.abundance)

    x <- .data %>% mp_extract_abundance(.abundance=!!.abundance, byRow=FALSE)

    pca <- prcomp(x, ...)
    

    if (action=="get"){
        sampleda <- mp_extract_sample(.data) %>% 
                    tibble::column_to_rownames(var="Sample")
        res <- new("pcasample", pca=pca, sampleda=sampleda)
        return(res)
    }
    da <- .data %>%  
          mp_extract_sample() %>%
          dplyr::left_join(
                 pca$x[, seq_len(.dim)] %>% 
                 as_tibble(rownames="Sample"),
                 by="Sample"
          ) 
    if (action=="only"){
        da %<>%
             add_internals_attr(object=pca, name="PCA") 
        return(da)
    }else if (action=="add"){
        .data@colData <- da %>% 
                         tibble::column_to_rownames(var="Sample") %>% 
                         S4Vectors::DataFrame() 
        .data %<>% add_internals_attr(object=pca, name="PCA")
        return(.data)    
    }
})

.internal_cal_pca <- function(.data, .abundance, .dim=3, action="add", ...){
    
    action %<>% match.arg(c("add", "only", "get"))

    .abundance <- rlang::enquo(.abundance)

    x <- .data %>% mp_extract_abundance(.abundance=!!.abundance, byRow=FALSE)
    
    pca <- prcomp(x, ...)

    if (action=="get"){
        sampleda <- .data %>% 
                    mp_extract_sample() %>%
                    tibble::column_to_rownames(var="Sample")
        res <- new("pcasample", pca=pca, sampleda=sampleda)
        return(res)
    }else if (action=="only"){
        da <- .data %>%
              mp_extract_sample() %>%
              dplyr::left_join(
                  pca$x[,seq_len(.dim)] %>%
                  as_tibble(rownames="Sample"),
                 by="Sample"
              ) %>%
              add_internals_attr(object=pca, name="PCA")
        return(da)
    }else if (action=="add"){
        .data %<>% 
            dplyr::left_join(
                pca$x[,seq_len(.dim)] %>%
                as_tibble(rownames="Sample"),
                by="Sample"
            ) %>%
            add_internals_attr(object=pca, name="PCA")

        return(.data)
    }
}

#' @rdname mp_cal_pca-methods
#' @aliases mp_cal_pca,tbl_mpse
#' @exportMethod mp_cal_pca
setMethod("mp_cal_pca", signature(.data="tbl_mpse"), .internal_cal_pca)

#' @rdname mp_cal_pca-methods
#' @aliases mp_cal_pca,grouped_df_mpse
#' @exportMethod mp_cal_pca
setMethod("mp_cal_pca", signature(.data="grouped_df_mpse"), .internal_cal_pca)

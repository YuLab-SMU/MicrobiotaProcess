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
#' \dontrun{
#' library(phyloseq)
#' data(GlobalPatterns)
#' subGlobal <- subset_samples(GlobalPatterns, 
#'          SampleType %in% c("Feces", "Mock", "Ocean", "Skin"))
#' pcares <- get_pca(subGlobal, method="hellinger")
#' pcaplot <- ggordpoint(pcares, biplot=TRUE, 
#'                       speciesannot=TRUE,
#'                       factorNames=c("SampleType"), ellipse=TRUE)
#' }
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
#' @author Shuangbin Xu
#' @examples
#' data(mouse.time.mpse)
#' library(ggplot2)
#' mpse <- mouse.time.mpse %>% 
#'           mp_decostand(.abundance=Abundance) %>% 
#'           mp_cal_pca(.abundance=hellinger)
#' # action = "only" to extract the non-redundant tibble to visualize
#' tbl <- mouse.time.mpse %>%
#'           mp_decostand(.abundance=Abundance) %>%
#'           mp_cal_pca(.abundance=hellinger, action="only")
#' tbl
#' x <- names(tbl)[grepl("PC1 ", names(tbl))] %>% as.symbol()
#' y <- names(tbl)[grepl("PC2 ", names(tbl))] %>% as.symbol()
#' ggplot(tbl) + 
#'  geom_point(aes(x=!!x, y=!!y, color=time))
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
        #sampleda <- mp_extract_sample(.data) %>% 
        #            tibble::column_to_rownames(var="Sample")
        #res <- new("pcasample", pca=pca, sampleda=sampleda)
        #return(res)
        return(pca)
    }
    dat <- pca %>% tidydr(display=c("sites", "features")) 

    da <- .data %>%
          mp_extract_sample() %>%
          dplyr::left_join(
                 dat %>%
                 select(seq_len(.dim+1)),
                 by=c("Sample"="sites")
          )
    if (action=="only"){
        da %<>%
             add_attr(dat %>% attr("features_tb"), name="features_tb") %>%
             add_internal_attr(object=pca, name="PCA") 
        return(da)
    }else if (action=="add"){
        .data@colData <- da %>% 
                         tibble::column_to_rownames(var="Sample") %>% 
                         S4Vectors::DataFrame(check.names=FALSE) 
        .data %<>% add_internal_attr(object=pca, name="PCA")
        return(.data)    
    }
})

.internal_cal_pca <- function(.data, .abundance, .dim=3, action="add", ...){
    
    action %<>% match.arg(c("add", "only", "get"))

    .abundance <- rlang::enquo(.abundance)

    x <- .data %>% mp_extract_abundance(.abundance=!!.abundance, byRow=FALSE)
    
    pca <- prcomp(x, ...)

    if (action=="get"){
        #sampleda <- .data %>% 
        #            mp_extract_sample() %>%
        #            tibble::column_to_rownames(var="Sample")
        #res <- new("pcasample", pca=pca, sampleda=sampleda)
        return(pca)
    }
    
    dat <- pca %>% 
           tidydr(display=c("sites", "features"))

    if (action=="only"){
        da <- .data %>%
              mp_extract_sample() %>%
              dplyr::left_join(
                 dat[, seq_len(.dim+1)],
                 by=c("Sample"="sites")
              ) %>%
              add_attr(attribute=dat %>% attr("features_tb"), 
                       name="features_tb") %>%
              add_internal_attr(object=pca, name="PCA")
        return(da)
    }else if (action=="add"){
        .data %<>% 
            dplyr::left_join(
                dat[,seq_len(.dim+1)],
                by=c("Sample"="sites")
            ) %>%
            add_internal_attr(object=pca, name="PCA")

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


#' Detrended Correspondence Analysis with MPSE or tbl_mpse object
#' @rdname mp_cal_dca-methods
#' @param .data MPSE or tbl_mpse object
#' @param .abundance the name of abundance to be calculated.
#' @param .dim integer The number of dimensions to be returned, default is 3.
#' @param action character "add" joins the 'decorana' result to the object, "only" return
#' a non-redundant tibble with the 'decorana' result. "get" return 'decorana' object can
#' be processed with related vegan function.
#' @param origin logical Use true origin even in detrended correspondence analysis.
#' default is TRUE.
#' @param ... additional parameters see also 'vegan::decorana'
#' @return update object or tbl according to the action.
#' @export
setGeneric("mp_cal_dca", function(.data, .abundance, .dim=3, action="add", origin=TRUE, ...)standardGeneric("mp_cal_dca"))

#' @rdname mp_cal_dca-methods
#' @aliases mp_cal_dca,MPSE
#' @exportMethod mp_cal_dca
setMethod("mp_cal_dca", signature(.data="MPSE"), function(.data, .abundance, .dim=3, action="add", origin=TRUE, ...){

    action %<>% match.arg(c("add", "only", "get"))

    .abundance <- rlang::enquo(.abundance)

    x <- .data %>% mp_extract_abundance(.abundance=!!.abundance, byRow=FALSE)

    dca <- vegan::decorana(x, ...)

    if (action=="get"){
        #sampleda <- mp_extract_sample(.data) %>%
        #            tibble::column_to_rownames(var="Sample")
        return(dca)
    }
    
    dat <- dca %>% 
           tidydr(display="features", origin=origin) %>%
           select(seq_len(.dim + 1))

    da <- .data %>%
          mp_extract_sample() %>%
          dplyr::left_join(
                 dat,
                 by=c("Sample"="sites")
          )

    if (action=="only"){
        da %<>%
             add_attr(dat %>% attr("features_tb"), name="features_tb") %>%
             add_internal_attr(object=dca, name="DCA")
        return(da)
    }else if (action=="add"){
        .data@colData <- da %>%
                         tibble::column_to_rownames(var="Sample") %>%
                         S4Vectors::DataFrame(check.names=FALSE)
        .data %<>% add_internal_attr(object=dca, name="DCA")
        return(.data)
    }            
            
})

.internal_cal_dca <- function(.data, .abundance, .dim=3, action="add", origin=TRUE, ...){

    action %<>% match.arg(c("add", "only", "get"))

    .abundance <- rlang::enquo(.abundance)

    x <- .data %>% mp_extract_abundance(.abundance=!!.abundance, byRow=FALSE)

    dca <- vegan::decorana(x, ...)

    if (action=="get"){
        #sampleda <- .data %>%
        #            mp_extract_sample() %>%
        #            tibble::column_to_rownames(var="Sample")
        #res <- new("pcasample", pca=pca, sampleda=sampleda)
        return(dca)
    }

    dat <- dca %>%
           tidydr(display="features", origin=origin) %>%
           select(seq_len(.dim + 1))

    if (action=="only"){
        da <- .data %>%
              mp_extract_sample() %>%
              dplyr::left_join(
                  dat,
                  by=c("Sample"="sites")
              ) %>%
              add_internal_attr(object=dca, name="DCA")
        return(da)
    }else if (action=="add"){
        .data %<>%
            dplyr::left_join(
                dat,
                by=c("Sample"="sites")
            ) %>%
            add_internal_attr(object=dca, name="DCA")

        return(.data)
    }
}

#' @rdname mp_cal_dca-methods
#' @aliases mp_cal_dca,tbl_mpse
#' @exportMethod mp_cal_dca
setMethod("mp_cal_dca", signature(.data="tbl_mpse"), .internal_cal_dca)

#' @rdname mp_cal_dca-methods
#' @aliases mp_cal_dca,grouped_df_mpse
#' @exportMethod mp_cal_dca
setMethod("mp_cal_dca", signature(.data="grouped_df_mpse"), .internal_cal_dca)

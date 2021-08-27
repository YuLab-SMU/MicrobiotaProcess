##' generate the result of rare curve.
##'
##' This function is designed to calculate the rare curve result of otu table
##' the result can be visualized by `ggrarecurve`.
##'
##' @title obtain the result of rare curve
##' @param obj phyloseq class or data.frame
##' shape of data.frame (nrow sample * ncol feature)
##' @param sampleda data.frame, (nrow sample * ncol factor)
##' @param chunks integer, the number of subsample in a sample, 
##' default is 400.
##' @param factorLevels list, the levels of the factors, default is NULL,
##' if you want to order the levels of factor, you can set this.
##' @param ..., additional parameters.
##' @return rarecurve class, which can be visualized by ggrarecurve
##' @author Shuangbin Xu
##' @export
##' @examples
##' \dontrun{
##'     data(test_otu_data)
##'     set.seed(1024)
##'     res <- get_rarecurve(test_otu_data, chunks=200)
##'     p <- ggrarecurve(obj=res, 
##'                      indexNames=c("Observe","Chao1","ACE"),
##'                      shadow=FALSE,
##'                      factorNames="group")
##' }
setGeneric("get_rarecurve", function(obj, ...)standardGeneric("get_rarecurve"))

#' @aliases get_rarecurve,data.frame
#' @rdname get_rarecurve
#' @export
setMethod("get_rarecurve", "data.frame", function(obj, sampleda, factorLevels=NULL, chunks=400){
    res <- stat_rare(data=obj, sampleda=sampleda, 
                     chunks=chunks, factorLevels=factorLevels, 
                     plotda=TRUE)
    res <- structure(list(data=res), class="rarecurve")
    return(res)

})

#' @aliases get_rarecurve,phyloseq
#' @rdname get_rarecurve
#' @export
setMethod("get_rarecurve", "phyloseq", function(obj, ...){
    otuda <- checkotu(obj)
    sampleda <- checksample(obj)
    res <- get_rarecurve(obj=otuda,
                         sampleda=sampleda,
                         ...)
    return(res)
})


#' Calculating the different alpha diversities index with different depth
#' @rdname mp_cal_rarecurve-methods
#' @param .data MPSE or tbl_mpse object
#' @param .abundance the name of otu abundance to be calculated.
#' @param action character it has three options, "add" joins the new
#' information to the input tbl (default), "only" return a non-redundant tibble 
#' with the just new information, ang 'get' return a 'rarecurve' object.
#' @param chunks numeric the split number of each sample to 
#' calculate alpha diversity, default is 400. eg. A sample has total 40000 reads,
#' if chunks is 400, it will be split to 100 sub-samples (100, 200, 300,..., 40000),
#' then alpha diversity index was calculated based on the sub-samples.
#' @param seed a random seed to make the result reproducible, default is 123.
#' @param force logical whether calculate rarecurve forcibly when the
#' '.abundance' is not be rarefied, default is FALSE
#' @param ... additional parameters.
#' @return update rarecurce calss
#' @seealso [mp_plot_rarecurve()]
#' @author Shuangbin Xu
#' @export
#' @examples
#' data(mouse.time.mpse)
#' mouse.time.mpse %>% 
#' mp_rrarefy() -> mpse
#' mpse
#' # bigger chunks means more robust, but it will become slower.
#' mpse %<>% mp_cal_rarecurve(.abundance=RareAbundance, chunks=100, action="add")
#' mpse
#' p1 <- mpse %>% mp_plot_rarecurve(.rare=RareAbundanceRarecurve, .alpha="Observe")
setGeneric("mp_cal_rarecurve", function(.data, .abundance=NULL, action="add", chunks=400, seed=123, force=FALSE, ...)standardGeneric("mp_cal_rarecurve"))

#' @rdname mp_cal_rarecurve-methods
#' @aliases mp_cal_rarecurve,MPSE
#' @exportMethod mp_cal_rarecurve
setMethod("mp_cal_rarecurve", signature(.data="MPSE"), function(.data, .abundance=NULL, action="add", chunks=400, seed=123, force=FALSE, ...){

    action %<>% match.arg(c("add", "get", "only"))
    .abundance <- rlang::enquo(.abundance)

    if (rlang::quo_is_null(.abundance)){
        .abundance <- as.symbol("RareAbundance")
    }

    if (!valid_rare(.data, .abundance=.abundance) && !force){
        glue::glue("The rarefied abundance of species might not be provided. Rarefaction of all
                    observations is performed automatically. If you still want to calculate the
                    alpha index with the '.abundance', you can set 'force=TRUE'. ")
        .data <- mp_rrarefy(.data=.data, ...)
        .abundance <- as.symbol("RareAbundance")
    }

    clnm <- paste0(rlang::as_name(.abundance), "Rarecurve")

    if (clnm %in% colnames(.data@colData)){
        message("The rarecurve has been found in the colData of the MPSE object ")
        return (.data)
    }

    da <- SummarizedExperiment::assays(.data)@listData[[rlang::as_name(.abundance)]] %>% as.data.frame(check.names=FALSE)
    
    sampleda <- .data %>%
                mp_extract_sample()

    dat <- .internal_apply_cal_rarecurve(da=da, sampleda=sampleda, chunks=chunks, seed=seed)

    if (action=="get"){
        dat <- structure(list(data=dat), class="rarecurve")
        return(dat)
    }
    
    nestnm <- colnames(dat)[!colnames(dat) %in% "sample"]
    params <- list(.data=dat, rarecurve=nestnm)
    dat <- do.call("nest", params) %>%
           dplyr::rename(!!clnm:="rarecurve")
    
    if (action=="only"){
        return(dat)
    }else if (action=="add"){
        .data@colData <- .data@colData %>%
            as_tibble(rownames="Sample") %>%
            left_join(dat, by=c("Sample"="sample"), suffix=c("", ".y")) %>% 
            column_to_rownames(var="Sample") %>%
            S4Vectors::DataFrame(check.names=FALSE)
        return(.data)
    }
    
})


.internal_mp_cal_rarecurve <- function(.data, .abundance=NULL, action="add", chunks=400, seed=123, force=FALSE, ...){
    
    action %<>% match.arg(c("add", "get", "only"))
    .abundance <- rlang::enquo(.abundance)

    if (rlang::quo_is_null(.abundance)){
        .abundance <- as.symbol("RareAbundance")
    }    

    if (!valid_rare(.data, .abundance=.abundance) && !force){
        glue::glue("The rarefied abundance of species might not be provided. Rarefaction of all
                    observations is performed automatically. If you still want to calculate the
                    alpha index with the '.abundance', you can set 'force=TRUE'. ")
        .data <- mp_rrarefy(.data=.data, ...)
        .abundance <- as.symbol("RareAbundance")
    }

    clnm <- paste0(rlang::as_name(.abundance), "Rarecurve")

    if (paste0(rlang::as_name(.abundance), "Rarecurve") %in% colnames(.data) &&
        .data %>%
        dplyr::pull("rarecurve") %>%
        class == "list" ||
        all(c("readsNums", "Alpha") %in% colnames(.data))){
        message("The rarecurve has been found in the object !")
        return (.data)
    }    

    da <- .data %>% 
           dplyr::ungroup() %>% 
           select(c("OTU", "Sample", rlang::as_name(.abundance))) %>%
           tidyr::pivot_wider(id_cols="OTU", names_from="Sample", values_from=rlang::as_name(.abundance)) %>%
           tibble::column_to_rownames(var="OTU")

    samplevar <- .data %>% attr("samplevar")
    sampleda <- .data %>% dplyr::ungroup() %>% select(samplevar) %>% distinct()

    dat <- .internal_apply_cal_rarecurve(da=da, sampleda=sampleda, chunks=chunks, seed=seed)
    if (action=="get"){
        dat <- structure(list(data=dat), class="rarecurve")
        return(dat)
    }
    
    nestnm <- colnames(dat)[!colnames(dat) %in% "sample"]
    params <- list(.data=dat, rarecurve=nestnm)
    dat <- do.call("nest", params) %>%
           dplyr::rename(!!clnm:="rarecurve")

    if (action=="only"){
        return(dat)
    }else if (action=="add"){
        dat <- .data %>% 
               dplyr::left_join(dat, by=c("Sample"="sample"), suffix=c("", ".y"))
        return(dat)
    }
}

.internal_apply_cal_rarecurve <- function(da, sampleda, chunks, seed){
    dat <- apply(da, 2, samplealpha, chunks=chunks, seed=seed) %>% 
           dplyr::bind_rows(.id="sample") %>%
           tidyr::pivot_longer(cols=!c("sample", "readsNums"), names_to="Alpha")
    sampleda <- sampleda[, !vapply(sampleda, is.list, logical(1))] 
    sampleda <- sampleda[, !colnames(sampleda) %in% colnames(dat)]
    if (ncol(sampleda) > 1){
        dat %<>%
            left_join(sampleda, by=c("sample"="Sample"), suffix=c("", ".y"))
    }
    #dat <- structure(list(data=dat), class="rarecurve")
    return(dat)
}

#' @rdname mp_cal_rarecurve-methods
#' @aliases mp_cal_rarecurve,tbl_mpse
#' @exportMethod mp_cal_rarecurve
setMethod("mp_cal_rarecurve", signature(.data="tbl_mpse"), .internal_mp_cal_rarecurve)

#' @rdname mp_cal_rarecurve-methods
#' @aliases mp_cal_rarecurve,grouped_df_mpse
#' @exportMethod mp_cal_rarecurve
setMethod("mp_cal_rarecurve", signature(.data="grouped_df_mpse"), .internal_mp_cal_rarecurve)

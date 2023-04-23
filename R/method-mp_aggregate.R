#' aggregate the assays with the specific group of sample and fun.
#' @rdname mp_aggregate-methods
#' @param .data MPSE object, required
#' @param .abundance the column names of abundance, default is Abundance.
#' @param .group the column names of sample meta-data, required
#' @param fun a function to compute the summary statistics, default is sum.
#' @param keep_colData logical whether to keep the sample meta-data with \code{.group} as row names,
#' default is TRUE.
#' @param ... additional parameters, see also \code{\link[stats]{aggregate}}.
#' @return a new object with .group as column names in assays
#' @export
#' @examples
#' \dontrun{
#' data(mouse.time.mpse)
#' newmpse <- mouse.time.mpse %>%
#'            mp_aggregate(.group = time)
#' newmpse
#' }
setGeneric("mp_aggregate", function(.data, .abundance, .group, fun=sum, keep_colData=TRUE, ...)standardGeneric("mp_aggregate"))

.internal_mp_aggregate <- function(.data, .abundance, .group, fun=sum, keep_colData=TRUE, ...){
    .abundance <- rlang::enquo(.abundance)
    .group <- rlang::enquo(.group)

    if (rlang::quo_is_missing(.abundance)){
        .abundance <- as.symbol("Abundance")
    }

    assayda <- .data %>% 
               mp_extract_assays(.abundance=!!.abundance, byRow=FALSE) %>%
               tibble::as_tibble(rownames="Sample")
    sampleda <- .data %>% mp_extract_sample()
    sampleda1 <- sampleda %>% select(!!rlang::sym("Sample"), !!.group)
    assayda %<>% 
        left_join(sampleda1, by="Sample") %>% 
        dplyr::select(-!!rlang::sym("Sample"))
    colData(.data) <- NULL
    fma <- as.formula(paste0(". ~", rlang::as_name(.group)))
    assayda <- stats::aggregate(fma, data=assayda, FUN=fun, ...)
    assayda %<>% 
             tibble::column_to_rownames(var=rlang::as_name(.group)) %>% 
             t()
    newda <- MPSE(assays=list(Abundance=assayda))
    taxatree(newda) <- taxatree(.data)
    otutree(newda) <- otutree(.data)
    refsequence(newda) <- refsequence(.data)
    SummarizedExperiment::rowData(newda) <- SummarizedExperiment::rowData(.data)

    if (keep_colData){
        sampleda %<>% dplyr::rename(.Old.Sample=!!rlang::sym("Sample"), Sample=!!.group) %>%
                      .internal_nest_tibble(columns="Sample")
        if (ncol(sampleda)>1){
            colData(newda) <- sampleda %>% 
                              tibble::column_to_rownames(var="Sample") %>% 
                              S4Vectors::DataFrame(check.names=FALSE)
        }
    }
    return (newda)
}

#' @rdname mp_aggregate-methods
#' @aliases mp_aggregate,MPSE
#' @exportMethod mp_aggregate
setMethod("mp_aggregate", signature(.data="MPSE"), .internal_mp_aggregate)

.internal_nest_tibble <- function(x, columns){
    nm <- colnames(x) 
    nm %<>% stats::setNames(nm)
    nm <- nm[!nm %in% columns]
    nm %<>% lapply(., function(x)x)
    x <- do.call(tidyr::nest, c(list(x), nm))
    x <- check_single_nrow_in_nest(x, columns)
    return(x)
}

check_single_nrow_in_nest <- function(da, columns){
    indnm <- da %>% 
             apply(., 2, function(x)all(lapply(x, function(i)nrow(unique(i))==1) %>% unlist())) 
    indnm <- indnm[indnm] %>% names()
    indnm <- indnm[!indnm %in% columns]
    if (length(indnm)>0){
        for( i in indnm){
            da %<>% dplyr::mutate(!!rlang::sym(i):=as.vector(unlist(lapply(!!rlang::sym(i), function(x)unique(x)))))
        }
    }
    return(da)
}



##' @title This Function Provideds Several Standardization Methods for Community Data
##' @docType methods
##' @name mp_decostand
##' @rdname mp_decostand-methods
##' @param .data MPSE or tbl_mpse object
##' @param .abundance the names of otu abundance to be applied standardization. 
##' @param method character the name of standardization method, it can one of 
##' 'total', 'max', 'frequency', 'normalize', 'range', 'rank', 'rrank', 'standardize'
##' 'pa', 'chi.square', 'hellinger' and 'log', see also \code{\link[vegan]{decostand}}
##' @param logbase numeric The logarithm base used in 'method=log', default is 2.
##' @param ... additional parameters, see also \code{\link[vegan]{decostand}}
##' @return update object
##' @author Shuangbin Xu
##' @export
##' @examples
##' data(mouse.time.mpse)
##' mouse.time.mpse %>% 
##' mp_decostand(.abundance=Abundance, method="hellinger")
setGeneric("mp_decostand", 
           function(.data, .abundance=NULL, method="hellinger", logbase=2, ...){
               standardGeneric("mp_decostand")}
)

##' @rdname mp_decostand-methods
##' @aliases mp_decostand,data.frame
##' @exportMethod mp_decostand
##' @source
##' mp_decostand for data.frame object is a wrapper method of vegan::decostand from the vegan package
##' @seealso
##' \link[vegan]{decostand}
setMethod("mp_decostand", signature(.data="data.frame"), 
          function(.data, method="hellinger", logbase=2, ...){
   x <- vegan::decostand(x=.data, method=method, logbase=logbase, ...)
   return(x)
})

##' @rdname mp_decostand-methods
##' @aliases mp_decostand,MPSE
##' @exportMethod mp_decostand
setMethod("mp_decostand", signature(.data="MPSE"),function(.data, .abundance=NULL, method="hellinger", logbase=2, ...){
    .abundance <- rlang::enquo(.abundance)

    if (method=="log"){
        newnm <- paste0(method, logbase)
    }else{
        newnm <- method
    }

    assaysvar <- .data %>% 
                  SummarizedExperiment::assayNames() 

    #if (newnm %in% assaysvar){
    #    message(paste("The ", newnm, " has exited in assays of the MPSE object"))
    #    return (.data)
    #}

    if (rlang::quo_is_null(.abundance)){
        .abundance <- as.symbol("Abundance")
    }
    
    xx <- SummarizedExperiment::assays(.data)@listData

    da <- .data %>%
          mp_extract_assays(.abundance=!!.abundance, byRow=FALSE)

    newda <- da %>% mp_decostand(method=method, logbase=logbase, ...) %>% t()

    if (method=="log"){
       newnm <- paste0(method, logbase)
    }else{
       newnm <- method
    }
    
    SummarizedExperiment::assays(.data)@listData <- c(xx, list(newda)) %>% 
                                                    setNames(c(assaysvar, newnm))
    return(.data)
})

# #' @rdname mp_decostand-methods
# #' @aliases mp_decostand,tbl_mpse
# #' @exportMethod mp_decostand

.internal_decostand <- function(.data, .abundance=NULL, method="hellinger", logbase=2, ...){
    .abundance = rlang::enquo(.abundance)
    

    if (method=="log"){
        newnm <- paste0(method, logbase)
    }else{
        newnm <- method
    }

    assaysvar <- .data %>% attr("assaysvar")

    if (rlang::quo_is_null(.abundance)){
        .abundance <- as.symbol("Abundance")
    }

    othernms <- colnames(.data)[!colnames(.data) %in% c("OTU", "Sample", assaysvar)]

    newda <- .data %>% 
             mp_extract_assays(.abundance=!!.abundance, byRow=FALSE) %>%
             mp_decostand(method=method, logbase=logbase, ...) %>%
             tibble::as_tibble(rownames="Sample") %>% 
             tidyr::pivot_longer(!as.symbol("Sample"), values_to=newnm, names_to="OTU")

    res <- .data %>% 
            dplyr::left_join(newda, by=c("OTU", "Sample"), suffix=c("", ".y")) %>%
            select(c("OTU", "Sample", assaysvar, newnm, othernms))

    res <- add_attr.tbl_mpse(x1 = res, x2 = .data)
    attr(res, "assaysvar") <- c(assaysvar, newnm)
    return(res)
}

##' @rdname mp_decostand-methods
##' @aliases mp_decostand,tbl_mpse
##' @exportMethod mp_decostand
setMethod("mp_decostand", signature(.data="tbl_mpse"), .internal_decostand)


#' @rdname mp_decostand-methods
#' @aliases mp_decostand,grouped_df_mpse
#' @exportMethod mp_decostand
setMethod("mp_decostand", signature(.data="grouped_df_mpse"), .internal_decostand)

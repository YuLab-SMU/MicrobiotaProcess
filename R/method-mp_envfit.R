#' Fits an Environmental Vector or Factor onto an Ordination With MPSE or tbl_mpse Object
#' @rdname mp_envfit-methods
#' @param .data MPSE or tbl_mpse object
#' @param .ord a name of ordination, option it is DCA, NMDS, RDA, CCA.
#' @param .env the names of columns of sample group or environment information.
#' @param .dim integer The number of dimensions to be returned, default is 3. 
#' @param action character "add" joins the envfit result to internal attributes of the object, 
#' "only" return a non-redundant tibble with the envfit result. "get" return 'envfit' object can
#' be analyzed using the related vegan funtion. 
#' @param permutations the number of permutations required, default is 999.
#' @param seed a random seed to make the analysis reproducible, default is 123. 
#' @param ... additional parameters see also 'vegan::envfit'
#' @return update object according action
#' @export
setGeneric("mp_envfit", function(.data, .ord, .env, .dim=3, action="add", permutations=999, seed=123, ...)standardGeneric("mp_envfit"))

.internal_cal_envfit <- function(.data, .ord, .env, .dim, action="add", permutations=999, seed=123, ...){

    .ord <- rlang::enquo(.ord) %>%
            rlang::as_name() %>%
            toupper()

    .env <- rlang::enquo(.env)

    .ord %<>% match.arg(c("NMDS", "RDA", "CCA", "DCA"))

    ordobj <- .data %>% 
              mp_extract_internal_attr(name=.ord)

    if (is.null(ordobj)){
        ordfun <- switch(.ord,
                         NMDS = "mp_cal_nmds",
                         RDA  = "mp_cal_rda",
                         CCA  = "mp_cal_cca",
                         DCA  = "mp_cal_dca")
        rlang::abort(paste0("The ", rlang::as_name(.ord), " is not present in the object, please run ", ordfun, "before performing mp_envfit."))
    }

    tmpX <- vegan::scores(ordobj, display="site")
    envda <- .data %>%
             mp_extract_sample() %>%
             arrange(match(rownames(tmpX), !!as.symbol("Sample"))) %>%
             tibble::column_to_rownames(var="Sample") %>%
             select(!!.env)

    res <- withr::with_seed(seed, vegan::envfit(ordobj, envda, permutations=permutations, choices=seq_len(.dim), ...))

    if (action=="get"){
        return(res)
    }else if (action=="only"){
        da <- .data %>% 
              mp_extract_sample() %>%
              add_attr(res %>% mp_fortify(), name=paste0(.ord, "_ENVFIT_tb"))
        return(da)
    }else if (action=="add"){
        attrnames <- paste0(.ord, "_ENVFIT")
        message("The result of mp_envfit has been saved to the internal attribute of the object !")
        message(paste0("It can be extracted using this-object %>% mp_extract_internal_attr(name='",attrnames, "')"))
        .data %<>%
            add_internal_attr(object=res, name=attrnames)
        return(.data)
    }
}

#' @rdname mp_envfit-methods
#' @aliases mp_envfit,MPSE
#' @exportMethod mp_envfit
setMethod("mp_envfit", signature(.data="MPSE"), .internal_cal_envfit)

#' @rdname mp_envfit-methods
#' @aliases mp_envfit,tbl_mpse
#' @exportMethod mp_envfit
setMethod("mp_envfit", signature(.data="tbl_mpse"), .internal_cal_envfit)

#' @rdname mp_envfit-methods
#' @aliases mp_envfit,grouped_df_mpse
#' @exportMethod mp_envfit
setMethod("mp_envfit", signature(.data="grouped_df_mpse"), .internal_cal_envfit)

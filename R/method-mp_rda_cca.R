#' [Partial] [Constrained] Correspondence Analysis with MPSE or tbl_mpse object
#' @rdname mp_cal_cca-methods
#' @param .data MPSE or tbl_mpse object
#' @param .abundance the name of abundance to be calculated.
#' @param .formula Model formula right hand side gives the constraining
#' variables, and conditioning variables can be given within a special 
#' function 'Condition' and keep left empty, such as ~ A + B or ~ A + Condition(B), 
#' default is NULL.
#' @param .dim integer The number of dimensions to be returned, default is 3.
#' @param action character "add" joins the cca result to the object, "only" return
#' a non-redundant tibble with the cca result. "get" return 'cca' object can
#' be analyzed using the related vegan funtion.
#' @param ... additional parameters see also 'cca' of vegan.
#' @return update object according action argument
#' @export
setGeneric("mp_cal_cca", function(.data, .abundance, .formula=NULL, .dim=3, action="add", ...)standardGeneric("mp_cal_cca"))

#' @rdname mp_cal_cca-methods
#' @aliases mp_cal_cca,MPSE
#' @exportMethod mp_cal_cca
setMethod("mp_cal_cca", signature(.data="MPSE"),function(.data, .abundance, .formula, .dim=3, action="add", ...){
    res <- .internal_cal_rda_cca1(.data=.data, 
                                  .abundance=!!rlang::enquo(.abundance), 
                                  .formula=.formula, 
                                  .dim=.dim,
                                  action=action,
                                  method="cca",
                                  ...
                                 )
    return(res)
})

.internal_cal_rda_cca1 <- function(.data, .abundance, .formula, .dim=3, action="add", method="cca", ...){
    
    cca_rda_method <- switch(method,
                            cca = vegan::cca,
                            rda = vegan::rda)

    action %<>% match.arg(c("add", "only", "get"))

    .abundance <- rlang::enquo(.abundance)

    x <- .data %>% mp_extract_abundance(.abundance=!!.abundance, byRow=FALSE)

    if (!is.null(.formula)){
        .formula <- paste0(c("x", .formula), collapse=" ") %>% as.formula()
        sampleda <- .data %>% 
                    mp_extract_sample() %>%
                    dplyr::arrange(match(rownames(x), !!as.symbol("Sample"))) %>%
                    tibble::column_to_rownames(var="Sample")
        ccares <- cca_rda_method(.formula, data=sampleda, ...)
    }else{
        ccares <- cca_rda_method(x, ...)
    }

    if (action == "get"){
        return(ccares)
    }

    dat <- vegan::scores(ccares, choices=seq_len(.dim))$sites %>%
           as_tibble(rownames="Sample")

    da <- .data %>% 
          mp_extract_sample() %>%
          dplyr::left_join(dat, by="Sample")

    if (action == "only"){
        da %<>%
            add_internal_attr(object=ccares,
                              name = switch(method, cca="CCA", rda="RDA"))
        return(da)
    }else if (action == "add"){
        .data@colData <- da %>% 
                         tibble::column_to_rownames(var="Sample") %>%
                         S4Vectors::DataFrame(check.names=FALSE)
        .data %<>% add_internal_attr(object=ccares,
                                     name = switch(method, cca="CCA", rda="RDA"))
        return(.data)
    }
}


.internal_cal_rda_cca2 <- function(.data, .abundance, .formula=NULL, .dim=3, action="add", method="cca", ...){
    
    cca_rda_method <- switch(method,
                            cca = vegan::cca,
                            rda = vegan::rda)

    action %<>% match.arg(c("add", "only", "get"))

    .abundance <- rlang::enquo(.abundance)

    x <- .data %>% mp_extract_abundance(.abundance=!!.abundance, byRow=FALSE)

    if (!is.null(.formula)){
        .formula <- paste0(c("x", .formula), collapse=" ") %>% as.formula()
        sampleda <- .data %>%
                    mp_extract_sample() %>%
                    dplyr::arrange(match(rownames(x), !!as.symbol("Sample"))) %>%
                    tibble::column_to_rownames(var="Sample")
        ccares <- cca_rda_method(.formula, data=sampleda, ...)
    }else{
        ccares <- cca_rda_method(x, ...)
    }

    dat <- vegan::scores(ccares, choices=seq_len(.dim))$sites %>%
           as_tibble(rownames="Sample")    

    if (action=="get"){
        return(ccares)
    }
    
    if (action=="only"){
        da <- .data %>%
              mp_extract_sample() %>%
              dplyr::left_join(
                  dat,
                  by="Sample"
              ) %>%
              add_internal_attr(object=ccares, 
                                name=switch(method, cca="CCA", rda="RDA"))
        return(da)
    }else if (action=="add"){
        .data %<>%
            dplyr::left_join(
                 dat,
                 by="Sample"
            ) %>%
            add_internal_attr(object=ccares, 
                              name=switch(method, cca="CCA", rda="RDA"))

        return(.data)
    }
}

#' @rdname mp_cal_cca-methods
#' @aliases mp_cal_cca,tbl_mpse
#' @exportMethod mp_cal_cca
setMethod("mp_cal_cca", signature(.data="tbl_mpse"), function(.data, .abundance, .formula=NULL, .dim=3, action="add", ...){
    res <- .internal_cal_rda_cca2(.data=.data,
                                  .abundance = !!rlang::enquo(.abundance),
                                  .formula = .formula,
                                  .dim = .dim,
                                  action = action,
                                  method = "cca",
                                  ...)
    return(res)    
})

#' @rdname mp_cal_cca-methods
#' @aliases mp_cal_cca,grouped_df_mpse
#' @exportMethod mp_cal_cca
setMethod("mp_cal_cca", signature(.data="grouped_df_mpse"), function(.data, .abundance, .formula=NULL, .dim=3, action="add", ...){
    res <- .internal_cal_rda_cca2(.data=.data,
                                  .abundance = !!rlang::enquo(.abundance),
                                  .formula = .formula,
                                  .dim = .dim,
                                  action = action,
                                  method = "cca",
                                  ...)
    return(res)
})


#' [Partial] [Constrained] Redundancy Analysis with MPSE or tbl_mpse object
#' @rdname mp_cal_rda-methods
#' @param .data MPSE or tbl_mpse object
#' @param .abundance the name of abundance to be calculated.
#' @param .formula Model formula right hand side gives the constraining
#' variables, and conditioning variables can be given within a special
#' function 'Condition' and keep left empty, such as ~ A + B or ~ A + Condition(B),
#' default is NULL.
#' @param .dim integer The number of dimensions to be returned, default is 3.
#' @param action character "add" joins the rda result to the object, "only" return
#' a non-redundant tibble with the rda result. "get" return 'rda' object can
#' be analyzed using the related vegan funtion.
#' @param ... additional parameters see also 'rda' of vegan.
#' @return update object according action argument
#' @export
setGeneric("mp_cal_rda", function(.data, .abundance, .formula=NULL, .dim=3, action="add", ...)standardGeneric("mp_cal_rda"))

#' @rdname mp_cal_rda-methods
#' @aliases mp_cal_rda,MPSE
#' @exportMethod mp_cal_rda
setMethod("mp_cal_rda", signature(.data="MPSE"),function(.data, .abundance, .formula=NULL, .dim=3, action="add", ...){
    res <- .internal_cal_rda_cca1(.data=.data,
                                  .abundance=!!rlang::enquo(.abundance),
                                  .formula=.formula,
                                  .dim=.dim,
                                  action=action,
                                  method="rda",
                                  ...
                                 )
    return(res)
})

#' @rdname mp_cal_rda-methods
#' @aliases mp_cal_rda,tbl_mpse
#' @exportMethod mp_cal_rda
setMethod("mp_cal_rda", signature(.data="tbl_mpse"), function(.data, .abundance, .formula=NULL, .dim=3, action="add", ...){
    res <- .internal_cal_rda_cca2(.data=.data,
                                  .abundance = !!rlang::enquo(.abundance),
                                  .formula = .formula,
                                  .dim = .dim,
                                  action = action,
                                  method = "rda",
                                  ...)
    return(res)
})

#' @rdname mp_cal_rda-methods
#' @aliases mp_cal_rda,grouped_df_mpse
#' @exportMethod mp_cal_rda
setMethod("mp_cal_rda", signature(.data="grouped_df_mpse"), function(.data, .abundance, .formula=NULL, .dim=3, action="add", ...){
    res <- .internal_cal_rda_cca2(.data=.data,
                                  .abundance = !!rlang::enquo(.abundance),
                                  .formula = .formula,
                                  .dim = .dim,
                                  action = action,
                                  method = "rda",
                                  ...)
    return(res)    
})




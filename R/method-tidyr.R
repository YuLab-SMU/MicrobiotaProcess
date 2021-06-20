#' @method unnest tbl_ps_nest
#' @importFrom tidyr unnest
#' @export
unnest.tbl_ps_nest <- function(data, cols, ..., keep_empty = FALSE, ptype = NULL,
    names_sep = NULL, names_repair = "check_unique"){
    x <- data
    class(x) <- "data.frame"
    if (missing(cols)||missing(...)){
        cols <- colnames(data)[unlist(lapply(data, function(i)rlang::is_list(i)))]
    }
    res <- unnest(x, {{cols}}, keep_empty=keep_empty, 
                  ptype=ptype, names_sep=names_sep, names_repair=names_repair)
    res <- add_attr.tbl_ps(x1=res, x2=data)
    return (res)
}

#' @method nest tbl_ps
#' @importFrom tidyr nest
#' @export
nest.tbl_ps <- function(.data, ..., .names_sep = NULL){
    x <- nest_internal(x=.data, ..., .names_sep=.names_sep)
    x <- add_attr.tbl_ps(x1=x, x2=.data)
    class(x) <- c("tbl_ps_nest", class(x))
    return(x)
}

#' @method nest grouped_df_ps
#' @export
nest.grouped_df_ps <- function(.data, ..., .names_sep = NULL){
    x <- nest_internal(x=.data, ..., .names_sep = .names_sep)
    x <- add_attr.tbl_ps(x1=x, x2=.data)
    class(x) <- c("tbl_ps_nest", class(x))
    return(x)
}

#' @importFrom tidyr nest
nest_internal <- function(x, ..., .names_sep = NULL){
    class(x) <- class(x)[!class(x) %in% c("grouped_df_ps", "tbl_ps_nest", "tbl_ps")]
    if (missing(...)){
        clnm <- colnames(x)
        taxavar <- attr(x, "taxavar")
        clnm <- clnm[!clnm %in% c("OTU", taxavar)]
        params <- c(list(x), lapply(clnm, function(x)x))
        names(params) <- c(".data", clnm)
    }else{
        res <- nest(.data=x, ..., .names_sep=.names_sep)
        return(res)
    }
    if (!is.null(.names_sep)){
        params <- c(params, .names_sep=.names_sep)
    }
    res <- do.call("nest", params)
    return(res)
}


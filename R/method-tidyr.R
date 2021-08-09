#' @method unnest tbl_mpse_nest
#' @importFrom tidyr unnest
#' @export
unnest.tbl_mpse_nest <- function(data, cols, ..., keep_empty = FALSE, ptype = NULL,
    names_sep = NULL, names_repair = "check_unique"){
    x <- data
    class(x) <- "data.frame"
    if (missing(cols)||missing(...)){
        cols <- colnames(data)[unlist(lapply(data, function(i)rlang::is_list(i)))]
    }
    res <- unnest(x, {{cols}}, keep_empty=keep_empty, 
                  ptype=ptype, names_sep=names_sep, names_repair=names_repair)
    res <- add_attr.tbl_mpse(x1=res, x2=data)
    return (res)
}

#' @method nest tbl_mpse
#' @importFrom tidyr nest
#' @export
nest.tbl_mpse <- function(.data, ..., .names_sep = NULL){
    x <- nest_internal(x=.data, ..., .names_sep=.names_sep)
    x <- add_attr.tbl_mpse(x1=x, x2=.data)
    class(x) <- c("tbl_mpse_nest", class(x))
    return(x)
}

#' @method nest grouped_df_mpse
#' @export
nest.grouped_df_mpse <- function(.data, ..., .names_sep = NULL){
    x <- nest_internal(x=.data, ..., .names_sep = .names_sep)
    x <- add_attr.tbl_mpse(x1=x, x2=.data)
    class(x) <- c("tbl_mpse_nest", class(x))
    return(x)
}

#' @importFrom tidyr nest
nest_internal <- function(x, ..., .names_sep = NULL){
    class(x) <- class(x)[!class(x) %in% c("grouped_df_mpse", "tbl_mpse_nest", "tbl_mpse")]
    if (missing(...)){
        idx <- x %>% vapply(is.list, logical(1))
        clnm <- colnames(x)
		clnm <- clnm[!idx]
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


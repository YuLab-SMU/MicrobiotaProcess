formatted_out <- function(x){
    #x1 <- "\033[38;5;246m"
    #x2 <- "\033[39m"
    #mg <- paste0(c(x1, x, x2), collapse="") %>% 
    mg <- pillar::style_subtle(x)
    return(mg)
}

tbl_mpse_return_message <- function(flag){
    x1 <- "# Note: MPSE object is converted to a tibble data "
    x2 <- "(tbl_mpse object) "
    x3 <- "for independent data analysis."
    if (flag){
        x <- paste0(c(x1, x2, x3), collapse="")
    }else{
        x <- paste0(c(x1, x3), collapse="")
    }
    formatted_out(x)
}

keep_mpse_message <- function(x = '\nA new MPSE object can be returned by setting keep.mpse = TRUE.'){
    formatted_out(x)
}


drop_class <- function(x, class){
    old <- class(x)
    class(x) <- old[!old %in% class]
    return (x)
}

add_internal_attr <- function(data, object, name){
    if(!"internal_attr" %in% (data %>% attributes() %>% names())){
        data %<>% add_attr(list(), name="internal_attr")
    }
    internals <- data %>% attr("internal_attr")
    internals[[name]] <- object
    data %<>% add_attr(internals, "internal_attr")
    return(data)
}


add_attr <- function(x, attribute, name){
    attr(x, name) <- attribute
    return(x)
}

left_join_msg <- paste0("The 'by' has some conditions, it should be specified 'OTU' or 'Sample' via by='OTU' or by=c('OTU'='id'),\n",
                       "by='Sample' or by=c('Sample'='sid').\n",
                       "Or the 'OTU' and 'Sample' can also be specified simultaneously via by=c('OTU', 'Sample') or\n",
                       "by = c('OTU'='id', 'Sample'='sid').")

.internal_nest <- function(x, keepnm, ..., .names_sep = NULL){
    nest <- utils::getFromNamespace("nest", "tidyr")
    if (missing(...)){
        clnm <- colnames(x)
        clnm <- clnm[!clnm %in% keepnm]
        params <- c(list(x), lapply(clnm, function(x)x))
        names(params) <- c(".data", clnm)
    }else{
        res <- nest(.data=x, ..., .names_sep=.names_sep)
        return(res)
    }
    if (!is.null(.names_sep)){
        params <- c(params, .names_sep=.names_sep)
    }
    res <- do.call(nest, params)
    return(res)
}

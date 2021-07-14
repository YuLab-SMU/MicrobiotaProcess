formatted_out <- function(x){
    x1 <- "\033[38;5;246m"
    x2 <- "\033[39m"
    paste0(c(x1, x, x2), collapse="")
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

keep_mpse_message <- formatted_out("#       A new MPSE object can be returned by setting .returnMPSE = TRUE.")


drop_class <- function(x, class){
    old <- class(x)
    class(x) <- old[!old %in% class]
    return (x)
}

add_internals_attr <- function(data, object, name){
    if(!"internals_attr" %in% (data %>% attributes() %>% names())){
        data %<>% add_attr(list(), name="internals_attr")
    }
    internals <- data %>% attr("internals_attr")
    internals[[name]] <- object
    data %<>% add_attr(internals, "internals_attr")
    return(data)
}


add_attr <- function(x, attribute, name){
    attr(x, name) <- attribute
    return(x)
}

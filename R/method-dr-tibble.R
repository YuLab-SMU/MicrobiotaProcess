#' @title convert the results of multiple dimensionality reduction to tbl_df
#' @rdname tidydr
#' @param x the R object, options 'prcomp', 'princomp', 'pcoa'
#' @param display character  
#' @param digits integer indicating the number of decimal places,
#' default is 2.
#' @param ... additional parameters
#' @return tbl_df object
#' @noRd
tidydr <- function(x, ...){
    UseMethod("tidydr")
}

#' @method tidydr prcomp
#' @importFrom stats setNames
#' @rdname tidydr
#' @noRd
tidydr.prcomp <- function(x, display="sites", digits=2, ...){
    vars <- x$sdev^2 
    vars <- (100 * vars/sum(vars)) %>%
             round(digits=digits) %>%
             paste0("(%)")
    
    if (inherits(x, "prcomp")){
        sites <- "x"
        features <- "rotation"
    }else if (inherits(x, "princomp")){
        sites <- "scores"
        features <- "loadings"
    }

    da <- x[[sites]] %>% 
          as.data.frame() %>%
          setNames(
                   x[[sites]] %>% 
                   colnames() %>% 
                   paste(vars)
          ) %>%
          tibble::as_tibble(rownames="sites")

    if ("features" %in% display){
        da %<>% 
            add_attr(
                     attribute=x[[features]] %>% 
                     tibble::as_tibble(rownames="features"),
                     name="features_tb"
                )
    }

    return(da)
}

#' @method tidydr princomp
#' @rdname tidydr
#' @noRd
tidydr.princomp <- tidydr.prcomp

#' @method tidydr pcoa
#' @rdname tidydr
#' @noRd
tidydr.pcoa <- function(x, digits=2, ...){
    vars <- (100 * x$values$Relative_eig) %>% 
             round(digits=digits) %>% 
             paste0("(%)")

    da <- x$vectors %>%
          as.data.frame() %>%
          setNames(
              paste0("PCo", 
                     x$vectors %>% 
                         ncol %>% 
                         seq_len
                     ) %>% 
              paste(vars)
          ) %>%
          tibble::as_tibble(rownames="sites")

    return(da)
}

#' @method tidydr decorana
#' @rdname tidydr
#' @noRd
tidydr.decorana <- function(x, display="sites", digits=2, ...){
    da <- vegan::scores(x, display="sites", ...) %>%
          as.data.frame() %>%
          tibble::as_tibble(rownames="sites")

    if ("features" %in% display){
        dat <- vegan::scores(x, display="species", ...) %>%
               tibble::as_tibble(rownames="features")
        da %<>% 
            add_attr(dat, name="features_tb")
    }
    return(da)
}

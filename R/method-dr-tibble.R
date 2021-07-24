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
    vars <- (100 * x$values$Relative_eig[seq_len(ncol(x$vectors))]) %>% 
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

#' @method tidydr metaMDS
#' @rdname tidydr
#' @noRd
tidydr.metaMDS <- function(x, display="sites", ...){
    da <- x %>% 
          vegan::scores(display="sites",...) %>%
          as.data.frame() %>%
          tibble::as_tibble(rownames="sites") 

    if ("features" %in% display){
        dat <- tryCatch(
                   vegan::scores(x, display="species", ...),
                   error=function(i){NULL}
               ) %>%
               tibble::as_tibble(rownames="features")
        da %<>% 
            add_attr(dat, name="features_tb")
    }
    return(da)
}

#' @method tidydr cca
#' @rdname tidydr
#' @noRd
tidydr.cca <- function(x, display=c("sp", "wa", "lc", "bp", "cn"), digits=2, ...){
    vars <- x %>% 
            vegan::eigenvals() %>% 
            summary %>% 
            data.frame %>% 
            dplyr::slice(2) %>% 
            as.numeric() * 100
    vars %<>% 
        round(digits=digits) %>% 
        paste0("(%)")

    if (length(display)<=1){
        display <- "sp"
    }
    dalist <- x %>% 
              vegan::scores(display, choices=seq_len(length(vars)), ...)

    da <- dalist$sites %>%
          as.data.frame() %>%
          setNames(
                   dalist$sites %>%
                   colnames() %>%
                   paste(vars)
                  ) %>%
          tibble::as_tibble(rownames="sites")

    dalist <- dalist[-match("sites", names(dalist))]
    dalist <- dalist[vapply(dalist, function(i)!all(is.na(i)), logical(1))]
    if (length(dalist)>0){
        leftnm <- names(dalist)
        for (i in leftnm){
            da %<>%
                add_attr(dalist[[i]] %>% 
                         tibble::as_tibble(rownames=i),
                         name=i)
        }
    }
    return(da)
}

#' @method tidydr rda
#' @rdname tidydr
#' @noRd
tidydr.rda <- tidydr.cca

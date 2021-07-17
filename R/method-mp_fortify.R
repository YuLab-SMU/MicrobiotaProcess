##' Fortify a model with data in MicrobiotaProcess
##' @title mp_fortify
##' @param model object
##' @param ... additional parameters
##' @return data frame or tbl_df object
##' @rdname mp_fortify
##' @export
mp_fortify <- function(model, ...){
    UseMethod("mp_fortify")
}

#' @method mp_fortify envfit
#' @export
mp_fortify.envfit <- function(model, ...){
    da <- lapply(c("vectors", "factors"), function(i){
               res <- vegan::scores(model, display=i)
               res <- do.call("cbind", c(list(res), model[[i]][c(2,4)])) %>%
                      as_tibble(rownames="label") %>%
                      mutate(type=i)
               return(res)
            }) %>% 
          dplyr::bind_rows()
    return(da)
}

#' @method mp_fortify adonis
#' @export
mp_fortify.adonis <- function(model, ...){
    da <- model$aov.tab %>% 
          as.matrix() %>% 
          tibble::as_tibble(rownames="factors")
    return(da)
}

#' @method mp_fortify anosim
#' @export
mp_fortify.anosim <- function(model, ...){
    cat("ANOSIM statistic R: ")
    cat(formatC(model$statistic, digits = max(3, getOption("digits") - 3)), "\n")
    nperm <- model$permutations
    if (nperm) {
        cat("      Significance:", format.pval(model$signif), "\n\n")
        cat(howHead(model$control))
    }
    cat("\n")
    da <- tibble::tibble(class=model$class.vec, rank=model$dis.rank) %>%
          arrange(!!as.symbol("class"))
    return(da)
}


howHead <- getFromNamespace("howHead", "vegan")

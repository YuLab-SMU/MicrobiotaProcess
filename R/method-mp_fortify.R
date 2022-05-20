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

#' @method mp_fortify anova.cca
#' @export
mp_fortify.anova.cca <- function(model, ...){
    da <- model %>% 
          base::as.data.frame() %>% 
          tibble::as_tibble(rownames = "factors")
    return(da)
}

#' @method mp_fortify anosim
#' @export
mp_fortify.anosim <- function(model, verbose=FALSE, ...){
    if (verbose){
        cat("ANOSIM statistic R: ")
        cat(formatC(model$statistic, digits = max(3, getOption("digits") - 3)), "\n")
        nperm <- model$permutations
        if (nperm) {
            cat("      Significance:", format.pval(model$signif), "\n\n")
            cat(howHead(model$control))
        }
        cat("\n")
    }
    da <- tibble::tibble(class=model$class.vec, rank=model$dis.rank) %>%
          arrange(!!as.symbol("class"))
    return(da)
}

#' @method mp_fortify mrpp
#' @export
mp_fortify.mrpp <- function(model, verbose=FALSE, ...){
    if (verbose){
        cat("Dissimilarity index:", model$distance, "\n")
        cat("Weights for groups: ", switch(model$weight.type, "n", "n-1",
            "n(n-1)", "n(n-1)/2"), "\n\n")
        cat("Class means and counts:\n\n")
        print(noquote(rbind(delta = formatC(model$classdelta, digits = 4),
            n = formatC(model$n, digits = 0))))
        cat("\n")
        if (!is.na(model$CS)) {
            cat("Classification strength: ")
            cat(formatC(model$CS, digits = 4), "\n")
        }
        cat("Chance corrected within-group agreement A: ")
        if (!is.na(model$A))
            cat(formatC(model$A, digits = 4), "\n")
        else cat("NA\n")
        cat("Based on observed delta", formatC(model$delta), "and expected delta",
            formatC(model$E.delta), "\n\n")
        nperm <- model$permutations
        if (nperm) {
            cat("Significance of delta:", format.pval(model$Pvalue),"\n")
        }
        cat(howHead(model$control))
        cat("\n")    
    }
    da <- model[names(model) %in% c("classdelta", "E.delta", "delta", "Pvalue", "A")]
    da <- do.call("cbind", da) %>%
          as_tibble(rownames="group")
    return(da)
}

#' @method mp_fortify mantel
#' @export
mp_fortify.mantel <- function(model, ...){
    da <- model[names(model) %in% c("statistic", "signif")]
    da[["QuantilesOfPerm"]] <- stats::quantile(model$perm, c(0.9, 0.95, 0.975, 0.99))
    da <- do.call("cbind", da) %>%
          as_tibble(rownames="quantile")
    return(da)
    
}

#' @method mp_fortify mantel.partial
#' @export
mp_fortify.mantel.partial <- function(model, ...){
    da <- NextMethod()
    return(da)
}

howHead <- getFromNamespace("howHead", "vegan")

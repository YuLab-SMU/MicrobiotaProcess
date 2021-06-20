#' @method as_tibble phyloseq
#' @importFrom tibble rownames_to_column 
#' @importFrom reshape melt
#' @importFrom dplyr full_join
#' @importFrom tidyr pivot_longer
#' @export
as_tibble.phyloseq <- function(x, ...){
    otuda <- checkotu(x) %>% 
             rownames_to_column(var="Sample") %>%
             pivot_longer(!.data$Sample, names_to = "OTU", values_to = "Abundance") 
    sampleda <- get_sample(x) %>% rownames_to_column(var="Sample.idy")
    if (!is.null(sampleda)){
        otuda <- otuda %>% full_join(sampleda, by=c("Sample"="Sample.idy"))
        samplevar <- colnames(sampleda)
        samplevar[samplevar == "Sample.idy"] <- "Sample"
    }else{
        samplevar <- NULL
    }
    taxada <- as.data.frame(x@tax_table)
    if (!all(dim(taxada)==0)){
        taxavar <- colnames(taxada)
        taxada %<>% fillNAtax() %>% rownames_to_column(var="OTU")
        otuda <- otuda %>% full_join(taxada, by=c("OTU"="OTU"))
        fillNAtax <- TRUE
    }else{
        taxavar <- NULL
        fillNAtax <- NULL
    }
    attr(otuda, "samplevar") <- samplevar
    attr(otuda, "taxavar") <- taxavar
    attr(otuda, "fillNAtax") <- fillNAtax
    attr(otuda, "tree") <- x@phy_tree
    attr(otuda, "refseq") <- x@refseq
    class(otuda) <- c("tbl_ps", class(otuda))
    return(otuda)
}

#' @method as_tibble grouped_df_ps
#' @export
as_tibble.grouped_df_ps <- function(x, ...){
    res <- NextMethod()
    res <- drop_class(x=res, class=c("grouped_df_ps", "grouped_df"))
    return(res)
}

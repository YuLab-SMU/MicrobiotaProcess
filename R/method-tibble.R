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
    taxada <- as.data.frame(x@tax_table) %>% rownames_to_column(var="OTU")
    if (!is.null(taxada)){
        taxavar <- colnames(taxada)
        otuda <- otuda %>% full_join(taxada, by=c("OTU"="OTU"))
    }else{
        taxavar <- NULL
    }
    attr(otuda, "samplevar") <- samplevar
    attr(otuda, "taxavar") <- taxavar
    attr(otuda, "tree") <- x@phy_tree
    attr(otuda, "refseq") <- x@refseq
    class(otuda) <- c("tbl_ps", class(otuda))
    return(otuda)
}

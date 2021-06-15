#' @method as_tibble phyloseq
#' @importFrom tibble rownames_to_column 
#' @importFrom reshape melt
#' @importFrom dplyr full_join
#' @export
as_tibble.phyloseq <- function(x, ...){
    otuda <- checkotu(x) %>% 
             rownames_to_column(var="Sample") %>%
             melt(id="Sample", variable_name="OTU") %>%
             rename(Abundance="value") %>%
             as_tibble()
    sampleda <- get_sample(x) %>% rownames_to_column(var="Sample.idy")
    if (!is.null(sampleda)){
        otuda <- otuda %>% full_join(sampleda, by=c("Sample"="Sample.idy"))
        samplevar <- colnames(sampleda)
        samplevar[samplevar == "Sample.idy"] <- "Sample"
    }else{
        samplevar <- NULL
    }
    taxada <- as.data.frame(tax_table(x)) %>% rownames_to_column(var="OTU")
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


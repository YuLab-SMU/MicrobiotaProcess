#' @title convert to phyloseq object. 
#' @param x object, tbl_ps object, which the result of 
#' as_tibble for phyloseq objcet.
#' @param ... additional params
#' @return phyloseq object.
#' @export
as.phyloseq <- function(x, ...){
    UseMethod("as.phyloseq")
}

#' @rdname as.phyloseq
#' @export
as_phyloseq <- as.phyloseq

#' @method as.phyloseq tbl_ps
#' @rdname as.phyloseq
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr distinct
#' @export
as.phyloseq.tbl_ps <- function(x, ...){
    samplevar <- attr(x, "samplevar")
    taxavar <- attr(x, "taxavar")
    tree <- attr(x, "tree")
    refseq <- attr(x, "refseq")
    otuda <- x[,colnames(x) %in% c("Sample", "OTU", "Abundance")]
    otuda <- otuda %>% 
             distinct() %>%
             pivot_wider(names_from=.data$OTU, values_from=.data$Abundance) %>% 
             column_to_rownames(var="Sample")
    if (!is.null(samplevar)){
        sampleda <- x[, colnames(x) %in% samplevar] %>%
                    distinct() %>% 
                    column_to_rownames(var="Sample")
        sampleda <- sample_data(sampleda)
    }else{
        sampleda <- NULL
    }

    if (!is.null(taxavar)){
        taxada <- x[,colnames(x) %in% taxavar] %>% 
                  distinct() %>%
                  column_to_rownames(var="OTU") %>% 
                  as.matrix()
        taxada <- tax_table(taxada)
    }else{
        taxada <- NULL
    }
    
    if (!is.null(tree)){
        tree <- keep.tip(tree, .data$OTU)
    }
    if (!is.null(refseq)){
        refseq <- refseq[.data$OTU]
    }
    res <- phyloseq(otu_table(otuda, taxa_are_rows=FALSE), sampleda, taxada, tree, refseq)
    return (res)
}

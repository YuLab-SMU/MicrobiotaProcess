#' @title convert to phyloseq object. 
#' @param x object, tbl_mpse object, which the result of 
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

#' @method as.phyloseq tbl_mpse
#' @rdname as.phyloseq
#' @importFrom tidyr pivot_wider
#' @importFrom dplyr distinct
#' @export
as.phyloseq.tbl_mpse <- function(x, ...){
    samplevar <- attr(x, "samplevar")
    taxatree <- attr(x, "taxatree")
    tree <- attr(x, "otutree")
    refseq <- attr(x, "refseq")
    otuda <- x %>% 
             select(c("Sample", "OTU", "Abundance")) %>%
             distinct() %>%
             pivot_wider(names_from="OTU", values_from="Abundance") %>% 
             column_to_rownames(var="Sample")
    if (!is.null(samplevar)){
        sampleda <- x[, colnames(x) %in% samplevar] %>%
                    distinct() %>% 
                    column_to_rownames(var="Sample")
        sampleda <- phyloseq::sample_data(sampleda)
    }else{
        sampleda <- NULL
    }

    if (!is.null(taxatree)){
        #taxada <- x[,colnames(x) %in% c("OTU",taxavar)] %>% 
        #          distinct() %>%
        #          column_to_rownames(var="OTU") %>% 
        #          as.matrix()
        taxada <- taxatree_to_tb(taxatree)
        taxada <- phyloseq::tax_table(taxada)
    }else{
        taxada <- NULL
    }
    
    if (!is.null(tree)){
        tree <- keep.tip(tree, colnames(otuda))
    }
    if (!is.null(refseq)){
        refseq <- refseq[colnames(otuda)]
    }
    res <- phyloseq::phyloseq(phyloseq::otu_table(otuda, taxa_are_rows=FALSE), sampleda, taxada, tree, refseq)
    return (res)
}

#' @method as.phyloseq grouped_df_mpse
#' @export
as.phyloseq.grouped_df_mpse <- function(x, ...){
    x %<>% ungroup()
    res <- as.phyloseq(x=x,...)
    return(res)
}

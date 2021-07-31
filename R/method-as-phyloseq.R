#' @title convert to phyloseq object. 
#' @param x object, tbl_mpse object, which the result of 
#' as_tibble for phyloseq objcet.
#' @param .abundance the column name to be as the abundance 
#' of otu table, default is Abundance.
#' @param ... additional params
#' @return phyloseq object.
#' @export
as.phyloseq <- function(x, .abundance, ...){
    UseMethod("as.phyloseq")
}

#' @rdname as.phyloseq
#' @export
as_phyloseq <- as.phyloseq

#' @method as.phyloseq MPSE
#' @rdname as.phyloseq
#' @export

as.phyloseq.MPSE <- function(x, .abundance, ...){
    .abundance <- rlang::enquo(.abundance)
    
    if (rlang::quo_is_missing(.abundance)){
        .abundance <- as.symbol("Abundance")
    }
    otuda <- x %>% 
              mp_extract_assays(.abundance=!!.abundance) %>%
              phyloseq::otu_table(taxa_are_rows=TRUE)     
    
    otutree <- x %>% mp_extract_tree(type="otutree")
    sampleda <- x %>% mp_extract_sample()
    taxada <- x %>% mp_extract_taxonomy()

    if (inherits(x, "MPSE")){
        refseq <- x@refseq
    }else{
        refseq <- x %>% attr("refseq")
    }

    if (ncol(sampleda) > 1){
        sampleda <- sampleda %>% 
                    tibble::column_to_rownames(var="Sample") %>%
                    phyloseq::sample_data()
    }else{
        sampleda <- NULL
    }

    if (!is.null(taxada)){
        taxada <- taxada %>%
                as.matrix() %>% 
                phyloseq::tax_table()
    }else{
        taxada <- NULL
    }

    if (!is.null(otutree)){
        otutree <- otutree@phylo
    }

    res <- phyloseq::phyloseq(otuda, sampleda, taxada, otutree, refseq)
    return (res)    
}

#' @method as.phyloseq tbl_mpse
#' @rdname as.phyloseq
#' @export
as.phyloseq.tbl_mpse <- as.phyloseq.MPSE

#' @method as.phyloseq grouped_df_mpse
#' @export
as.phyloseq.grouped_df_mpse <- as.phyloseq.MPSE

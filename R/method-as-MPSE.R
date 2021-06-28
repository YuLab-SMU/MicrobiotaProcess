##' @docType methods
##' @name as.MPSE
##' @rdname as.MPSE-methods
##' @title as.MPSE method
##' @param .data tbl_mpse object
##' @param ... additional parameters.
##' @return MPSE object
##' @export
setGeneric("as.MPSE", function(.data, ...){standardGeneric("as.MPSE")})



##' @rdname as.MPSE-methods
##' @aliases as.MPSE,tbl_mpse
##' @exportMethod as.MPSE
setMethod("as.MPSE", signature(.data="tbl_mpse"),
    function(.data, ...){

        otutree <- attr(.data, "otutree")
        refseq <- attr(.data, "refseq")
        samplevar <- attr(.data, "samplevar")
        taxavar <- attr(.data, "taxavar")
        assaysvar <- attr(.data, "assaysvar")

        .data %<>% as_tibble()

        assaysda <- lapply(assaysvar, function(x) 
            .data %>%
            select(c("Sample", "OTU", x)) %>%
            distinct() %>%
            pivot_wider(names_from="Sample", values_from=x) %>%
            column_to_rownames(var="OTU")
        ) %>% stats::setNames(assaysvar)

        if (!is.null(samplevar)){
            sampleda <- .data %>% 
                        select(samplevar) %>%
                        distinct() %>%
                        column_to_rownames(var="Sample") 
                        
        }else{
            samplelist <- .data %>% select("Sample") %>% 
                           unlist(use.names=FALSE) %>% as.vector()
            sampleda <- S4Vectors::DataFrame(NULL, row.names=samplelist) 
        }

        if (!is.null(otutree)){
            rmtip <- setdiff(otutree@phylo$tip.label, rownames(assaysda[[1]]))
            otutree <- treeio::drop.tip(otutree, tip=rmtip)
        }

        if (!is.null(refseq)){
            refseq <- refseq[rownames(assaysda[[1]])]
        }

        mpse <- MPSE(
                    assays  = assaysda,
                    colData = sampleda,
                    otutree = otutree,
                    refseq  = refseq,
                )

        if (!is.null(taxavar)){
            taxada <- .data %>%
                      select(c("OTU",taxavar)) %>%
                      distinct() %>%
                      column_to_rownames(var="OTU")
            SummarizedExperiment::rowData(mpse) <- taxada
        }

        return (mpse)
})


##' @rdname as.MPSE-methods
##' @aliases as.MPSE,phyloseq
##' @exportMethod as.MPSE
setMethod("as.MPSE", signature(.data="phyloseq"), 
    function(.data, ...){

        otuda <- checkotu(.data) %>% as.matrix() %>% t()
        sampleda <- get_sample(.data) 
        taxada <- as.data.frame(.data@tax_table)
        otutree <- NULL
        #taxatree <- NULL

        if (ncol(taxada)!=0){
            taxada %<>% fillNAtax()
            taxada$OTU <- rownames(taxada)
            taxatree <- convert_to_treedata(data=taxada)
            taxada$OTU <- NULL
        }

        if (!is.null(.data@phy_tree)){
            otutree <- .data@phy_tree %>% as.treedata()
        }

        if (is.null(sampleda)){
            sampleda <- S4Vectors::DataFrame(row.names=colnames(otuda))
        }

        mpse <- MPSE(
                   colData  = sampleda,
                   assays   = list(Abundance=otuda),
                   otutree  = otutree,
                   #taxatree = taxatree,
                   refseq   = .data@refseq
                )

        if (ncol(taxada)!=0){
            SummarizedExperiment::rowData(mpse) <- taxada
        }

        return(mpse)
})



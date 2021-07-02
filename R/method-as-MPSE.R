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
        #taxavar <- attr(.data, "taxavar")
        otumetavar <- attr(.data, "otumetavar")
        taxatree <- attr(.data, "taxatree")
        assaysvar <- attr(.data, "assaysvar")
        mutatevar <- attr(.data, "mutatevar")

        .data %<>% as_tibble()
        assaysda <- lapply(assaysvar, function(x) 
            .data %>%
            select(c("Sample", "OTU", x)) %>%
            distinct() %>%
            pivot_wider(names_from="Sample", values_from=x) %>%
            column_to_rownames(var="OTU")
        ) %>% stats::setNames(assaysvar)

        nsample <- ncol(assaysda[[1]])
        if (!is.null(mutatevar)){
            indx <- lapply(mutatevar, function(i)
                        .data %>% 
                        select(i) %>% 
                        distinct() %>% 
                        nrow() == nsample
                        ) %>% unlist()
            mutatevar <- mutatevar[indx]
            samplevar <- c(samplevar, mutatevar)
        }
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
            if (length(rmtip)>0){
                otutree <- treeio::drop.tip(otutree, tip=rmtip)
            }
        }
        if (!is.null(taxatree)){
            rmtip <- setdiff(taxatree@phylo$tip.label, rownames(assaysda[[1]]))
            if (length(rmtip)>0){
                taxatree <- treeio::drop.tip(taxatree, 
                                             tip=rmtip, 
                                             collapse.singles=FALSE) 
            }
        }

        if (!is.null(refseq)){
            refseq <- refseq[rownames(assaysda[[1]])]
        }

        mpse <- MPSE(
                    assays  = assaysda,
                    colData = sampleda,
                    otutree = otutree,
                    refseq  = refseq,
                    taxatree = taxatree
                )

        if (!is.null(otumetavar)){
            otumeta <- .data %>%
                      select(c("OTU", otumetavar)) %>%
                      distinct() %>%
                      column_to_rownames(var="OTU")
            SummarizedExperiment::rowData(mpse) <- otumeta
        }
        methods::validObject(mpse)
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
        taxatree <- NULL

        if (ncol(taxada)!=0){
            taxada %<>% fillNAtax()
            taxatree <- convert_to_treedata2(x=taxada)
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
                   taxatree = taxatree,
                   refseq   = .data@refseq
                )

        #if (ncol(taxada)!=0){
        #    SummarizedExperiment::rowData(mpse) <- taxada
        #}
        methods::validObject(mpse)
        return(mpse)
})

##' @rdname as.MPSE-methods
##' @aliases as.MPSE,grouped_df_mpse
##' @exportMethod as.MPSE
setMethod("as.MPSE", signature(.data="grouped_df_mpse"), 
          function(.data, ...){
    mpse <- .data %>% ungroup() %>% as.MPSE()
    methods::validObject(mpse)
    return(mpse)          
})

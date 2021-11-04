##' @title as.MPSE method
##' @description 
##' convert the .data object to MPSE object
##' @param .data one type of tbl_mpse, phyloseq, biom, SummarizedExperiment or TreeSummarizedExperiment class
##' @param ... additional parameters, meaningless now.
##' @return MPSE object
##' @export
##' @author Shuangbin Xu
##' @examples
##' data(test_otu_data)
##' test_otu_data %>% as.MPSE -> mpse
##' mpse
as.MPSE <- function(.data, ...){
    if (inherits(.data, "MPSE")){
        return (.data)
    }else if (inherits(.data, "grouped_df_mpse")){
        res <- .as.MPSE.grouped_df_mpse(.data, ...)
    }else if (inherits(.data, "tbl_mpse")){
        res <- .as.MPSE.tbl_mpse(.data, ...)
    }else if (inherits(.data, "phyloseq")){
        res <- .as.MPSE.phyloseq(.data, ...)
    }else if (inherits(.data, "SummarizedExperiment")){
        res <- .as.MPSE.TSE(.data, ...)
    }else if (inherits(.data, "biom")){
        res <- .as.MPSE.biom(.data, ...)
    }else {
        message("The as.MPSE now only works for tbl_mpse, grouped_df_mpse, phyloseq, biom, SummarizedExperiment, and TreeSummarizedExperiment class.")
        return()
    }
    return (res)
}

.as.MPSE.tbl_mpse <- function(.data, ...){

        otutree <- attr(.data, "otutree")
        refseq <- attr(.data, "refseq")
        samplevar <- attr(.data, "samplevar")
        #taxavar <- attr(.data, "taxavar")
        otumetavar <- attr(.data, "otumetavar")
        taxatree <- attr(.data, "taxatree")
        assaysvar <- attr(.data, "assaysvar")
        mutatevar <- attr(.data, "mutatevar")

        internal_attr <- attr(.data, "internal_attr")

        .data %<>% as_tibble()
        assaysda <- .internal_build_assay(da=.data, x=assaysvar)

        if (length(mutatevar)>0){
            res.var <- .internal_check_mutate(da=.data, 
                                              x=mutatevar, 
                                              assay=assaysda[[1]])
            samplevar <- union(samplevar, res.var$samplevar)
            otumetavar <- union(otumetavar, res.var$otumetavar)
            newassay <- setdiff(res.var$assaysvar, assaysvar)
            if (length(newassay) > 0) {
                assaysda <- c(assaysda, .internal_build_assay(da=.data, x=newassay))
            }
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

        otutree <- .internal_drop.tip(tree=otutree, newnm=rownames(assaysda[[1]]))

        taxatree <- .internal_drop.tip(
                                       tree=taxatree, 
                                       newnm=rownames(assaysda[[1]]),
                                       collapse.singles = FALSE
                                   )

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
        if (!is.null(internal_attr)){
            mpse %<>% add_attr(attribute=internal_attr, name="internal_attr")
        }
        methods::validObject(mpse)
        return (mpse)
}

.internal_build_assay <- function(da, x){
    assayda <- lapply(x, function(i)
        da %>%
            dplyr::select(c("Sample", "OTU", i)) %>%
            dplyr::distinct() %>%
            tidyr::pivot_wider(names_from="Sample", values_from=i) %>%
            tibble::column_to_rownames(var="OTU") %>% as.matrix()
        ) %>% stats::setNames(x)
    return(assayda)
}

.internal_check_mutate <- function(da, x, assay){
    indx1 <- lapply(x, function(i) 
               da %>%
               select(c("OTU", i)) %>%
               distinct() %>% 
               nrow() %>%
               magrittr::equals(nrow(assay))
             ) %>% unlist()
    otumetavar <- x[indx1]
    x <- x[!indx1]
    if (length(x)==0){
        return(list(otumetavar=otumetavar))
    }

    indx2 <- lapply(x, function(i)
               da %>%
               select(c("Sample", i)) %>%
               distinct() %>%
               nrow() %>%
               magrittr::equals(ncol(assay))
             ) %>% unlist()

    samplevar <- x[indx2]
    x <- x[!indx2]
    if (length(x)==0){
        return(list(otumetavar=otumetavar, samplevar=samplevar))
    }

    indx3 <- lapply(x, function(i)
               da %>% 
               select(c("OTU", "Sample", i)) %>%
               distinct() %>%
               nrow() %>%
               magrittr::equals(nrow(assay)*ncol(assay)) &&
               da %>% pull(i) %>% is.numeric()
             ) %>% unlist()

    return(list(otumetavar=otumetavar, samplevar=samplevar, assaysvar=x[indx3]))
}

.as.MPSE.phyloseq <- function(.data, ...){

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
}

.as.MPSE.grouped_df_mpse <- function(.data, ...){
    mpse <- .data %>% ungroup() %>% as.MPSE()
    methods::validObject(mpse)
    return(mpse)          
}


.as.MPSE.TSE <- function(.data, ...){

    assaysvar <- SummarizedExperiment::assayNames(.data)
    SummarizedExperiment::assayNames(.data) <- c("Abundance", assaysvar[-1])
    
    if (is.null(rownames(.data))){
        flag <- NULL
    }else{
        flag <- rownames(.data) %>% base::strsplit("\\|") %>% lapply(length) %>% unlist
    }
    rowda <- SummarizedExperiment::rowData(.data)
    res <- .internal_check_taxonomy(rowda, flag)
    taxatab <- res$taxatab
    newrowda <- res$newrowda
    if (all(flag>5)){
        #`rownames<-` <- methods::getMethod("rownames<-", "TreeSummarizedExperiment")
        if (inherits(.data, "TreeSummarizedExperiment")){
            rownames(.data) <- rownames(taxatab)
        }
    }

    if (inherits(.data, "TreeSummarizedExperiment")){
        rowTree <- getFromNamespace("rowTree", "TreeSummarizedExperiment")
        otu.tree <- rowTree(.data)
    }else if (inherits(.data, "SummarizedExperiment")){
        meta.da <- S4Vectors::metadata(.data)
        indx <- which(meta.da %>% lapply(class)=="phylo")
        if (length(indx) == 0){
            otu.tree <- NULL
        }else{
            otu.tree <- meta.da[[indx[1]]]
        }
    }
    
    if (!is.null(otu.tree)){
        if (all(flag>5)){
            otu.tree$tip.label %<>%
              base::strsplit("\\|") %>%
              lapply(., function(x)x[length(x)]) %>%
              unlist()
        }
        keepnms <- intersect(otu.tree$tip.label, rownames(.data))
        otu.tree <- ape::keep.tip(otu.tree, tip=keepnms)
        .data <- .data[rownames(.data) %in% keepnms, , drop=FALSE]
        taxatab <- taxatab[rownames(taxatab) %in% keepnms, , drop=FALSE]
        newrowda <- newrowda[rownames(newrowda) %in% keepnms, , drop=FALSE]
        otu.tree %<>% treeio::as.treedata()
    }

    if (!is.null(taxatab) && ncol(taxatab)>0){
        taxa.tree <- convert_to_treedata2(x=data.frame(taxatab))
    }else{
        taxa.tree <- NULL
    }

    if (inherits(.data, "TreeSummarizedExperiment")){
        referenceSeq <- getFromNamespace("referenceSeq", "TreeSummarizedExperiment")
        refseq <- referenceSeq(.data)[[1]]
    }else{
        refseq <- NULL
    }

    mpse <- MPSE(assays=SummarizedExperiment::assays(.data),
                 colData=SummarizedExperiment::colData(.data),
                 refseq=refseq,
                 otutree=otu.tree,
                 taxatree=taxa.tree
            ) 
    
    if (!is.null(newrowda) && ncol(newrowda)>0){
        SummarizedExperiment::rowData(mpse) <- newrowda
    }

    return(mpse)
}

.as.MPSE.biom <- function(.data, ...){
    x <- .internal_parse_biom(.data)
    if ( !is.null(x$taxatab)){
        taxa.tree <- convert_to_treedata2(x$taxatab)
    }else{
        taxa.tree <- NULL
    }
    
    mpse <- MPSE(assays = list(Abundance=x$otutab),
                 taxatree = taxa.tree
            )

    return(mpse)
}

.internal_check_taxonomy <- function(x, flag){
    #flag <- rownames(x) %>% strsplit("\\|") %>% lapply(length) %>% unlist 
    #col.ind <- colnames(x) %in% c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    if (is.null(flag)){
        return(NULL)
    }
    df <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    col.ind <- lapply(df, function(i)which(grepl(i, colnames(x), ignore.case=TRUE))) %>% unlist()
    if (all(flag>5)){
        max.sep <- max(flag)
        taxatab <- x %>% 
          rownames() %>% 
          data.frame() %>% 
          split_str_to_list(sep="\\|") %>% 
          magrittr::set_colnames(c(taxaclass[seq_len(max.sep-1)], "OTU")) %>%
          magrittr::set_rownames(paste0("row", seq_len(nrow(x)))) %>%
          fillNAtax() %>%
          magrittr::extract(paste0("row", seq_len(nrow(x))), ) %>%
          magrittr::set_rownames(NULL) %>%
          tibble::column_to_rownames(var="OTU")
    }else{
        taxatab <- x[, col.ind] %>% data.frame(check.names=FALSE)
        if (ncol(taxatab)>2){
          taxatab <- fillNAtax(taxatab)
        }
    }
    newrowda <- x[, !seq_len(ncol(x)) %in% col.ind, drop=FALSE]
    return(list(taxatab=taxatab, newrowda=newrowda))
}

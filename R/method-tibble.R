#' @method as_tibble phyloseq
#' @importFrom tibble rownames_to_column 
#' @importFrom dplyr left_join
#' @importFrom tidyr pivot_longer
#' @export
as_tibble.phyloseq <- function(x, ...){
    otuda <- checkotu(x) %>% 
             tibble::as_tibble(rownames="Sample") %>%
             pivot_longer(!.data$Sample, names_to = "OTU", values_to = "Abundance") 
    sampleda <- get_sample(x) %>%
                avoid_conflict_names() %>%
                rownames_to_column(var="Sample")
    if (nrow(sampleda)!=0){
        otuda <- otuda %>% left_join(sampleda, by="Sample")
        samplevar <- colnames(sampleda)
    }else{
        samplevar <- NULL
    }
    taxada <- as.data.frame(x@tax_table)
    if (!all(dim(taxada)==0)){
        taxavar <- colnames(taxada)
        taxada %<>% fillNAtax() %>%
                    avoid_conflict_names() %>% 
                    tibble::as_tibble(rownames="OTU")
        otuda <- otuda %>% left_join(taxada, by="OTU")
        fillNAtax <- TRUE
    }else{
        taxavar <- NULL
        fillNAtax <- NULL
    }
    if (!is.null(x@phy_tree)){
        otutree <- x@phy_tree %>% as.treedata() 
    }else{
        otutree <- NULL
    }
    attr(otuda, "samplevar") <- samplevar
    attr(otuda, "taxavar") <- taxavar
    attr(otuda, "assaysvar") <- "Abundance"
    attr(otuda, "fillNAtax") <- fillNAtax
    attr(otuda, "otutree") <- otutree
    attr(otuda, "refseq") <- x@refseq
    class(otuda) <- c("tbl_mpse", class(otuda))
    return(otuda)
}

#' @method as_tibble grouped_df_mpse
#' @export
as_tibble.grouped_df_mpse <- function(x, ...){
    res <- NextMethod()
    res <- add_attr.tbl_mpse(x1=res, x2=x)
    res <- drop_class(x=res, class=c("grouped_df_mpse", "grouped_df"))
    return(res)
}

#' @method as_tibble MPSE
#' @export
as_tibble.MPSE <- function(x, ...){
    res <- .as_tibble_MPSE(x=x, ...)
    return (res)
}

.as_tibble_MPSE <- function(x, .subset = NULL){
    .subset = rlang::enquo(.subset)
    otuda <- extract_count_data(x)
    sampleda <- 
        SummarizedExperiment::colData(x) %>%
        avoid_conflict_names() %>%
        tibble::as_tibble(rownames="Sample")
    if (nrow(sampleda)!=0){
        otuda <- otuda %>% left_join(sampleda, by="Sample")
        samplevar <- colnames(sampleda)
    }else{
        samplevar <- NULL
    }
    taxada <-
        SummarizedExperiment::rowData(x) %>%
        avoid_conflict_names() 
    if (ncol(taxada)!=0){
        taxavar <- colnames(taxada)
        taxada %<>% tibble::as_tibble(rownames="OTU")
        otuda <- otuda %>% left_join(taxada, by="OTU")
        fillNAtax <- TRUE
    }else{
        taxavar <- NULL
        fillNAtax <- TRUE
    }
    
    if (!rlang::quo_is_null(.subset)){
        otuda <- otuda %>% dplyr::select(!!.subset)
        samplevar <- samplevar[samplevar %in% colnames(otuda)]
        taxavar <- taxavar[taxavar %in% colnames(otuda)]
    }
    attr(otuda, "samplevar") <- samplevar
    attr(otuda, "taxavar") <- taxavar
    attr(otuda, "assaysvar") <- names(x@assays)
    attr(otuda, "fillNAtax") <- fillNAtax
    attr(otuda, "otutree") <- x@otutree
    #attr(otuda, "taxatree") <- x@taxatree
    attr(otuda, "refseq") <- x@refseq
    class(otuda) <- c("tbl_mpse", class(otuda))
    return(otuda)
}

avoid_conflict_names <- function(.data){
    cls <- colnames(.data)
    cls <- gsub("^OTU$", "OTU.x", cls)
    cls <- gsub("^Sample$", "Sample.x", cls)
    stats::setNames(.data, cls)
    return(.data)
}

extract_count_data <- function(SE_object){
    allda <- SummarizedExperiment::assays(SE_object) %>% as.list()
    clnm <- names(allda)
    if (any(is.na(clnm)) || length(clnm) != length(allda)){
        stop ("The names of assays should be provided !", call. = FALSE)
    }
    da <- mapply(trans_to_longer, 
              allda, 
              clnm, 
              SIMPLIFY=FALSE) %>%
          Reduce(left_join, .) %>%
          suppressMessages()
    return (da)
}

trans_to_longer <- function(.data, name){
    .data %<>% 
        tibble::as_tibble(rownames = "OTU", .name_repair = "minimal") %>%
        tidyr::pivot_longer(!.data$OTU, names_to = "Sample", values_to=name)
    return(.data)
}


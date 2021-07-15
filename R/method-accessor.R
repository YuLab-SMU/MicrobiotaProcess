##' @method [ diffAnalysisClass
##' @export
`[.diffAnalysisClass` <- function(x, i, j, asis = FALSE, ...) {
    result <- x@result
    y <- result[i, j, ...]
    if (!asis)
        return(y)
    x@result <- y
    return(x)
}

##' @method [[ diffAnalysisClass
##' @export
`[[.diffAnalysisClass` <- function(x, i) {
    result <- x@result
    if (!i %in% names(result))
        stop("input term not found...")
    result[[i]]
}

##' @method $ diffAnalysisClass
##' @export
`$.diffAnalysisClass` <-  function(x, name) {
    x <- x@result
    x[, name]
}

##' @method dim diffAnalysisClass
##' @export
dim.diffAnalysisClass <- function(x) {
    dim(x@result)
}

##' @importFrom utils tail
##' @method tail diffAnalysisClass
##' @export
tail.diffAnalysisClass <- function(x, n=6L, ...) {
    tail(x@result, n, ...)
}

#' @importFrom utils head
#' @method head diffAnalysisClass
#' @export
head.diffAnalysisClass <- function(x, n=6L, ...){
    head(as.data.frame(x), n=n, ...)
}

#' @method head alphasample
#' @export
head.alphasample <- function(x, n=6L, ...){
    head(as.data.frame(x), n=n, ...)
}

#' MPSE accessors
#' @param x MPSE object
#' @param i,j,... Indices specifying elements to extract or replace.
#' Indices are 'numeric' or 'character' vectors or empty (missing) or
#' NULL.  Numeric values are coerced to integer as by 'as.integer' 
#' (and hence truncated towards zero).  Character vectors will be matched 
#' to the 'names' of the object (or for matrices/arrays, the 'dimnames')
#' @param drop logical If 'TRUE' the result is coerced to the lowest 
#' possible dimension (see the examples).  This only works for extracting 
#' elements, not for the replacement.
#' @param value character vector, Either ‘NULL’ or a character vector equal 
#' of length equal to the appropriate dimension.
#' @name MPSE-accessors
NULL

#' @rdname MPSE-accessors
#' @export
setMethod("[", signature(x="MPSE"),
          function(x, i, j, ..., drop=TRUE){

    otutree <- x@otutree
    refseq <- x@refseq
    taxatree <- x@taxatree
    nx <- methods::callNextMethod()
    newotus <- rownames(nx)
   
    otutree <- .internal_drop.tip(tree=otutree, newnm=newotus)
 
    taxatree <- .internal_drop.tip(
                                       tree=taxatree, 
                                       newnm=newotus, 
                                       collapse.singles=FALSE
                                )

    if (!is.null(refseq)){
        refseq <- refseq[newotus]
    }
 
    nx@taxatree <- taxatree
    nx@otutree <- otutree
    nx@refseq <- refseq
    methods::validObject(nx)
    return(nx)
})

#' @method [ tbl_mpse
#' @export
`[.tbl_mpse` <- function(x, i, j, ..., drop=TRUE){
    res <- NextMethod()
    if (valid_names(res, type="tbl_mpse")){
        res <- add_attr.tbl_mpse(x1=res, x2=x)
    }else{
        res <- drop_class(res, class="tbl_mpse")
    }
    return (res)
}

#' extract the abundance matrix from MPSE object or tbl_mpse object
#' @rdname mp_extract_abundance-methods
#' @param x MPSE or tbl_mpse object
#' @param .abundance the name of abundance to be extracted.
#' @param byRow logical if it is set TRUE, 'otu X sample' shape will return,
#' else 'sample X otu' will return.
#' @return otu abundance a data.frame object
#' @param ... additional parameters.
#' @export
setGeneric("mp_extract_abundance", function(x, .abundance, byRow=TRUE, ...)standardGeneric("mp_extract_abundance"))

#' @rdname mp_extract_abundance-methods
#' @aliases mp_extract_abundance,MPSE
#' @exportMethod mp_extract_abundance
setMethod("mp_extract_abundance", signature(x="MPSE"), function(x, .abundance, byRow=TRUE, ...){

    .abundance <- rlang::enquo(.abundance)
    if (rlang::quo_is_missing(.abundance)){
        rlang::abort("The abundance name is required !")
    }

    assaysvar <- SummarizedExperiment::assayNames(x)
    if (!rlang::as_name(.abundance) %in% assaysvar){
        rlang::abort("The abundance is not in the object!")
    }

    xx <- SummarizedExperiment::assays(x)@listData
    dat <- xx[[rlang::as_name(.abundance)]]

    if (byRow){
        dat %<>% as.data.frame(check.names=FALSE)
    }else{
        dat %<>% t() %>% as.data.frame(check.names=FALSE)
    }
    return(dat)
})

.internal_extract_abundance <- function(x, .abundance, byRow=TRUE, ...){
    
    .abundance <- rlang::enquo(.abundance)

    if (rlang::quo_is_missing(.abundance)){
        rlang::abort("The abundance name is required !")
    }

    assaysvar <- x %>% attr("assaysvar")

    if (!rlang::as_name(.abundance) %in% assaysvar){
        rlang::abort("The abundance is not in the object!")
    }

    dat <- x %>% 
        ungroup() %>%
        select(c("OTU", "Sample", rlang::as_name(.abundance))) 

    if (byRow){
        dat %<>%
            tidyr::pivot_wider(id_cols="OTU", 
                               names_from="Sample", 
                               values_from=rlang::as_name(.abundance)) %>%
            tibble::column_to_rownames(var="OTU")
    }else{
        dat %<>%
            tidyr::pivot_wider(id_cols="Sample",
                               names_from="OTU",
                               values_from=rlang::as_name(.abundance)) %>%
            tibble::column_to_rownames(var="Sample")   
    }
    
    return(dat)
}

#' @rdname mp_extract_abundance-methods
#' @aliases mp_extract_abundance,tbl_mpse
#' @exportMethod mp_extract_abundance
setMethod("mp_extract_abundance", signature(x="tbl_mpse"), .internal_extract_abundance)

#' @rdname mp_extract_abundance-methods
#' @aliases mp_extract_abundance,grouped_df_mpse
#' @exportMethod mp_extract_abundance
setMethod("mp_extract_abundance", signature(x="grouped_df_mpse"), .internal_extract_abundance)


#' @title extract the taxonomy tree in MPSE object
#' @docType methods
#' @rdname mp_extract_tree-methods
#' @param x MPSE object
#' @param type character taxatree or otutree
#' @param ... additional arguments
#' @return taxatree treedata object
#' @export
setGeneric("mp_extract_tree", function(x, type="taxatree", ...){standardGeneric("mp_extract_tree")})


#' @rdname mp_extract_tree-methods
#' @aliases mp_extract_tree,MPSE
#' @exportMethod mp_extract_tree
setMethod("mp_extract_tree", signature(x="MPSE"), function(x, type="taxatree", ...){
    type %<>% match.arg(c("taxatree", "otutree"))
    tree <- slot(x, type)
    if (!is.null(tree)){
        return(tree)
    }else{
        message(tree_empty(type=type))
    }
})

#' @rdname mp_extract_tree-methods
#' @aliases mp_extract_tree,tbl_mpse
#' @exportMethod mp_extract_tree
setMethod("mp_extract_tree", signature(x="tbl_mpse"),function(x, type="taxatree", ...){
    .internal_tree(x=x, type=type)
})

#' @rdname mp_extract_tree-methods
#' @aliases mp_extract_tree,grouped_df_mpse
#' @exportMethod mp_extract_tree
setMethod("mp_extract_tree", signature(x="grouped_df_mpse"),function(x, type="taxatree", ...){
    .internal_tree(x=x, type=type)
})

#' @title extract the taxonomy annotation in MPSE object
#' @docType methods
#' @rdname mp_extract_taxonomy-methods
#' @param x MPSE object
#' @param ... additional arguments
#' @return data.frame contained taxonomy annotation.
#' @export
setGeneric("mp_extract_taxonomy", function(x, ...)standardGeneric("mp_extract_taxonomy"))

#' @rdname mp_extract_taxonomy-methods
#' @aliases mp_extract_taxonomy,MPSE
#' @exportMethod mp_extract_taxonomy
setMethod("mp_extract_taxonomy", signature(x="MPSE"), function(x, ...){
    da <- .internal_extract_taxonomy(taxatree=x@taxatree, classnm=class(x)[1])
    return(da)
})

#' @rdname mp_extract_taxonomy-methods
#' @aliases mp_extract_taxonomy,tbl_mpse
#' @exportMethod mp_extract_taxonomy
setMethod("mp_extract_taxonomy", signature(x="tbl_mpse"), function(x, ...){
    taxatree <- x %>% attr("taxatree")
    classnm <- class(x)[1]
    da <- .internal_extract_taxonomy(taxatree = taxatree, 
                                     classnm = classnm)
    return(da)
})

#' @rdname mp_extract_taxonomy-methods
#' @aliases mp_extract_taxonomy,grouped_df_mpse
#' @exportMethod mp_extract_taxonomy
setMethod("mp_extract_taxonomy", signature(x="grouped_df_mpse"), function(x, ...){
    taxatree <- x %>% attr("taxatree")
    classnm <- class(x)[1]
    da <- .internal_extract_taxonomy(taxatree = taxatree,
                                     classnm = classnm)
    return (da)
})

#' @title extract the sample information in MPSE object
#' @docType methods
#' @rdname mp_extract_sample-methods
#' @param x MPSE object
#' @param ... additional arguments
#' @return tbl_df contained sample information.
#' @export
setGeneric("mp_extract_sample", function(x, ...)standardGeneric("mp_extract_sample"))

#' @rdname mp_extract_sample-methods
#' @aliases mp_extract_sample,MPSE
#' @exportMethod mp_extract_sample
setMethod("mp_extract_sample", signature(x="MPSE"), function(x, ...){
     da <- x@colData %>%
        avoid_conflict_names() %>%
        tibble::as_tibble(rownames="Sample")
     return(da)
})

.internal_extract_sample <- function(x, ...){
    samplevar <- x %>% attr("samplevar")
    da <- x %>%
        dplyr::ungroup() %>%
        select(samplevar) %>%
        distinct()
    return(da)
}

#' @rdname mp_extract_sample-methods
#' @aliases mp_extract_sample,tbl_mpse
#' @exportMethod mp_extract_sample
setMethod("mp_extract_sample", signature(x="tbl_mpse"), .internal_extract_sample)

#' @rdname mp_extract_sample-methods
#' @aliases mp_extract_sample,grouped_df_mpse
#' @exportMethod mp_extract_sample
setMethod("mp_extract_sample", signature(x="grouped_df_mpse"), .internal_extract_sample)


#' @title extract the feature (OTU) information in MPSE object
#' @docType methods
#' @rdname mp_extract_feature-methods
#' @param x MPSE object
#' @param addtaxa logical whether adding the taxonomy information
#' default is FALSE.
#' @param ... additional arguments
#' @return tbl_df contained feature (OTU) information.
#' @export
setGeneric("mp_extract_feature", function(x, addtaxa=FALSE, ...)standardGeneric("mp_extract_feature"))

#' @rdname mp_extract_feature-methods
#' @aliases mp_extract_feature,MPSE
#' @exportMethod mp_extract_feature
setMethod("mp_extract_feature", signature(x="MPSE"), function(x, addtaxa=FALSE, ...){
    da <- SummarizedExperiment::rowData(x) %>%
          avoid_conflict_names() %>%
          tibble::as_tibble(rownames="OTU")
    if (!is.null(x@taxatree)){
        taxanm <- x@taxatree@data %>%
                  dplyr::filter(!.data$nodeClass %in% c("OTU", "Root")) %>%
                  dplyr::pull(.data$nodeClass) %>% unique
        trda <- x@taxatree %>%
                taxatree_to_tb() %>%
                tibble::as_tibble(rownames="OTU")
        if (!addtaxa){
            trda %<>% dplyr::select(-taxanm)
        }
        da %<>% dplyr::left_join(trda, by="OTU") 
                
    }
    return(da)
})

.internal_extract_feature <- function(x, addtaxa=FALSE, ...){
    otumetavar <- x %>% attr("otumetavar")
    taxatree <- x %>% attr("taxatree")
    da <- x %>%
        dplyr::ungroup() %>%
        select(c("OTU", otumetavar)) %>%
        distinct()
    if (!is.null(taxatree)){
        taxanm <- taxatree@data %>%
                  dplyr::filter(!.data$nodeClass %in% c("OTU", "Root")) %>%
                  dplyr::pull(.data$nodeClass) %>% unique
        trda <- taxatree %>%
                taxatree_to_tb() %>%
                tibble::as_tibble(rownames="OTU")
        if (!addtaxa){
            trda %<>% dplyr::select(-taxanm)
        }
        da %<>% dplyr::left_join(trda, by="OTU")
    }
    return(da)
}

#' @rdname mp_extract_feature-methods
#' @aliases mp_extract_feature,tbl_mpse
#' @exportMethod mp_extract_feature
setMethod("mp_extract_feature", signature(x="tbl_mpse"), .internal_extract_feature)

#' @rdname mp_extract_feature-methods
#' @aliases mp_extract_feature,grouped_df_mpse
#' @exportMethod mp_extract_feature
setMethod("mp_extract_feature", signature(x="grouped_df_mpse"), .internal_extract_feature)


.internal_extract_taxonomy <- function(taxatree, classnm){
    if (is.null(taxatree)){
        message(paste0("The taxonomy annotation is empty in the ", classnm," object"))
        return(NULL)
    }
    taxanm <- taxatree@data %>%
              dplyr::filter(.data$nodeClass != "Root") %>%
              dplyr::pull(.data$nodeClass) %>% unique
    taxada <- taxatree_to_tb(taxatree) %>%
              tibble::as_tibble(rownames="OTU") %>%
              dplyr::select(taxanm) %>%
              tibble::column_to_rownames(var="OTU")
    return(taxada)

}

#' @title extract the dist object from MPSE or tbl_mpse object
#' @docType methods
#' @rdname mp_extract_dist-methods
#' @param x MPSE object or tbl_mpse object
#' @param distmethod character the method of calculated distance.
#' @return dist object.
#' @export
setGeneric("mp_extract_dist", function(x, distmethod)standardGeneric("mp_extract_dist"))

#' @rdname mp_extract_dist-methods
#' @aliases mp_extract_dist,MPSE
#' @exportMethod mp_extract_dist
setMethod("mp_extract_dist", signature(x="MPSE"), function(x, distmethod){
    data <- x@colData %>% avoid_conflict_names() %>% as_tibble(rownames="Sample")
    res <- .internal_extract_dist(data=data, distmethod=distmethod)
    return(res)
})

#' @rdname mp_extract_dist-methods
#' @aliases mp_extract_dist,tbl_mpse
#' @exportMethod mp_extract_dist
setMethod("mp_extract_dist", signature(x="tbl_mpse"), function(x, distmethod){
    res <- .internal_extract_dist(data=x, distmethod=distmethod)
    return(res)
})

#' @rdname mp_extract_dist-methods
#' @aliases mp_extract_dist,grouped_df_mpse
#' @exportMethod mp_extract_dist
setMethod("mp_extract_dist", signature(x="grouped_df_mpse"), function(x, distmethod){
    data <- x %>% ungroup()
    res <- .internal_extract_dist(data=data, distmethod=distmethod)
    return(res)
})


#' Extracting the PCA, PCoA, etc results from MPSE or tbl_mpse object
#' @rdname mp_extract_internal_attr-methods
#' @param x MPSE or tbl_mpse object
#' @param name character 'PCA' or 'PCoA'
#' @param ... additional parameters
#' @return prcomp or pcoa etc object
#' @export
setGeneric("mp_extract_internal_attr", function(x, name, ...)standardGeneric("mp_extract_internal_attr"))

.internal_extract_internal_attr <- function(x, name, ...){
    dat <- x %>% attr("internal_attr")
    message(paste0("The object contained internal attribute: ",paste0(names(dat), collapse=" ")))
    return(dat[[name]])
}

#' @rdname mp_extract_internal_attr-methods
#' @aliases mp_extract_internal_attr,MPSE
#' @exportMethod mp_extract_internal_attr
setMethod("mp_extract_internal_attr", signature(x="MPSE"), .internal_extract_internal_attr)

#' @rdname mp_extract_internal_attr-methods
#' @aliases mp_extract_internal_attr,tbl_mpse
#' @exportMethod mp_extract_internal_attr
setMethod("mp_extract_internal_attr", signature(x="tbl_mpse"), .internal_extract_internal_attr)

#' @rdname mp_extract_internal_attr-methods
#' @aliases mp_extract_internal_attr,grouped_df_mpse
#' @exportMethod mp_extract_internal_attr
setMethod("mp_extract_internal_attr", signature(x="grouped_df_mpse"), .internal_extract_internal_attr)

.internal_extract_dist <- function(data, distmethod){
    if (!distmethod %in% colnames(data)){
        rlang::abort(paste0("There is not ", distmethod, 
                            " distance in the object, please check whether the mp_cal_dist has been performed!"))
    }
    distobj <- data %>%
            select(c("Sample", distmethod)) %>%
            distinct() %>%
            tidyr::unnest() %>%
            suppressWarnings() %>%
            rename(x="Sample", y=paste0(distmethod,"Sampley"), r=distmethod) %>%
            corrr::retract() %>%
            tibble::column_to_rownames(var=colnames(.)[1])
    distobj <- distobj[colnames(distobj), ] %>% 
               stats::as.dist() %>%
               add_attr(distmethod, "method")
    return(distobj)
}

.internal_tree <- function(x, type){
    type %<>% match.arg(c("taxatree", "otutree"))
    tree <- x %>% attr(type)
    if (!is.null(tree)){
        return(tree)
    }else{
        message(tree_empty(type=type))
    }
}

tree_empty <- function(type){
    x <- paste0("The ", type," is empty in the MPSE object!")
    return(x)
}

#' @rdname MPSE-accessors
#' @aliases rownames<-,MPSE
#' @export
setReplaceMethod("rownames", signature(x="MPSE"), function(x, value){
    nx <- methods::callNextMethod()
    oldnm <- rownames(x)
    if (!is.null(x@otutree)){
        nx@otutree <- rename_tiplab(x@otutree, oldname=oldnm, newname=value)
    }
    if (!is.null(x@taxatree)){
        nx@taxatree <- rename_tiplab(x@taxatree, oldname=oldnm, newname=value)
    }
    if (!is.null(x@refseq)){
        if (is.null(value)){
            nx@refseq <- NULL
        }else{
            names(nx@refseq) <- value[match(names(x@refseq), oldnm)]
        }
    }

    if (!is.null(value)){
        old2new <- data.frame(new=value, oldrowname=oldnm) %>% 
                   column_to_rownames(var="new")
        SummarizedExperiment::rowData(nx) <- old2new
    }
    methods::validObject(nx)
    return(nx)
})


rename_tiplab <- function(treedata, oldname, newname){
    if (is.null(newname)){
        return(NULL)
    }
    tip.label <- treedata@phylo$tip.label
    treedata@phylo$tip.label <- newname[match(tip.label, oldname)]
    return(treedata)
}

.internal_drop.tip <- function(tree, newnm=NULL, collapse.singles=TRUE, rmotus=NULL){
    if (is.null(tree)){
        return (NULL)
    }
    if (is.null(rmotus)){
        rmotus <- setdiff(tree@phylo$tip.label, newnm)
    }
    if (length(rmotus) > 0 && length(rmotus) != treeio::Ntip(tree)){
        otutree <- treeio::drop.tip(tree, tip=rmotus, collapse.singles=collapse.singles)
    }else{
        otutree <- tree
    }
    return (otutree)
}

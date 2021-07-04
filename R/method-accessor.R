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

#' @title extract the taxonomy tree in MPSE object
#' @docType methods
#' @rdname mp_taxatree-methods
#' @param x MPSE object
#' @param ... additional arguments
#' @return taxatree treedata object
#' @export
setGeneric("mp_taxatree", function(x, ...){standardGeneric("mp_taxatree")})


#' @rdname mp_taxatree-methods
#' @aliases mp_taxatree,MPSE
#' @exportMethod mp_taxatree
setMethod("mp_taxatree", signature(x="MPSE"),
          function(x, ...){
    if (!is.null(x@taxatree)){
        return(x@taxatree)
    }else{
        message(taxatree_empty)
    }
})

#' @rdname mp_taxatree-methods
#' @aliases mp_taxatree,tbl_mpse
#' @exportMethod mp_taxatree
setMethod("mp_taxatree", signature(x="tbl_mpse"),function(x, ...){
    .internal_taxatree(x=x)
})

#' @rdname mp_taxatree-methods
#' @aliases mp_taxatree,grouped_df_mpse
#' @exportMethod mp_taxatree
setMethod("mp_taxatree", signature(x="grouped_df_mpse"),function(x, ...){
    .internal_taxatree(x=x)
})

.internal_taxatree <- function(x){
    taxatree <- attr(x, "taxatree")
    if (!is.null(taxatree)){
        return(taxatree)
    }else{
        message(taxatree_empty)
    } 
}

taxatree_empty <- "The taxatree is empty in the MPSE object!"

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

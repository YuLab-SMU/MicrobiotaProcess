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
#' @param x R object, MPSE class in here
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
#' @param object parameter of tax_table, R object, MPSE class in here.
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
    rowda <- SummarizedExperiment::rowData(nx)

    otutree <- .internal_drop.tip(tree=otutree, newnm=newotus)
 
    if (length(newotus) == 1){
        if (!is.null(taxatree)){
            taxatb <- taxatree %>% 
                      taxatree_to_tb() %>% 
                      tibble::as_tibble(rownames="OTU") %>%
                      dplyr::filter(.data$OTU==newotus)
            rowda <- rowda %>% 
                     tibble::as_tibble(rownames="OTU") %>% 
                     dplyr::filter(.data$OTU==newotus) %>%
                     dplyr::left_join(taxatb, by="OTU") %>%
                     tibble::column_to_rownames(var="OTU") %>%
                     S4Vectors::DataFrame(check.names=FALSE)
            taxatree <- NULL
        }
    }else{
        taxatree <- .internal_drop.tip(
                                       tree=taxatree, 
                                       newnm=newotus, 
                                       collapse.singles=FALSE
                                )
    }

    if (!is.null(refseq)){
        refseq <- refseq[newotus]
    }
 
    nx@taxatree <- taxatree
    nx@otutree <- otutree
    nx@refseq <- refseq
    SummarizedExperiment::rowData(nx) <- rowda
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

#' @rdname MPSE-accessors
#' @export
setReplaceMethod("colData", c("MPSE", "DataFrame"), function(x, ..., value){
    res <- methods::callNextMethod()
})

#' @rdname MPSE-accessors
#' @export
setReplaceMethod("colData", c("MPSE", "NULL"), function(x, ..., value){
    res <- methods::callNextMethod()
})

#' extract the abundance matrix from MPSE object or tbl_mpse object
#' @rdname mp_extract_assays-methods
#' @param x MPSE or tbl_mpse object
#' @param .abundance the name of abundance to be extracted.
#' @param byRow logical if it is set TRUE, 'otu X sample' shape will return,
#' else 'sample X otu' will return.
#' @return otu abundance a data.frame object
#' @param ... additional parameters.
#' @export
setGeneric("mp_extract_assays", function(x, .abundance, byRow=TRUE, ...)standardGeneric("mp_extract_assays"))

#' @rdname mp_extract_assays-methods
#' @aliases mp_extract_assays,MPSE
#' @exportMethod mp_extract_assays
setMethod("mp_extract_assays", signature(x="MPSE"), function(x, .abundance, byRow=TRUE, ...){
    .abundance <- rlang::enquo(.abundance)
    if (rlang::quo_is_missing(.abundance)){
        rlang::abort("The abundance name is required !")
    }

    assaysvar <- SummarizedExperiment::assayNames(x)
    if (!rlang::as_name(.abundance) %in% assaysvar){
        message_wrap(paste0("The assays slot of object does not contain ", 
                            rlang::as_name(.abundance)))
        return(NULL)
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

#' @rdname mp_extract_assays-methods
#' @aliases mp_extract_assays,tbl_mpse
#' @exportMethod mp_extract_assays
setMethod("mp_extract_assays", signature(x="tbl_mpse"), .internal_extract_abundance)

#' @rdname mp_extract_assays-methods
#' @aliases mp_extract_assays,grouped_df_mpse
#' @exportMethod mp_extract_assays
setMethod("mp_extract_assays", signature(x="grouped_df_mpse"), .internal_extract_abundance)


#' @title extract the taxonomy tree in MPSE object
#' @docType methods
#' @rdname mp_extract_tree-methods
#' @param x MPSE object
#' @param type character taxatree or otutree
#' @param tip.level character This argument will keep the nodes 
#' belong to the tip.level as tip nodes when type is taxatree, default is OTU, 
#' which will return the taxa tree with OTU level as tips.
#' @param ... additional arguments
#' @return taxatree treedata object
#' @export
setGeneric("mp_extract_tree", function(x, type="taxatree", tip.level="OTU", ...){standardGeneric("mp_extract_tree")})


#' @rdname mp_extract_tree-methods
#' @aliases mp_extract_tree,MPSE
#' @exportMethod mp_extract_tree
setMethod("mp_extract_tree", signature(x="MPSE"), function(x, type="taxatree", tip.level="OTU", ...){
    type %<>% match.arg(c("taxatree", "otutree"))
    tree <- methods::slot(x, type)
    tip.level <- rlang::enquo(tip.level)
    if (!is.null(tree)){
        if (type == "taxatree"){
            tree <- .extract_tree_at_tiplevel(tree, tip.level=!!tip.level)
        }
        return(tree)
    }else{
        message(tree_empty(type=type))
    }
})

#' @rdname mp_extract_tree-methods
#' @aliases mp_extract_tree,tbl_mpse
#' @exportMethod mp_extract_tree
setMethod("mp_extract_tree", signature(x="tbl_mpse"),function(x, type="taxatree", tip.level="OTU", ...){
    tip.level <- rlang::enquo(tip.level)
    .internal_tree(x = x, type = type, tip.level = !!tip.level)
})

#' @rdname mp_extract_tree-methods
#' @aliases mp_extract_tree,grouped_df_mpse
#' @exportMethod mp_extract_tree
setMethod("mp_extract_tree", signature(x="grouped_df_mpse"),function(x, type="taxatree", tip.level="OTU", ...){
    tip.level <- rlang::enquo(tip.level)
    .internal_tree(x=x, type=type, tip.level=!!tip.level)
})

#' @rdname mp_extract_tree-methods
#' @export
mp_extract_taxatree <- function(x, tip.level = "OTU", ...){
    tip.level <- rlang::enquo(tip.level)
    x <- mp_extract_tree(x = x, type="taxatree", tip.level = !!tip.level, ...)
    return(x)
}

#' @rdname mp_extract_tree-methods
#' @export
mp_extract_otutree <- function(x, ...){
    x <- mp_extract_tree(x = x, type = 'otutree', ...)
    return(x)
}

#' @rdname mp_extract_taxonomy-methods
#' @param x MPSE object
#' @param ... additional parameters, now is meaningless.
#' @return data.frame contained taxonomy information
#' @export
setGeneric("taxonomy", function(x, ...){standardGeneric("taxonomy")})

.internal_taxonomy <- function(x, ...){
    x <- x %>% mp_extract_taxonomy(...) %>%
         tibble::column_to_rownames(var="OTU")
    return(x)
}

#' @rdname mp_extract_taxonomy-methods
#' @aliases taxonomy,MPSE
#' @exportMethod taxonomy
setMethod("taxonomy", signature(x = "MPSE"), .internal_taxonomy)

#' @rdname mp_extract_taxonomy-methods
#' @aliases taxonomy,tbl_mpse
#' @exportMethod taxonomy
setMethod("taxonomy", signature(x = "tbl_mpse"), .internal_taxonomy)

#' @rdname mp_extract_taxonomy-methods
#' @aliases taxonomy,grouped_df_mpse
#' @export
setMethod("taxonomy", signature(x = "grouped_df_mpse"), .internal_taxonomy)



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

.internal_extract_taxonomy_ <- function(x, ...){
    taxatree <- x %>% attr("taxatree")
    classnm <- class(x)[1]
    da <- .internal_extract_taxonomy(taxatree = taxatree,
                                     classnm = classnm)
    return(da)
}

#' @rdname mp_extract_taxonomy-methods
#' @aliases mp_extract_taxonomy,tbl_mpse
#' @exportMethod mp_extract_taxonomy
setMethod("mp_extract_taxonomy", signature(x="tbl_mpse"), .internal_extract_taxonomy_)

#' @rdname mp_extract_taxonomy-methods
#' @aliases mp_extract_taxonomy,grouped_df_mpse
#' @exportMethod mp_extract_taxonomy
setMethod("mp_extract_taxonomy", signature(x="grouped_df_mpse"), .internal_extract_taxonomy_)

#' @rdname MPSE-accessors
#' @return taxonomyTable class
#' @export
setGeneric("tax_table", function(object)standardGeneric("tax_table"))

.internal_tax_table <- function(object){
    da <- mp_extract_taxonomy(object) %>%
          tibble::column_to_rownames(var="OTU") %>% 
          as.matrix()
    phyloseq::tax_table(da)
}

#' @rdname MPSE-accessors
#' @aliases tax_table,MPSE
#' @exportMethod tax_table
setMethod("tax_table", signature(object = "MPSE"), .internal_tax_table)

#' @rdname MPSE-accessors
#' @aliases tax_table,tbl_mpse
#' @exportMethod tax_table
setMethod("tax_table", signature(object = "tbl_mpse"), .internal_tax_table)

#' @rdname MPSE-accessors
#' @aliases tax_table,grouped_df_mpse 
#' @exportMethod tax_table
setMethod("tax_table", signature(object = "grouped_df_mpse"), .internal_tax_table)

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
       data.frame(check.names=FALSE) %>%
       avoid_conflict_names() %>%
       tibble::as_tibble(rownames="Sample") %>%
       modify_AsIs_list()
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

modify_AsIs_list <- function(x, ...){
nms <- lapply(x, function(x)inherits(x, "AsIs") && typeof(x)=="list")
nms <- names(nms[unlist(nms)])
x %<>% dplyr::mutate_at(dplyr::vars(nms), ~unclass(.))
x
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
        tip.level <- x@taxatree %>% 
                     dplyr::filter(.data$isTip, keep.td=FALSE) %>% 
                     pull(.data$nodeClass) %>% 
                     unique()
        taxanm <- x@taxatree %>%
                  dplyr::filter(!.data$nodeClass %in% c(tip.level, "Root"), keep.td=FALSE) %>%
                  dplyr::pull(.data$nodeClass) %>% unique()
        trda <- x@taxatree %>%
                taxatree_to_tb() #%>%
                #tibble::as_tibble(rownames=tip.level)
        if (!addtaxa){
            trda %<>% dplyr::select(-taxanm)
        }
        trda %<>% dplyr::select(setdiff(colnames(trda), colnames(da))) %>%
            tibble::as_tibble(rownames=tip.level)
        da %<>% dplyr::left_join(trda, by=c("OTU"=tip.level)) 
                
    }
    da %<>% modify_AsIs_list()
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
        tip.level <- taxatree %>%
                     dplyr::filter(.data$isTip, keep.td=FALSE) %>%
                     pull(.data$nodeClass) %>%
                     unique()        
        taxanm <- taxatree %>%
                  dplyr::filter(!.data$nodeClass %in% c(tip.level, "Root"), 
                                keep.td = FALSE) %>%
                  dplyr::pull(.data$nodeClass) %>% unique
        trda <- taxatree %>%
                taxatree_to_tb() %>%
                tibble::as_tibble(rownames = tip.level)
        if (!addtaxa){
            trda %<>% dplyr::select(-taxanm)
        }
        da %<>% dplyr::left_join(trda, by=c("OTU"= tip.level))
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

#' Extract the result of mp_cal_rarecurve with action="add" from MPSE or tbl_mpse object
#' @rdname mp_extract_rarecurve-methods
#' @param x MPSE object or tbl_mpse object
#' @param .rarecurve the column name of rarecurve after run mp_cal_rarecurve with action="add".
#' @param ... additional parameter
#' @return rarecurve object that be be visualized by ggrarecurve
#' @export
setGeneric("mp_extract_rarecurve", function(x, .rarecurve, ...)standardGeneric("mp_extract_rarecurve"))

.internal_extract_rarecurve <- function(x, .rarecurve, ...){
.rarecurve <- rlang::enquo(.rarecurve)
if (rlang::quo_is_missing(.rarecurve)){
    .rarecurve <- as.symbol("RareAbundanceRarecurve")
}
dat <- x %>% 
    mp_extract_sample() %>%
    dplyr::select("Sample", !!.rarecurve) %>% 
    dplyr::rename(sample="Sample") %>%
    tidyr::unnest() %>% 
    suppressWarnings()
return(structure(list(data=dat), class="rarecurve"))
}

#' @rdname mp_extract_rarecurve-methods
#' @aliases mp_extract_rarecurve,MPSE
#' @exportMethod mp_extract_rarecurve
setMethod("mp_extract_rarecurve", signature(x="MPSE"), .internal_extract_rarecurve)

#' @rdname mp_extract_rarecurve-methods
#' @aliases mp_extract_rarecurve,tbl_mpse
#' @exportMethod mp_extract_rarecurve
setMethod("mp_extract_rarecurve", signature(x="tbl_mpse"), .internal_extract_rarecurve)

#' @rdname mp_extract_rarecurve-methods
#' @aliases mp_extract_rarecurve,grouped_df_mpse
#' @exportMethod mp_extract_rarecurve
setMethod("mp_extract_rarecurve", signature(x="grouped_df_mpse"), .internal_extract_rarecurve)

#' Extracting the abundance metric from MPSE or tbl_mpse object
#' @description 
#' Extracting the abundance metric from the MPSE or tbl_mpse,
#' the 'mp_cal_abundance' must have been run with action='add'.
#' @rdname mp_extract_abundance-methods
#' @param x MPSE or tbl_mpse object
#' @param taxa.class character the name of taxonomy 
#' class level what you want to extract
#' @param topn integer the number of the top most abundant, default
#' is NULL.
#' @param rmun logical whether to remove the unknown taxa, such as "g__un_xxx",
#' default is FALSE (the unknown taxa class will be considered as 'Others').
#' @param ... additional parameters
#' @author Shuangbin Xu
#' @export
setGeneric("mp_extract_abundance", function(x, taxa.class="all", topn=NULL, rmun=FALSE, ...)standardGeneric("mp_extract_abundance"))

#' @importFrom dplyr all_of
.internal_extract_abundance <- function(x, taxa.class="all", topn = NULL, rmun = FALSE, ...){
    taxa.class <- rlang::enquo(taxa.class)
    
    taxatree <-  x %>% 
                 mp_extract_tree()
    if (inherits(x, "MPSE")){
        assaysvar <- x %>% SummarizedExperiment::assayNames()
    }else{
        assaysvar <- x %>% attr("assaysvar")
    }
    
    if (is.null(taxatree)){
        taxa.class <- rlang::sym("OTU")
        da <- x %>% mp_extract_feature() %>% 
              dplyr::rename(label="OTU") %>% 
              dplyr::mutate(nodeClass="OTU") 
        if (ncol(da)==1){
            message_wrap("Please make sure the mp_cal_abundance(..., action='add') has been run.
                          Or you can extract the assay via mp_extract_assays since the taxonomy is NULL")
            return(NULL)
        }
        #return(da)
    }else{
        flag <-c(colnames(taxatree@data), colnames(taxatree@extraInfo)) %>% 
               unique() %in% c("node", "nodeClass", "nodeDepth")
        if (length(flag)==3 && all(flag)){
            message("Please make sure the mp_cal_abundance(..., action='add') has been run.")
            return(NULL)        
        }
        da <- taxatree %>% 
              as_tibble %>%    
              dplyr::select(-c("parent", "node", "nodeDepth")) %>%
              dplyr::filter(.data$nodeClass != "Root")
        taxa.class <- rlang::as_name(taxa.class)
        if (taxa.class!="all"){
            da <- da %>% 
                  dplyr::filter(.data$nodeClass == taxa.class)
        }
    }
    
    if (taxa.class!="all"){
        AbundBy <- colnames(da)[vapply(da, is.list, logical(1))]
        dat <- da %>% tidyr::unnest(cols=AbundBy[1])
        clnm <- colnames(dat)[vapply(dat, is.numeric, logical(1))]
        if (rmun){
            dat %<>% dplyr::mutate(label=ifelse(grepl("__un_", .data$label), "Others", .data$label))
        }
        totallabel <- dat %>%
              dplyr::group_by(.data$label) %>%
              dplyr::summarize(TotalByLabel=sum(!!as.symbol(clnm[1]))) %>%
              dplyr::arrange(dplyr::desc(.data$TotalByLabel)) %>% 
              dplyr::filter(.data$label != "Others") %>%
              dplyr::pull(.data$label)
        if (is.null(topn)){topn <- length(totallabel)}
        keepn <- min(topn, length(totallabel))
        if (keepn < length(totallabel)){
            keeplabel <- totallabel[seq_len(keepn)]
        }else{
            keeplabel <- totallabel
        }
        dtf <- list()
        for (i in AbundBy){
            df <- da %>% 
                  select(c("label", "nodeClass", i)) %>%
                  tidyr::unnest(cols=i)
            nms <- colnames(df)
            abunnm <- nms[vapply(df, is.numeric, logical(1))]
            gpnm <- nms[!nms %in% c("label", "nodeClass", abunnm)]
            if (gpnm[1]=="Sample"){
                gpnm <- "Sample"
            }else{
                gpnm <- gpnm
            }
            df1 <- df %>% 
                   dplyr::filter(.data$label %in% keeplabel)
            df2 <- df %>%
                   dplyr::filter(!.data$label %in% keeplabel) %>%
                   dplyr::mutate(label="Others") %>% 
                   dplyr::group_by(across(all_of(gpnm))) %>%
                   dplyr::mutate(across(all_of(abunnm), sum)) %>%
                   dplyr::ungroup() %>%
                   distinct() 
            dtf[[i]] <- dplyr::bind_rows(df1, df2) %>%
                  dplyr::mutate(label=factor(.data$label, levels=c(keeplabel, "Others"))) %>%
                  tidyr::nest(!!i:=nms[!nms %in% c("label", "nodeClass")])
    
        }
        da <- dtf %>% purrr::reduce(left_join, by=c("label", "nodeClass"))
    }
    return(da)
}

#' @rdname mp_extract_abundance-methods
#' @aliases mp_extract_abundance,MPSE
#' @export mp_extract_abundance
setMethod("mp_extract_abundance", signature(x="MPSE"), .internal_extract_abundance)

#' @rdname mp_extract_abundance-methods
#' @aliases mp_extract_abundance,tbl_mpse
#' @export mp_extract_abundance
setMethod("mp_extract_abundance", signature(x="tbl_mpse"), .internal_extract_abundance)

#' @rdname mp_extract_abundance-methods
#' @aliases mp_extract_abundance,grouped_df_mpse
#' @export mp_extract_abundance
setMethod("mp_extract_abundance", signature(x="grouped_df_mpse"), .internal_extract_abundance)


.internal_extract_taxonomy <- function(taxatree, classnm){
    if (is.null(taxatree)){
        message(paste0("The taxonomy annotation is empty in the ", classnm," object"))
        return(NULL)
    }
    tip.level <- taxatree %>% 
                 dplyr::filter(.data$isTip, keep.td = FALSE) %>% 
                 pull(.data$nodeClass) %>% 
                 unique()
    taxanm <- taxatree@data %>%
              dplyr::filter(!.data$nodeClass %in% c(tip.level, "Root")) %>%
              dplyr::pull(.data$nodeClass) %>% unique
    taxada <- taxatree_to_tb(taxatree) %>%
              tibble::as_tibble(rownames = "OTU") %>%
              dplyr::select(c("OTU", taxanm)) #%>%
              #tibble::column_to_rownames(var="OTU")
    return(taxada)

}

#' Extracting the PCA, PCoA, etc results from MPSE or tbl_mpse object
#' @rdname mp_extract_internal_attr-methods
#' @param x MPSE or tbl_mpse object
#' @param name character 'PCA' or 'PCoA'
#' @param ... additional parameters
#' @return prcomp or pcoa etc object
#' @export
setGeneric("mp_extract_internal_attr", function(x, name, ...)standardGeneric("mp_extract_internal_attr"))

.internal_extract_internal_attr <- function(x, name, ...){
name <- rlang::enquo(name) %>% rlang::as_name()
dat <- x %>% attr("internal_attr")
message(paste0("The object contained internal attribute: ",paste0(names(dat), collapse=" ")))
indx <- grep(paste0(name, "$"), names(dat), ignore.case=TRUE)
return(dat[[indx]])
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

#' @title extract the dist object from MPSE or tbl_mpse object
#' @docType methods
#' @rdname mp_extract_dist-methods
#' @param x MPSE object or tbl_mpse object
#' @param distmethod character the method of calculated distance.
#' @param type character, which type distance to be extracted, 'sample' represents
#' the distance between the samples based on feature abundance matrix, 'feature' represents 
#' the distance between the features based on feature abundance matrix, 'env' represents the
#' the distance between the samples based on continuous environment factors, default is 'sample'.
#' @param .group the column name of sample information, which only work with type='sample' or 
#' type='env', default is NULL, when it is provided, a tibble that can be visualized via ggplot2 
#' will return.
#' @param ... additional parameters
#' @return dist object or tbl_df object when .group is provided.
#' @export
setGeneric("mp_extract_dist", function(x, distmethod, type='sample', .group=NULL, ...)standardGeneric("mp_extract_dist"))

.internal_extract_dist <- function(x, distmethod, type='sample', .group=NULL, ...){
    .group <- rlang::enquo(.group)
    dots <- list(...)
    type %<>% match.arg(c('sample', 'feature', 'env'))

    if ('env.flag' %in% names(dots)){
        if (dots$env.flag){
            type <- 'env'
        }else{
            type <- 'sample'
        }
    }

    if(type == 'feature'){
        data <- x %>% mp_extract_feature(addtaxa = TRUE)
        distname <- paste0(distmethod, 'Featurey') %>% as.symbol()
        distmethod <- paste0('Feature_', distmethod)
        .group <- rlang::quo(NULL)
        prefix <- 'OTU'
    }else{
        data <- x %>% mp_extract_sample()
        if (type == 'env'){
           distmethod <- paste0('Env_', distmethod) 
        }
        distname <- paste0(distmethod, "Sampley") %>% as.symbol()
        prefix <- 'Sample'
    }
           
    #data <- x %>% mp_extract_sample()
    #distmethod <- switch(as.character(type),
    #                     "TRUE" = paste0("Env_", distmethod),
    #                     "FALSE" = distmethod)
    
    if (!distmethod %in% colnames(data)){
        rlang::abort(paste0("There is not ", distmethod, 
                            " distance in the object, please check whether the mp_cal_dist has been performed!"))
    }
    
    #distname <- paste0(distmethod, "Sampley") %>% as.symbol()
    
    if (rlang::quo_is_null(.group)){
        distobj <- data %>%
                select(c(prefix, distmethod)) %>%
                distinct() %>%
                tidyr::unnest() %>%
                suppressWarnings() %>%
                rename(x=prefix, y=distname, r=distmethod) %>%
                corrr::retract() %>%
                tibble::column_to_rownames(var=colnames(.)[1]) %>%
                magrittr::extract(,rownames(.))
        #distobj <- distobj[colnames(distobj), ] 
        distobj[lower.tri(distobj)] <- t(distobj)[lower.tri(t(distobj))]
        distobj %<>% stats::as.dist() %>%
                     add_attr(distmethod, "method")
        return(distobj)
    }else{
        group.y <- paste0(rlang::as_name(.group), ".tmp") %>% as.symbol()
        dist.tb <- data %>%
                   dplyr::select(c(prefix, distmethod, !!.group)) %>%
                   tidyr::unnest(cols=distmethod) %>%
                   dplyr::mutate(dplyr::across(!!.group, 
                                               ~.x[match(!!as.symbol(distname), !!as.symbol(prefix))], 
                                               .names=rlang::as_name(group.y))) %>% 
                   dplyr::rowwise() %>% 
                   dplyr::mutate(GroupsComparison=paste0(sort(c(!!.group, !!group.y)),collapse="-vs-")) %>%
                   dplyr::filter(!!rlang::sym(prefix) != !!as.symbol(distname)) %>%
                   dplyr::select(c(prefix, distmethod, "GroupsComparison", distname))
        return(dist.tb)
    }
}

#' @rdname mp_extract_dist-methods
#' @aliases mp_extract_dist,MPSE
#' @exportMethod mp_extract_dist
setMethod("mp_extract_dist", signature(x="MPSE"), .internal_extract_dist)

#' @rdname mp_extract_dist-methods
#' @aliases mp_extract_dist,tbl_mpse
#' @exportMethod mp_extract_dist
setMethod("mp_extract_dist", signature(x="tbl_mpse"), .internal_extract_dist)

#' @rdname mp_extract_dist-methods
#' @aliases mp_extract_dist,grouped_df_mpse
#' @exportMethod mp_extract_dist
setMethod("mp_extract_dist", signature(x="grouped_df_mpse"), .internal_extract_dist)


.internal_tree <- function(x, type, tip.level){
    tip.level <- rlang::enquo(tip.level)
    type %<>% match.arg(c("taxatree", "otutree"))
    tree <- x %>% attr(type)
    if (!is.null(tree)){
        if (type == "taxatree"){
            tree <- .extract_tree_at_tiplevel(tree, tip.level=!!tip.level)
        }
        return(tree)
    }else{
        message(tree_empty(type=type))
    }
}

#' @importFrom treeio drop.tip
.extract_tree_at_tiplevel <- function(tree, tip.level){
    tip.level <- rlang::enquo(tip.level) %>% rlang::as_name()
    if (tip.level == "OTU"){
        return(tree)
    }
    rmnms <- tree %>% as_tibble %>% 
             dplyr::filter(.data$nodeDepth > .data$nodeDepth[match(tip.level, .data$nodeClass)]) %>%
             select(.data$label, .data$nodeDepth) %>% 
             group_by(.data$nodeDepth) %>% 
             dplyr::summarize(label=list(.data$label)) %>% 
             arrange(dplyr::desc(.data$nodeDepth)) %>% 
             pull(.data$label) 
    if (length(rmnms) > 0){
        for ( i in rmnms){
            tree <- drop.tip(tree, tip=i, collapse.singles=FALSE, trim.internal=FALSE)
        }
    }
    return(tree)
}

tree_empty <- function(type){
    x <- paste0("The ", type," is empty in the MPSE object!")
    return(x)
}

#' @rdname MPSE-accessors
#' @param x MPSE object
#' @export
setGeneric("otutree", function(x, ...)standardGeneric("otutree"))

#' @rdname MPSE-accessors
#' @aliases otutree,MPSE
#' @export
setMethod("otutree", signature(x="MPSE"), function(x,...){
    tree <- x %>% mp_extract_otutree(...)
    return(tree)
})


#' @rdname MPSE-accessors
#' @aliases otutree,tbl_mpse
#' @export
setMethod("otutree", signature(x="tbl_mpse"), function(x,...){
    tree <- x %>% mp_extract_otutree(...)
    return(tree)
})

#' @rdname MPSE-accessors
#' @aliases otutree,group_df_mpse
#' @export
setMethod("otutree", signature(x="MPSE"), function(x,...){
    tree <- x %>% mp_extract_otutree(...)
    return(tree)
})

#' @rdname MPSE-accessors 
#' @param x MPSE object
#' @param value treedata class, phylo class or NULL
#' @export
setGeneric("otutree<-", function(x, ..., value)standardGeneric("otutree<-"))

.internal_otutree_replace <- function(x, ..., value){
    if (inherits(value, "treedata")){
        newnms <- intersect(rownames(x), value@phylo$tip.label)
    }else if(inherits(value, "phylo")){
        newnms <- intersect(rownames(x), value$tip.label)
    }
    if (length(newnms)==0){
        stop_wrap("There are not the same labels between tip labels of the tree 
                   and the rownames of mpse, please check the tree is correct")
    }
    if (length(newnms) != nrow(x)){
        dropnms <- setdiff(rownames(x), newnms)
        message_wrap("droping rows without the input tree matches:")
        message_wrap(paste0(dropnms, collapse="\n"))
        x <- x[rownames(x) %in% newnms,]
    }
    x@otutree <- .internal_drop.tip(tree=value, newnm=rownames(x), collapse.singles=FALSE) %>%
                 as.treedata() 
    methods::validObject(x)
    return(x)
}

#' @rdname MPSE-accessors
#' @aliases otutree<-,MPSE
#' @export
setReplaceMethod("otutree", signature(x="MPSE", value="treedata"), .internal_otutree_replace)

#' @rdname MPSE-accessors
#' @aliases otutree<-,MPSE
#' @export
setReplaceMethod("otutree", signature(x="MPSE", value="phylo"), .internal_otutree_replace)

#' @rdname MPSE-accessors
#' @aliases otutree<-,MPSE
#' @export
setReplaceMethod("otutree", signature(x="MPSE", value="NULL"), function(x, ..., value){
    x@otutree <- NULL
    methods::validObject(x)
    return(x)
})

.internal_replace_otutree2 <- function(x, ..., value){
    nms <- x %>% dplyr::pull(.data$OTU) %>% unique()
    value <- .internal_drop.tip(tree = value, newnm = nms, collapse.singles = FALSE)
    x %<>%
          add_attr(value, name="otutree")
    return(x)
}

#' @rdname MPSE-accessors
#' @aliases otutree<-,tbl_mpse
#' @export
setReplaceMethod("otutree", signature(x="tbl_mpse", value="treedata"), .internal_replace_otutree2)

#' @rdname MPSE-accessors
#' @aliases otutree<-,grouped_df_mpse
#' @export
setReplaceMethod("otutree", signature(x="grouped_df_mpse", value="treedata"), .internal_replace_otutree2)

#' @rdname MPSE-accessors
#' @aliases otutree<-,tbl_mpse
#' @export
setReplaceMethod("otutree", signature(x="tbl_mpse", value="NULL"), .internal_replace_otutree2)

#' @rdname MPSE-accessors
#' @aliases otutree<-,grouped_df_mpse
#' @export
setReplaceMethod("otutree", signature(x="grouped_df_mpse", value="NULL"), .internal_replace_otutree2)

#' @rdname MPSE-accessors
#' @param x MPSE object
#' @export
setGeneric("taxatree", function(x, ...)standardGeneric("taxatree"))

#' @rdname MPSE-accessors
#' @aliases taxatree,MPSE
#' @export
setMethod("taxatree", signature(x="MPSE"),function(x, ...){
    tree <- x %>% mp_extract_taxatree(...)
    return(tree)
})

#' @rdname MPSE-accessors
#' @aliases taxatree,tbl_mpse
#' @export
setMethod("taxatree", signature(x = 'tbl_mpse'), function(x, ...){
    x %>% mp_extract_taxatree(...)
})

#' @rdname MPSE-accessors
#' @aliases taxatree,grouped_df_mpse
#' @export
setMethod('taxatree', signature(x = 'grouped_df_mpse'), function(x, ...){
    x %>% mp_extract_taxatree(...)
})

#' @rdname MPSE-accessors
#' @param x MPSE object
#' @param value  treedata object or NULL
#' @export
setGeneric("taxatree<-", function(x, ..., value)standardGeneric("taxatree<-"))

#' @rdname MPSE-accessors
#' @aliases taxatree<-,MPSE
#' @export
setReplaceMethod("taxatree", signature(x="MPSE", value="treedata"), function(x, ..., value){
    newnms <- intersect(rownames(x), value@phylo$tip.label)
    if (length(newnms) != nrow(x)){
        dropnms <- setdiff(rownames(x), newnms)
        message_wrap("droping rows without the input taxonomy matches:")
        message_wrap(paste0(dropnms, collapse="\n"))
        x <- x[rownames(x) %in% newnms,]
    }
    if (length(newnms)==0){
        stop_wrap("There are not the same labels between rownames of the taxonomy information 
                   and the rownames of mpse, please check the taxonomy is correct.")
    }
    x@taxatree <- .internal_drop.tip(tree=value, newnm=rownames(x), collapse.singles=FALSE)
    methods::validObject(x)
    return(x)
})

#' @rdname MPSE-accessors
#' @aliases taxatree<-,MPSE
#' @export
setReplaceMethod("taxatree", signature(x="MPSE", value="NULL"), function(x, ..., value){
    x@taxatree <- NULL
    methods::validObject(x)
    return(x)
})

.internal_replace_taxatree <- function(x, ..., value){
    attr(x, 'taxatree') <- value
    return(x)
}

#' @rdname MPSE-accessors
#' @aliases taxatree<-,tbl_mpse
#' @export
setReplaceMethod('taxatree', signature(x = 'tbl_mpse', value = 'treedata'), .internal_replace_taxatree)

#' @rdname MPSE-accessors
#' @aliases taxatree<-,tbl_mpse
#' @export
setReplaceMethod('taxatree', signature(x = 'tbl_mpse', value = 'NULL'), .internal_replace_taxatree)

#' @rdname MPSE-accessors
#' @aliases taxatree<-,grouped_df_mpse
#' @export
setReplaceMethod('taxatree', signature(x = 'grouped_df_mpse', value = 'treedata'), .internal_replace_taxatree)

#' @rdname MPSE-accessors
#' @aliases taxatree<-,grouped_df_mpse
#' @export 
setReplaceMethod('taxatree', signature(x = 'grouped_df_mpse', value = 'NULL'), .internal_replace_taxatree)


#' @rdname MPSE-accessors
#' @param x MPSE object
#' @param value data.frame, matrix, taxonomyTable or NULL
#' @export
setGeneric("taxonomy<-", function(x, ..., value)standardGeneric("taxonomy<-"))

.internal_taxonomy_replace <- function(x, ..., value){
    if (is.null(value)){
        taxa.tree <- NULL
    }else{
        taxa.tree <- value %>% convert_to_treedata(include.rownames = TRUE, ...)
    }
    taxatree(x) <- taxa.tree
    return(x)
}

#' @rdname MPSE-accessors
#' @aliases taxonomy<-,MPSE
#' @export
setReplaceMethod("taxonomy", signature(x = "MPSE", value = "data.frame"), .internal_taxonomy_replace)

#' @rdname MPSE-accessors
#' @aliases taxonomy<-,MPSE
#' @export
setReplaceMethod("taxonomy", signature(x = "MPSE", value = "matrix"), .internal_taxonomy_replace)

#' @rdname MPSE-accessors
#' @aliases taxonomy<-,MPSE
#' @export
setReplaceMethod("taxonomy", signature(x = "MPSE", value = "taxonomyTable"), .internal_taxonomy_replace)

#' @rdname MPSE-accessors
#' @aliases taxonomy<-,MPSE
#' @export
setReplaceMethod("taxonomy", signature(x = "MPSE", value = "NULL"), .internal_taxonomy_replace)

#' @rdname MPSE-accessors
#' @param x MPSE object
#' @export
setGeneric("refsequence", function(x, ...)standardGeneric("refsequence"))

#' @rdname MPSE-accessors
#' @aliases refsequence,MPSE
#' @export
setMethod("refsequence", signature(x="MPSE"), function(x, ...){
    refseq <- x@refseq
    if (is.null(refseq)){
        message("The representative sequence is empty")
    }
    return(refseq)
})

#' @rdname MPSE-accessors
#' @param x MPSE object
#' @param value XStringSet object or NULL
#' @export
setGeneric("refsequence<-", function(x, ..., value)standardGeneric("refsequence<-"))

#' @rdname MPSE-accessors
#' @aliases refsequence<-,MPSE
#' @export
setReplaceMethod("refsequence", signature(x="MPSE", value="XStringSet"), function(x, ..., value){
    x@refseq <- value[rownames(x)]
    methods::validObject(x)
    return(x)
})

#' @rdname MPSE-accessors
#' @aliases refsequence<-,MPSE
#' @export
setReplaceMethod("refsequence", signature(x="MPSE", value="NULL"), function(x, ..., value){
    x@refseq <- NULL
    methods::validObject(x)
    return(x)
})

#' Extract the representative sequences from MPSE object
#' @param x MPSE object
#' @param ... additional parameters, meaningless now.
#' @rdname mp_extract_refseq-methods
#' @export
setGeneric("mp_extract_refseq", function(x, ...)standardGeneric("mp_extract_refseq"))

.internal_extract_refseq <- function(x, ...){
    if (inherits(x, "MPSE")){
        ref.seq <- x@refseq
    }else{
        ref.seq <- x %>% attr("refseq")
    }
    return(ref.seq)
}

#' @rdname mp_extract_refseq-methods
#' @aliases mp_extract_refseq,MPSE
#' @export mp_extract_refseq
setMethod("mp_extract_refseq", signature(x="MPSE"), .internal_extract_refseq)

#' @rdname mp_extract_refseq-methods
#' @aliases mp_extract_refseq,tbl_mpse
#' @export mp_extract_refseq
setMethod("mp_extract_refseq", signature(x="tbl_mpse"), .internal_extract_refseq)

#' @rdname mp_extract_refseq-methods
#' @aliases mp_extract_refseq,grouped_df_mpse
#' @export mp_extract_refseq
setMethod("mp_extract_refseq", signature(x="grouped_df_mpse"), .internal_extract_refseq)


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

    if (!is.null(value) && !is.null(oldnm)){
        old2new <- data.frame(.NEW=value, .OLDROWNAMES=oldnm) 
        old2new %<>% dplyr::left_join(x %>% mp_extract_feature(), by=c(.OLDROWNAMES="OTU")) %>%
                   tibble::column_to_rownames(var=".NEW")
        
        SummarizedExperiment::rowData(nx) <- old2new
    }
    methods::validObject(nx)
    return(nx)
})

#' select specific taxa level as rownames of MPSE
#' @param x MPSE object
#' @param tip.level the taxonomy level, default is 'OTU'.
#' @rdname mp_select_as_tip-methods
#' @export
#' @examples
#' \dontrun{
#' data(mouse.time.mpse)
#' newmpse <- mouse.time.mpse %>%
#'            mp_select_as_tip(tip.level = Species)
#' newmpse
#' }
setGeneric("mp_select_as_tip", function(x, tip.level = 'OTU')standardGeneric("mp_select_as_tip"))

.mp_select_tip <- function(x, tip.level="OTU"){
    tip.level <- rlang::enquo(tip.level)
    if (rlang::as_name(tip.level) == "OTU" || is.null(taxatree(x))){
        return(x)
    }
    newassay <- x %>% mp_cal_abundance(
       .abundance = "Abundance", 
       force = TRUE, 
       relative = FALSE,
       action = 'only'
    ) %>%
    suppressMessages() %>%
    dplyr::filter(.data$nodeClass == rlang::as_name(tip.level)) %>%
    tidyr::unnest(cols = "AbundanceBySample") %>%
    select("label", "Sample", "Abundance") %>%
    tidyr::pivot_wider(
       id_cols = "label", 
       names_from = "Sample", 
       values_from = "Abundance"
    ) %>%
    tibble::column_to_rownames(var = 'label')

    abund <- x %>% mp_extract_assays(.abundance="Abundance")
    if (is.null(x %>% mp_extract_assays(.abundance="RareAbundance")) %>% suppressMessages() &&
        identical(all.equal(abund, round(abund)), TRUE)
        ){
        message_wrap("
           The assays slot of original MPSE does not contain RareAbundance,
           The mp_rrarefy() will be run before selecting the specific tip level automatically.
                     ")
        x %<>% mp_rrarefy()
    }

    if ("RareAbundance" %in% SummarizedExperiment::assayNames(x)){
        newassay <- list(
          Abundance = newassay, 
          RareAbundance = x %>% mp_cal_abundance(
                .abundance = "RareAbundance",
                force = TRUE,
                relative = FALSE,
                action = "only"
            ) %>%
            suppressMessages() %>%
            dplyr::filter(.data$nodeClass == rlang::as_name(tip.level)) %>%
            tidyr::unnest(cols = "RareAbundanceBySample") %>%
            select("label", "Sample", "RareAbundance") %>%
            tidyr::pivot_wider(
               id_cols = "label",
               names_from = "Sample",
               values_from = "RareAbundance"
            ) %>%
            tibble::column_to_rownames(var = 'label') 
        )
    }

    taxa.tree <- x %>% mp_extract_taxatree(tip.level = !!tip.level)
    mpse <- MPSE(assays = newassay, taxatree = taxa.tree)
    if (inherits(x, "MPSE")){
        colData(mpse) <- SummarizedExperiment::colData(x)
    }else{
        mpse <- as_tibble(mpse) %>% 
                dplyr::left_join(x %>% mp_extract_sample, by="Sample")
    }
    return(mpse)
}

#' @rdname mp_select_as_tip-methods
#' @aliases mp_select_as_tip,MPSE
#' @export mp_select_as_tip
setMethod("mp_select_as_tip", signature(x = "MPSE"), .mp_select_tip)

#' @rdname mp_select_as_tip-methods
#' @aliases mp_select_as_tip,tbl_mpse
#' @export mp_select_as_tip
setMethod("mp_select_as_tip", signature(x = "tbl_mpse"), .mp_select_tip)

#' @rdname mp_select_as_tip-methods
#' @aliases mp_select_as_tip,grouped_df_mpse
#' @export mp_select_as_tip
setMethod("mp_select_as_tip", signature(x = "grouped_df_mpse"), .mp_select_tip)

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
        if (inherits(tree, "treedata")){
            rmotus <- setdiff(tree@phylo$tip.label, newnm)
        }else if (inherits(tree, "phylo")){
            rmotus <- setdiff(tree$tip.label, newnm)
        }
    }
    if (length(rmotus) > 0 && length(rmotus) != treeio::Ntip(tree)){
        otutree <- treeio::drop.tip(tree, tip=rmotus, collapse.singles=collapse.singles)
    }else{
        otutree <- tree
    }
    return (otutree)
}

.dist2tbl <- function(x, y){
    distmethod <- 'DistFromOutSide'
    if (inherits(x, 'list')){
        if (!is.null(names(x))){
            distmethod <- names(x)
        }
        x <- x[[1]]
    }
    x <- as.matrix(x)
    flag <- match(unique(dim(x)), dim(y))
    if (is.na(flag)){
        stop_wrap('The y is a dist object, but the dimension is different with the x (mpse).')
    }
    if (flag == 1){
        x.name <- 'OTU'
    }else{
        x.name <- 'Sample'
    }
    distsampley <- paste0(distmethod, x.name, 'y')
    x <- x %>%
        corrr::as_cordf(diagonal=0) %>%
        corrr::shave() %>%
        corrr::stretch(na.rm=TRUE) %>%
        dplyr::rename(!!distmethod:="r", !!distsampley:="y", !!x.name:='x') %>%
        tidyr::nest(!!distmethod:=c(!!as.symbol(distsampley), !!as.symbol(distmethod)))
    return(x)
}

##' @method filter diffAnalysisClass
##' @importFrom rlang quos
##' @export
filter.diffAnalysisClass <- function(.data, ..., .preserve = FALSE) {
    dots <- quos(...)
    .data@result %<>% filter(!!!dots, .preserve = .preserve)
    return(.data)
}

##' @method filter alphasample
##' @export
filter.alphasample <- function(.data, ..., .preserve = FALSE){
    .data <- internal_filter(x=.data, ..., .preserve = .preserve)
    return(.data)
}

##' @method filter phyloseq
##' @importFrom ape keep.tip
##' @export
filter.phyloseq <- function(.data, ..., .preserve = FALSE){
    .data <- internal_filter(x=.data, ..., .preserve = .preserve)
    return (.data)
}

#' @method filter tbl_ps
#' @export
filter.tbl_ps <- function(.data, ..., .preserve = FALSE){
    res <- NextMethod()
    res <- add_attr.tbl_ps(x1=res, x2=.data)
    return(res)
}

#' @method filter grouped_df_ps
#' @export
filter.grouped_df_ps <- function(.data, ..., .preserve=FALSE){
    res <- NextMethod()
    res <- add_attr.tbl_ps(x1=res, x2=.data, class="grouped_df_ps")
    return(res)
}

internal_filter <- function(x, ..., .preserve){
    dots <- quos(...)
    x <- as.data.frame(x)
    x %<>% filter(!!!dots, .preserve = .preserve)
    return(x)
}

##' @method select diffAnalysisClass 
##' @importFrom dplyr select
##' @export
select.diffAnalysisClass <- function(.data, ...) {
    dots <- quos(...)
    .data@result %<>% select(!!!dots, ...)
    return(.data)
}

##' @method select alphasample
##' @importFrom dplyr select
##' @export
select.alphasample <- function(.data, ...) {
    dots <- quos(...)
    .data <- internal_select(x=.data, dots=dots, ...)
    return(.data)
}

##' @method select phyloseq
##' @export
select.phyloseq <- function(.data, ...){
    dots <- quos(...)
    .data <- internal_select(x=.data, dots=dots, ...)
    if(any(!c("Sample", "OTU", "Abundance") %in% colnames(.data))){
        stop("The Sample,OTU,Abundance columns should not be removed !")
    }
    return (.data)
}

##' @method select tbl_ps
##' @export
select.tbl_ps <- function(.data, ...){
    loc <- tidyselect::eval_select(expr(c(...)), .data)
    loc <- loc[!names(loc) %in% c("Sample", "OTU", "Abundance")]
    .data <- check_attr.tbl_ps(x=.data, recol=loc, type="select")
    res <- NextMethod()
    res <- add_attr.tbl_ps(x1=res, x2=.data)
    return (res)
}

##' @method select grouped_df_ps
##' @export
select.grouped_df_ps <- function(.data, ...){
    loc <- tidyselect::eval_select(expr(c(...)), .data)
    loc <- loc[!names(loc) %in% c("Sample", "OTU", "Abundance")]
    .data <- check_attr.tbl_ps(x=.data, recol=loc, type="select")
    res <- NextMethod()
    res <- add_attr.tbl_ps(x1=res, x2=.data, class="grouped_df_ps")
    return (res)

}

internal_select <- function(x, dots, ...){
    x <- as.data.frame(x)
    x %<>% select(!!!dots, ...)
    return(x)
}

##' @importFrom dplyr group_by_drop_default
##' @importFrom dplyr group_by
##' @method group_by phyloseq
##' @export
group_by.phyloseq <- function(.data, ..., .add = FALSE, .drop = group_by_drop_default(.data)){
    message("A tibble is returned for independent data analysis.")
    .data %<>% as_tibble()
    res <- group_by(.data=.data, ..., .add = .add, .drop = .drop)
    res <- add_attr.tbl_ps(x1=res, x2=.data)
    return(res)
}

##' @method group_by tbl_ps
##' @export
group_by.tbl_ps <- function(.data, ..., .add = FALSE, .drop = group_by_drop_default(.data)){
    res <- NextMethod()
    res <- add_attr.tbl_ps(x1=res, x2=.data)
    class(res) <- c("grouped_df_ps", class(res))
    return(res)
}

##' @method ungroup grouped_df_ps
##' @export
ungroup.grouped_df_ps <- function(x, ...){
    res <- NextMethod()
    res <- add_attr.tbl_ps(x1=res, x2=x)
    res <- drop_class(res, class=c("grouped_df_ps", "grouped_df"))
    return(res)
}

##' @method mutate phyloseq
##' @export
mutate.phyloseq <- function(.data, ...){
    .data %<>% as_tibble()
    res <- mutate(.data=.data, ...)
    return (res)
}

##' @method mutate tbl_ps
##' @export
mutate.tbl_ps <- function(.data, ...){
    res <- NextMethod()
    res <- add_attr.tbl_ps(x1=res, x2=.data)
    res <- add_var(res, type="mutate")
    return(res)
}

##' @method mutate grouped_df_ps
##' @export
mutate.grouped_df_ps <- function(.data, ...){
    res <- NextMethod()
    res <- add_attr.tbl_ps(x1=res, x2=.data, class="grouped_df_ps")
    res <- add_var(res, type="mutate")
    return(res)
}

##' @method distinct phyloseq
##' @export
distinct.phyloseq <- function(.data, ..., .keep_all = FALSE){
    .data %<>%  as_tibble()
    res <- distinct(.data=.data, ..., .keep_all = .keep_all)
    return(res)
}

##' @method distinct tbl_ps
##' @export
distinct.tbl_ps <- function(.data, ..., .keep_all = FALSE){
    res <- NextMethod()
    res <- add_attr.tbl_ps(x1 = res, x2 = .data)
    return (res)
}

##' @method distinct grouped_df_ps
##' @export
distinct.grouped_df_ps <- function(.data, ..., .keep_all = FALSE){
    res <- NextMethod()
    res <- add_attr.tbl_ps(x1 = res, x2 = .data, class="grouped_df_ps")
    return (res)
}

##' @method rename phyloseq
##' @export
rename.phyloseq <- function(.data, ...){
    .data %<>% as_tibble()
    res <- rename(.data=.data, ...)
    return (res)
}

##' @method rename tbl_ps
##' @importFrom tidyselect eval_select
##' @importFrom rlang expr
##' @export
rename.tbl_ps <- function(.data, ...){
    cols <- tidyselect::eval_select(expr(c(...)), .data)
    .data <- check_attr.tbl_ps(x=.data, recol=cols)
    res <- NextMethod()
    res <- add_attr.tbl_ps(x1=res, x2=.data)
    return (res)
}

##' @method rename grouped_df_ps
##' @export
rename.grouped_df_ps <- function(.data, ...){
    cols <- tidyselect::eval_select(expr(c(...)), .data)
    .data <- check_attr.tbl_ps(x=.data, recol=cols)
    res <- NextMethod()
    res <- add_attr.tbl_ps(x1=res, x2=.data, class="grouped_df_ps")
    return (res)
}

##' @method arrange phyloseq
##' @export
arrange.phyloseq <- function(.data, ..., by_group = FALSE){
    .data %<>% as_tibble()
    res <- arrange(.data = .data, ..., by_group = by_group)
    return (res)
}

##' @method arrange tbl_ps
##' @export
arrange.tbl_ps <- function(.data, ..., by_group = FALSE){
    res <- NextMethod()
    res <- add_attr.tbl_ps(x1 = res, x2 = .data)
    return(res)
}

##' @method arrange grouped_df_ps
##' @export
arrange.grouped_df_ps <- function(.data, ..., by_group = FALSE){
    res <- NextMethod()
    res <- add_attr.tbl_ps(x1 = res, x2 = .data, class = "grouped_df_ps")
    return (res)
}

##' @method left_join tbl_ps
##' @export
left_join.tbl_ps <- function(x, y, by=NULL, copy=FALSE, suffix = c(".x", ".y")){
    res <- NextMethod()
    res <- add_attr.tbl_ps(x1 = res, x2 = x) 
    res <- add_var(res, type="join")
    return(res)
}

add_attr.tbl_ps <- function(x1, x2, class="tbl_ps"){
    attr(x1, "samplevar") <- attr(x2, "samplevar")
    attr(x1, "mutatevar") <- attr(x2, "mutatevar")
    attr(x1, "taxavar") <- attr(x2, "taxavar")
    attr(x1, "fillNAtax") <- attr(x2, "fillNAtax")
    attr(x1, "tree") <- attr(x2, "tree")
    attr(x1, "refseq") <- attr(x2, "refseq")
    class(x1) <- add_class(new=class, old=class(x1))
    return(x1)   
}

add_class <- function(new, old){
    if (!new %in% old){
        return(c(new, old))
    }else{
        return (old)
    }
}

check_attr.tbl_ps <- function(x, recol, type="rename"){
    clnm <- colnames(x)
    renm <- clnm[recol]
    samplevar <- attr(x, "samplevar")
    taxavar <- attr(x, "taxavar")
    if (any(renm %in% c("Sample", "OTU", "Abundance")) && type=="rename"){
        stop("The Sample, OTU, and Abundance are not be renamed !")
    }
    if (any(renm %in% samplevar)){
        item1 <- intersect(renm, samplevar)
        indx <- match(item1, samplevar)
        indy <- match(item1, renm)
        samplevar[indx] <- names(recol[indy])
    }
    if (any(renm %in% taxavar)){
        item2 <- intersect(renm, taxavar)
        indx <- match(item2, taxavar)
        indy <- match(item2, renm)
        taxavar[indx] <- names(recol[indy])
    }
    attr(x, "samplevar") <- samplevar
    attr(x, "taxavar") <- taxavar
    return(x)
}

add_var <- function(x, type){
    cl <- colnames(x)
    samplevar <- attr(x, "samplevar")
    taxavar <- attr(x, "taxavar")
    mutatevar <- attr(x, "mutatevar")
    newvar <- setdiff(cl, c("OTU", "Sample", "Abundance", samplevar, taxavar, mutatevar))
    if (type == "mutate"){
        attr(x, "mutatevar") <- c(mutatevar, newvar)
    }else{
        attr(x, "samplevar") <- c(samplevar, newvar)
    }
    return(x)
}

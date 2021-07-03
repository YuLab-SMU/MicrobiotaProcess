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

##' @method filter MPSE
##' @importFrom ape keep.tip
##' @export
filter.MPSE <- function(.data, ..., .preserve = FALSE, .returnMPSE = FALSE){
    .data %<>% as_tibble()
    dots <- quos(...)
    res <- .data %>% filter(!!!dots, .preserve = .preserve)
    res <- add_attr.tbl_mpse(x1=res, x2=.data)
    if (!.returnMPSE){
        flag <- valid_names(res, type="tbl_mpse")
        xm <- tbl_mpse_return_message(flag)
        if (flag){
           xm <- c(xm, keep_mpse_message)
        }else{
           res %<>% tibble::as_tibble()
        }
        writeLines(xm)
    }else{
        if (valid_names(res)){
            res %<>% as.MPSE()
        }
    }    
    return (res)
}

#' @method filter tbl_mpse
#' @export
filter.tbl_mpse <- function(.data, ..., .preserve = FALSE){
    res <- NextMethod()
    if(valid_names(res, type="tbl_mpse")){
        res <- add_attr.tbl_mpse(x1=res, x2=.data)
    }
    return(res)
}

#' @method filter grouped_df_mpse
#' @export
filter.grouped_df_mpse <- function(.data, ..., .preserve=FALSE){
    res <- NextMethod()
    res <- add_attr.tbl_mpse(x1=res, x2=.data, class="grouped_df_mpse")
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

##' @method select MPSE 
##' @export
select.MPSE <- function(.data, ..., .returnMPSE = FALSE){
    .data %<>% as_tibble()
    loc <- tidyselect::eval_select(expr(c(...)), .data)
    #loc <- loc[!names(loc) %in% c("Sample", "OTU", "Abundance")]
    .data <- check_attr.tbl_mpse(x=.data, recol=loc, type="select")
    res <- select(.data=.data, ...)
    res <- add_attr.tbl_mpse(x1=res, x2=.data)
    if (!.returnMPSE){
        flag <- valid_names(res, type="tbl_mpse")
        xm <- tbl_mpse_return_message(flag)
        if (flag){
           xm <- c(xm, keep_mpse_message)
        }else{
           res %<>% tibble::as_tibble()
        }
        writeLines(xm)
    }else{
        if (valid_names(res)){
            res %<>% as.MPSE()
        }
    }
    return(res)
}

##' @method select tbl_mpse
##' @export
select.tbl_mpse <- function(.data, ...){
    loc <- tidyselect::eval_select(expr(c(...)), .data)
    #loc <- loc[!names(loc) %in% c("Sample", "OTU", "Abundance")]
    .data <- check_attr.tbl_mpse(x=.data, recol=loc, type="select")
    res <- NextMethod()
    if (valid_names(res, type="tbl_mpse")){
        res <- add_attr.tbl_mpse(x1=res, x2=.data)
    }else{
        res <- drop_class(res, class="tbl_mpse")
    }
    return (res)
}

##' @method select grouped_df_mpse
##' @export
select.grouped_df_mpse <- function(.data, ...){
    loc <- tidyselect::eval_select(expr(c(...)), .data)
    #loc <- loc[!names(loc) %in% c("Sample", "OTU", "Abundance")]
    .data <- check_attr.tbl_mpse(x=.data, recol=loc, type="select")
    res <- NextMethod()
    if (valid_names(res, type="grouped_df_mpse")){
        res <- add_attr.tbl_mpse(x1=res, x2=.data, class="grouped_df_mpse")
    }else{
        res <- drop_class(res, class=c("grouped_df_mpse", "tbl_mpse"))
    }
    return (res)
}

internal_select <- function(x, dots, ...){
    x <- as.data.frame(x)
    x %<>% select(!!!dots, ...)
    return(x)
}

##' @importFrom dplyr group_by_drop_default
##' @importFrom dplyr group_by
##' @method group_by MPSE
##' @export
group_by.MPSE <- function(.data, ..., .add = FALSE, .drop = group_by_drop_default(.data)){
    writeLines(tbl_mpse_return_message(TRUE))
    .data %<>% as_tibble()
    res <- group_by(.data=.data, ..., .add = .add, .drop=.drop)
    res <- add_attr.tbl_mpse(x1=res, x2=.data)
    class(res) <- c("grouped_df_mpse", class(res))
    return(res)
}

##' @method group_by tbl_mpse
##' @export
group_by.tbl_mpse <- function(.data, ..., .add = FALSE, .drop = group_by_drop_default(.data)){
    res <- NextMethod()
    if (valid_names(res, type="tbl_mpse")){
        res <- add_attr.tbl_mpse(x1=res, x2=.data)
        class(res) <- c("grouped_df_mpse", class(res))
    }else{
        res <- drop_class(res, class="tbl_mpse")
    }
    return(res)
}

##' @method ungroup grouped_df_mpse
##' @export
ungroup.grouped_df_mpse <- function(x, ...){
    res <- NextMethod()
    if (valid_names(res, type="tbl_mpse")){
        res <- add_attr.tbl_mpse(x1=res, x2=x)
    }
    res <- drop_class(res, class=c("grouped_df_mpse", "grouped_df"))
    return(res)
}

##' @method mutate MPSE
##' @export
mutate.MPSE <- function(.data, ...){
    writeLines(tbl_mpse_return_message(TRUE))
    .data %<>% as_tibble()
    res <- mutate(.data=.data, ...)
    res <- add_attr.tbl_mpse(x1=res, x2=.data)
    res <- add_var(res, type="mutate")
    return (res)
}

##' @method mutate tbl_mpse
##' @export
mutate.tbl_mpse <- function(.data, ...){
    res <- NextMethod()
    if (valid_names(res, type="tbl_mpse")){
        res <- add_attr.tbl_mpse(x1=res, x2=.data)
        res <- add_var(res, type="mutate")
    }else{
        res <- drop_class(res, class="tbl_mpse")
    }
    return(res)
}

##' @method mutate grouped_df_mpse
##' @export
mutate.grouped_df_mpse <- function(.data, ...){
    res <- NextMethod()
    if (valid_names(res, type="grouped_df_mpse")){
        res <- add_attr.tbl_mpse(x1=res, x2=.data, class="grouped_df_mpse")
        res <- add_var(res, type="mutate")
    }else{
        res <- drop_class(res, class=c("group_df_mpse", "tbl_mpse"))
    }
    return(res)
}

##' @method distinct MPSE
##' @export
distinct.MPSE <- function(.data, ..., .keep_all = FALSE){
    writeLines(tbl_mpse_return_message(TRUE))
    .data %<>%  as_tibble()
    res <- distinct(.data=.data, ..., .keep_all = .keep_all)
    res <- add_attr.tbl_mpse(x1 = res, x2 = .data)
    return(res)
}

##' @method distinct tbl_mpse
##' @export
distinct.tbl_mpse <- function(.data, ..., .keep_all = FALSE){
    res <- NextMethod()
    if (valid_names(res, type="tbl_mpse")){
        res <- add_attr.tbl_mpse(x1 = res, x2 = .data)
    }else{
        res <- drop_class(res, class="tbl_mpse")
    }
    return (res)
}

##' @method distinct grouped_df_mpse
##' @export
distinct.grouped_df_mpse <- function(.data, ..., .keep_all = FALSE){
    res <- NextMethod()
    if (valid_names(res, type="grouped_df_mpse")){
        res <- add_attr.tbl_mpse(x1 = res, x2 = .data, class="grouped_df_mpse")
    }else{
        res <- drop_class(res, class=c("grouped_df_mpse", "tbl_mpse"))
    }
    return (res)
}

##' @method rename MPSE
##' @export
rename.MPSE <- function(.data, ..., .returnMPSE=FALSE){
    .data %<>% as_tibble()
    cols <- tidyselect::eval_select(expr(c(...)), .data)
    .data <- check_attr.tbl_mpse(x=.data, recol=cols)
    res <- rename(.data=.data, ...)
    res <- add_attr.tbl_mpse(x1=res, x2=.data)
    if (!.returnMPSE){
        flag <- valid_names(res, type="tbl_mpse")
        xm <- tbl_mpse_return_message(flag)
        if (flag){
           xm <- c(xm, keep_mpse_message)
        }
        writeLines(xm)    
    }else{
        if (valid_names(res)){
            res %<>% as.MPSE()
        }
    }    
    return (res)
}

##' @method rename tbl_mpse
##' @importFrom tidyselect eval_select
##' @importFrom rlang expr
##' @export
rename.tbl_mpse <- function(.data, ...){
    cols <- tidyselect::eval_select(expr(c(...)), .data)
    .data <- check_attr.tbl_mpse(x=.data, recol=cols)
    res <- NextMethod()
    if (valid_names(res, type="tbl_mpse")){
        res <- add_attr.tbl_mpse(x1=res, x2=.data)
    }else{
        res <- drop_class(res, class="tbl_mpse")
    }
    return (res)
}

##' @method rename grouped_df_mpse
##' @export
rename.grouped_df_mpse <- function(.data, ...){
    cols <- tidyselect::eval_select(expr(c(...)), .data)
    .data <- check_attr.tbl_mpse(x=.data, recol=cols)
    res <- NextMethod()
    if (valid_names(res, type="grouped_df_mpse")){
        res <- add_attr.tbl_mpse(x1=res, x2=.data, class="grouped_df_mpse")
    }else{
        res <- drop_class(res, class=c("grouped_df_mpse", "tbl_mpse"))
    }
    return (res)
}

##' @method arrange MPSE
##' @export
arrange.MPSE <- function(.data, ..., by_group = FALSE, .returnMPSE=FALSE){
    .data %<>% as_tibble()
    res <- arrange(.data=.data, ..., by_group = FALSE)
    res <- add_attr.tbl_mpse(x1 = res, x2 = .data)
    if (!.returnMPSE){
        flag <- valid_names(res, type="tbl_mpse")
        xm <- tbl_mpse_return_message(flag)
        if (flag){
           xm <- c(xm, keep_mpse_message)
        }
        writeLines(xm)    
    }else{
        if (valid_names(res)){
            res %<>% as.MPSE()
        }
    }    
    return (res)
}

##' @method arrange tbl_mpse
##' @export
arrange.tbl_mpse <- function(.data, ..., by_group = FALSE){
    res <- NextMethod()
    if (valid_names(res, type="tbl_mpse")){
        res <- add_attr.tbl_mpse(x1 = res, x2 = .data)
    }else{
        res <- drop_class(res, class="tbl_mpse")
    }
    return(res)
}

##' @method arrange grouped_df_mpse
##' @export
arrange.grouped_df_mpse <- function(.data, ..., by_group = FALSE){
    res <- NextMethod()
    if (valid_names(res, type="grouped_df_mpse")){
        res <- add_attr.tbl_mpse(x1 = res, x2 = .data, class = "grouped_df_mpse")
    }else{
        res <- drop_class(res, class="tbl_mpse")
    }
    return (res)
}

##' @method left_join tbl_mpse
##' @export
left_join.tbl_mpse <- function(x, y, by=NULL, copy=FALSE, suffix = c(".x", ".y")){
    res <- NextMethod()
    if (valid_names(res, type="tbl_mpse")){
        res <- add_attr.tbl_mpse(x1 = res, x2 = x) 
        res <- add_var(res, type="join")
    }else{
        res <- drop_class(res, class="tbl_mpse")
    }
    return(res)
}

add_attr.tbl_mpse <- function(x1, x2, class="tbl_mpse"){
    attr(x1, "samplevar") <- attr(x2, "samplevar")
    attr(x1, "mutatevar") <- attr(x2, "mutatevar")
    attr(x1, "otumetavar") <- attr(x2, "otumetavar")
    attr(x1, "assaysvar") <- attr(x2, "assaysvar")
    attr(x1, "taxavar") <- attr(x2, "taxavar")
    attr(x1, "fillNAtax") <- attr(x2, "fillNAtax")
    otutree <- attr(x2, "otutree")
    taxatree <- attr(x2, "taxatree")
    refseq <- attr(x2, "refseq")
    rmotus <- setdiff(unique(x2$OTU), unique(x1$OTU))
    otutree <- .internal_drop.tip(tree=otutree, rmotus=rmotus)
    taxatree <- .internal_drop.tip(tree=taxatree, rmotus=rmotus, collapse.singles=FALSE)
    if (!is.null(refseq)){
        refseq <- refseq[!names(refseq) %in% rmotus]
    }
    attr(x1, "otutree") <- otutree
    attr(x1, "taxatree") <- taxatree
    attr(x1, "refseq") <- refseq
    class(x1) <- add_class(new=class, old=class(x1))
    return(x1)   
}

add_class <- function(new, old){
    x <- setdiff(new, old)
    if (length(x)>0){
        return(c(x, old))
    }else{
        return (old)
    }
}

check_attr.tbl_mpse <- function(x, recol, type="rename"){
    clnm <- colnames(x)
    renm <- clnm[recol]
    samplevar <- attr(x, "samplevar")
    taxavar <- attr(x, "taxavar")
    assaysvar <- attr(x, "assaysvar")
    otumetavar <- attr(x, "otumetavar")
    if (any(renm %in% c("Sample", "OTU", "Abundance")) && type=="rename"){
        stop("The Sample, OTU, and Abundance are not be renamed !")
    }
    item1 <- intersect(samplevar, renm)
    item2 <- intersect(taxavar, renm)
    item3 <- intersect(assaysvar, renm)
    item4 <- intersect(otumetavar, renm)
    if (any(renm %in% samplevar) && type=="rename"){
        indx <- match(item1, samplevar)
        indy <- match(item1, renm)
        samplevar[indx] <- names(recol[indy])
    }
    if (any(renm %in% taxavar) && type=="rename"){
        indx <- match(item2, taxavar)
        indy <- match(item2, renm)
        taxavar[indx] <- names(recol[indy])
    }
	if (any(renm %in% assaysvar) && type=="rename"){
        indx <- match(item3, assaysvar)
        indy <- match(item3, renm)
        assaysvar[indx] <- names(recol[indy])
	}
    if (any(renm %in% otumetavar) && type=="rename"){
        indx <- match(item4, otumetavar)
        indy <- match(item4, renm)
        otumetavar[indx] <- names(recol[indy])
    }
    if (type != "rename"){
        samplevar <- item1
        taxavar <- item2
        assaysvar <- item3
        otumetavar <- item4
    }
    attr(x, "samplevar") <- samplevar
    attr(x, "taxavar") <- taxavar
	attr(x, "assaysvar") <- assaysvar
    attr(x, "otumetavar") <- otumetavar
    return(x)
}

add_var <- function(x, type){
    cl <- colnames(x)
    samplevar <- attr(x, "samplevar")
    taxavar <- attr(x, "taxavar")
    mutatevar <- attr(x, "mutatevar")
    assaysvar <- attr(x, "assaysvar")
    otumetavar <- attr(x, "otumetavar")
    newvar <- setdiff(cl, c("OTU", "Sample", "Abundance", samplevar, taxavar, mutatevar, assaysvar, otumetavar))
    if (type == "mutate"){
        attr(x, "mutatevar") <- c(mutatevar, newvar)
    }else{
        attr(x, "samplevar") <- c(samplevar, newvar)
    }
    return(x)
}

valid_names <- function(x, type="MPSE"){
   flag <- all(c("OTU", "Sample", "Abundance") %in% colnames(x)) 
   if (!flag && type=="MPSE"){
       rlang::abort("The OTU, Sample and Abundance must be present to convert to MPSE object!")
   }else{
       return (flag)
   }
}

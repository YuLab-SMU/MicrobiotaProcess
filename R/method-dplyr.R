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
##' @export
filter.MPSE <- function(.data, ..., .preserve = FALSE, keep.mpse = TRUE){
    .data %<>% as_tibble()
    dots <- quos(...)
    res <- .data %>% filter(!!!dots, .preserve = .preserve)
    res <- add_attr.tbl_mpse(x1=res, x2=.data)
    if (nrow(res)==0){keep.mpse <- FALSE}
    if (!keep.mpse){
        flag <- valid_names(res, type="tbl_mpse")
        xm <- tbl_mpse_return_message(flag)
        if (flag){
           xm <- c(xm, keep_mpse_message())
        }else{
           res %<>% tibble::as_tibble()
        }
        message_wrap(xm)
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
select.MPSE <- function(.data, ..., keep.mpse = FALSE){
    .data %<>% as_tibble()
    loc <- tidyselect::eval_select(expr(c(...)), .data)
    #loc <- loc[!names(loc) %in% c("Sample", "OTU", "Abundance")]
    .data <- check_attr.tbl_mpse(x=.data, recol=loc, type="select")
    res <- select(.data=.data, ...)
    res <- add_attr.tbl_mpse(x1=res, x2=.data)
    if (! keep.mpse){
        flag <- valid_names(res, type="tbl_mpse")
        xm <- tbl_mpse_return_message(flag)
        if (flag){
           xm <- c(xm, keep_mpse_message())
        }else{
           res %<>% tibble::as_tibble()
        }
        message_wrap(xm)
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
    message_wrap(tbl_mpse_return_message(TRUE))
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
mutate.MPSE <- function(.data, keep.mpse = TRUE, ...){
    #writeLines(tbl_mpse_return_message(TRUE))
    .data %<>% as_tibble()
    res <- mutate(.data=.data, ...)
    res <- add_attr.tbl_mpse(x1=res, x2=.data)
    res <- add_var(res, type="mutate")
    if (! keep.mpse){
        flag <- valid_names(res, type="tbl_mpse")
        xm <- tbl_mpse_return_message(flag)
        if (flag){
           xm <- c(xm, keep_mpse_message())
        }else{
           res %<>% tibble::as_tibble()
        }
        message_wrap(xm)
    }else{
        if (valid_names(res)){
            res %<>% as.MPSE()
        }
    }
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
    message_wrap(tbl_mpse_return_message(TRUE))
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
rename.MPSE <- function(.data, ..., keep.mpse = TRUE){
    .data %<>% as_tibble()
    cols <- tidyselect::eval_select(expr(c(...)), .data)
    .data <- check_attr.tbl_mpse(x=.data, recol=cols)
    res <- rename(.data=.data, ...)
    res <- add_attr.tbl_mpse(x1=res, x2=.data)
    if (! keep.mpse){
        flag <- valid_names(res, type="tbl_mpse")
        xm <- tbl_mpse_return_message(flag)
        if (flag){
           xm <- c(xm, keep_mpse_message())
        }
        message_wrap(xm)    
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
arrange.MPSE <- function(.data, ..., by_group = FALSE, keep.mpse=FALSE){
    .data %<>% as_tibble()
    res <- arrange(.data=.data, ..., by_group = FALSE)
    res <- add_attr.tbl_mpse(x1 = res, x2 = .data)
    if (! keep.mpse){
        flag <- valid_names(res, type="tbl_mpse")
        xm <- tbl_mpse_return_message(flag)
        if (flag){
           xm <- c(xm, keep_mpse_message())
        }
        message_wrap(xm)    
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

##' @method left_join MPSE
##' @export
left_join.MPSE <- function(x, y, by=NULL, copy=FALSE, ...){
    if (inherits(y, 'dist') || all(unlist(lapply(y, function(i)inherits(i, 'dist'))))){
        y <- .dist2tbl(y, x)
    }
    dots <- rlang::quos(...)
    suffix <- c("", ".y")
    if ("suffix" %in% names(dots)){
        dots <- dots[names(dots)!="suffix"]
    }
    if (is.null(by)){
        if (all(c("OTU", "Sample") %in% colnames(y))){
            by <- c("OTU", "Sample")
        }else if ("OTU" %in% colnames(y)){
            by <- "OTU"
        }else if ("Sample" %in% colnames(y)){
            by <- "Sample"
        }else{
            rlang::abort(left_join_msg)
        }
    }

    idx <- names(by)
    if (length(by)==1){
        idnm <- ifelse(!is.null(idx), idx, by)
        if (!idnm %in% c("OTU", "Sample")){
            rlang::abort(left_join_msg)
        }
        idnm %<>% match.arg(c("OTU", "Sample"))
        if (idnm=="Sample"){
            sampleda <- x %>% mp_extract_sample()
            ornm <- sampleda %>% colnames() 
            sampleda <- sampleda %>% left_join(y, by=by, copy=copy, suffix=suffix, !!!dots)
            if (any(duplicated(sampleda$Sample))){
                sampleda %<>% .internal_nest(keepnm=ornm)
            }
            x@colData <- sampleda %>% 
                         tibble::column_to_rownames(var="Sample") %>%
                         S4Vectors::DataFrame()
        }
        if (idnm=="OTU"){
            featureda <- x %>% mp_extract_feature()
            ornm <- featureda %>% colnames()
            featureda %<>% left_join(y, by=by, copy=copy, suffix=suffix, !!!dots)
            if (any(duplicated(featureda$OTU))){
                featureda %<>% .internal_nest(keepnm=ornm)
            }
            SummarizedExperiment::rowData(x) <- featureda %>% 
                         tibble::column_to_rownames(var="OTU") %>%
                         S4Vectors::DataFrame() 
        }
    }else if (length(by==2)){
       if (is.null(idx)){
           idnm <- by
       }else{
           idx[which(nchar(idx)==0)] <- by[which(nchar(idx)==0)]
           idnm <- idx
       }
       names(by) <- idnm
       if (!all(idnm %in% c("OTU", "Sample"))){
           rlang::abort(left_join_msg)
       }else{
           y %<>% dplyr::rename(by)
           ynm <- y %>% select(-c("OTU", "Sample")) %>% lapply(is.numeric) %>% unlist()
           ynm <- ynm[ynm] %>% names()
           y <- x %>%
                mp_extract_assays(.abundance="Abundance") %>% 
                tibble::as_tibble(rownames="OTU") %>% 
                tidyr::pivot_longer(cols=!.data$OTU, names_to="Sample", values_to="Abundance") %>%
                left_join(y, suffix=suffix, by=by) %>%
                dplyr::mutate_if(is.numeric, ~replace(., is.na(.), 0)) 
           newassay <- .internal_build_assay(y, ynm)
           SummarizedExperiment::assays(x) <- c(SummarizedExperiment::assays(x), newassay)
       }
    }
    return(x)
}

##' @method left_join tbl_mpse
##' @export
left_join.tbl_mpse <- function(x, y, by=NULL, copy=FALSE, suffix = c("", ".y"), ...){
    if (inherits(y, 'dist') || all(unlist(lapply(y, function(i)inherits(i, 'dist'))))){
        y <- .dist2tbl(y, x)
    }
    res <- NextMethod()
    if (valid_names(res, type="tbl_mpse")){
        res <- add_attr.tbl_mpse(x1 = res, x2 = x) 
        res <- add_var(res, type="join")
    }else{
        res <- drop_class(res, class="tbl_mpse")
    }
    return(res)
}

##' @method pull MPSE
##' @export
pull.MPSE <- function(.data, var = -1, name = NULL, ...){
    var <- rlang::enquo(var)
    name <- rlang::enquo(name)
    da <- .data %>% 
          as_tibble() %>%
          dplyr::pull(var= !!var, name = !!name, ...)
    return(da)
}

##' @method slice MPSE
##' @export
slice.MPSE <- function(.data, ..., .preserve = FALSE){
    message_wrap(tbl_mpse_return_message(TRUE))
    dots <- rlang::quos(...)
    da <- .data %>% 
          as_tibble() %>%
          slice(!!!dots, .preserve=.preserve)
    return(da)
}

add_attr.tbl_mpse <- function(x1, x2, class="tbl_mpse"){
    #attr(x1, "samplevar") <- attr(x2, "samplevar")
    #attr(x1, "mutatevar") <- attr(x2, "mutatevar")
    #attr(x1, "otumetavar") <- attr(x2, "otumetavar")
    #attr(x1, "assaysvar") <- attr(x2, "assaysvar")
    #attr(x1, "taxavar") <- attr(x2, "taxavar")
    #attr(x1, "fillNAtax") <- attr(x2, "fillNAtax")
    #attr(x1, "internal_attr") <- attr(x2, "internal_attr")
    taxavar <- attr(x2, "taxavar")
    otutree <- attr(x2, "otutree")
    taxatree <- attr(x2, "taxatree")
    if (!is.null(taxatree)){
        nodeClass <- taxatree %>% 
                     tidytree::pull("nodeClass") %>% 
                     unique()
        if (!all(setdiff(nodeClass, taxavar) %in% c("OTU", "Root")) && length(taxavar) > 0){
            taxatree %<>% taxatree_to_tb() %>% select(taxavar) %>% convert_to_treedata2() 
        }else if (length(taxavar)==0){
            taxatree <- NULL
        }
    }
    refseq <- attr(x2, "refseq")
    otumetavar <- attr(x2, "otumetavar")
    if ("OTU" %in% colnames(x1)){
        rmotus <- setdiff(unique(x2$OTU), unique(x1$OTU))
        otutree <- .internal_drop.tip(tree=otutree, rmotus=rmotus)
        if (length(unique(x1$OTU))==1){
            if (!is.null(taxatree)){
                taxatb <- taxatree %>%
                          taxatree_to_tb() %>%
                          tibble::as_tibble(rownames="OTU") %>%
                          dplyr::filter(! .data$OTU %in% rmotus)
                taxatree <- NULL
                x1 <- merge(x1, taxatb, by="OTU")
                otumetavar <- c(colnames(taxatb)[colnames(taxatb)!="OTU"], otumetavar)
            }         
        }else{
            taxatree <- .internal_drop.tip(tree=taxatree, rmotus=rmotus, collapse.singles=FALSE)
        }
        if (!is.null(refseq)){
            refseq <- refseq[!names(refseq) %in% rmotus]
        }
        attr(x1, "otutree") <- otutree
        attr(x1, "taxatree") <- taxatree
        attr(x1, "refseq") <- refseq
    }
    attr(x1, "samplevar") <- attr(x2, "samplevar")
    attr(x1, "mutatevar") <- attr(x2, "mutatevar")
    attr(x1, "otumetavar") <- otumetavar
    attr(x1, "assaysvar") <- attr(x2, "assaysvar")
    attr(x1, "taxavar") <- taxavar
    attr(x1, "fillNAtax") <- attr(x2, "fillNAtax")
    attr(x1, "internal_attr") <- attr(x2, "internal_attr")
    x1 <- add_class(x=x1, new=class)
    return(x1)   
}

add_class <- function(x, new){
    xx <- setdiff(new, class(x))
    if (length(xx)>0){
        class(x) <- base::union(xx, class(x))
    }
    return (x)
    
}

check_attr.tbl_mpse <- function(x, recol, type="rename"){
    clnm <- colnames(x)
    renm <- clnm[recol]
    samplevar <- attr(x, "samplevar")
    taxavar <- attr(x, "taxavar")
    assaysvar <- attr(x, "assaysvar")
    otumetavar <- attr(x, "otumetavar")
    if (any(renm %in% c("Sample", "OTU", "Abundance")) && type=="rename"){
        rlang::abort("The Sample, OTU, and Abundance do not be renamed !")
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
        indname <- lapply(names(recol[indy]), function(i) i) %>% 
                   stats::setNames(taxavar[indx])
        taxavar[indx] <- names(recol[indy])
        taxa.tree <- attr(x, "taxatree")
        taxa.tree %<>% dplyr::mutate(nodeClass=do.call(dplyr::recode, c(.x=rlang::sym("nodeClass"), indname)))
        attr(x, "taxatree") <- taxa.tree
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
       return (flag && nrow(x)>0)
   }
}

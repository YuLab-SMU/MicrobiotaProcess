#' Calculate the correlation between the features with specified abundance.
#'
#' @rdname mp_cal_cor-methods
#' @param .data MPSE or tbl_mpse object
#' @param .abundance the name of otu abundance to be calculated
#' @param method character indicating which correlation coefficient to be calculate,
#' must be one of "pearson", "spearman", "kendall", default is "spearman".
#' @param alternative character indicates the alternative hypothesis,  must be one of
#' "two.sided", "greater" or "less", default is "two.sided".
#' @param R.square logical whether calculate the R square, default is FALSE.
#' @param p.adjust.method character correction method, default is NULL, which will 
#' not correct p value
#' @param action character, "add" joins the correlation data to the object, "only" return
#' a non-redundant tibble with the correlation information, "get" return 'tbl_graph' object,
#' which can be visualized by 'ggraph'.
#' @param ... additional parameters.
#' @return tbl_df (action='only') or tbl_graph (action='get') or MPSE object (action='add')
#' @export
setGeneric("mp_cal_cor", 
           function(
             .data, 
             .abundance, 
             method = 'spearman', 
             alternative = c("two.sided", "less", "greater"), 
             R.square = FALSE, 
             p.adjust.method = NULL, 
             action='add',
             ...)standardGeneric('mp_cal_cor'))

.internal_cal_cor <- function(.data, 
                            .abundance, 
                            method = 'spearman', 
                            alternative = c("two.sided", "less", "greater"), 
                            R.square=FALSE, 
                            p.adjust.method = NULL,
                            action = 'add', 
                            ...){

    .abundance <- rlang::enquo(.abundance)
    action <- match.arg(action, c('add', 'only', 'get'))
    
    da <- .data %>% mp_extract_assays(.abundance=!!.abundance)
    
    da <- .run_cor(da, 
             method = method, 
             alternative = alternative, 
             p.adjust.method = p.adjust.method,
             R.square = R.square
          )

    if (action == 'only'){
        return(da)
    }

    if (action == 'get'){
        da <- .convert_to_graph(da)
        fe.da <- .data %>% mp_extract_feature(addtaxa = TRUE)
        da <- da %>% dplyr::left_join(fe.da, by=c('label'= 'OTU'))
        return(da)
    }

    if (action == 'add'){
        da %<>% tidyr::nest(.by='x', .key = method) 
        .data %<>% left_join(da, by=c('OTU'='x'))
        return(.data)
    }

}

.run_cor <- function(data,
                     method = 'spearman',
                     alternative = c("two.sided", "less", "greater"),
                     p.adjust.method = NULL,
                     R.square = FALSE,
                     ...){
    n <- nrow(data)
    ia <- match.arg(alternative)
    p.dat <- r.dat <- matrix(data = NA, nrow=n, ncol=n)
    for (i in seq_len(n)){
        for (j in seq_len(n)){
            x <- .run_cor_item(data[i,], data[j,], method = method, alternative = ia, ...)
            p.dat[i, j] <- x[1]
            r.dat[i, j] <- x[2]
        }
    }

    #if (!is.null(p.adjust.method)){
    #    p.ad.dat <- matrix(p.adjust(p.dat, p.adjust.method), ncol=n, nrow=n)
    #}

    #if (R.square){
    #    r.sq.dat <- r.dat^2
    #}
    p.dat <- data.frame(p.dat)
    r.dat <- data.frame(r.dat)

    rownames(p.dat) <- colnames(p.dat) <- rownames(r.dat) <- colnames(r.dat) <- rownames(data)
    
    r.dat <- r.dat %>% .convert_td_cor(nm='R')
    r.dat$type <- ifelse(r.dat$R > 0, 'positive', 'negative')

    p.dat <- p.dat %>% .convert_td_cor(nm='p.val')

    r.p.dat <- r.dat %>% left_join(p.dat, by=c('x'='x', 'y'='y'))
    
    if (R.square){
        r.p.dat$R2 <- r.p.dat$R ^ 2
    }
    
    if (!is.null(p.adjust.method)){
        r.p.dat[[p.adjust.method]] <- p.adjust(r.p.dat$p.val, p.adjust.method)
    }

    return(r.p.dat)
}

.run_cor_item <- function(x, y, method = 'spearman', alternative='two.sided', ...){
    res <- stats::cor.test(as.numeric(x), as.numeric(y), method = method)
    res <- c(res$p.value, as.numeric(res$estimate))
    return(res)
}

.convert_td_cor <- function(x, nm){
    x <- x %>% 
        corrr::as_cordf() %>% 
        corrr::shave() %>% 
        corrr::stretch(na.rm=TRUE) 
    colnames(x)[match('r', colnames(x))] <- nm
    return(x)
}

.convert_to_graph <- function(d){
    tmp <- unique(c(d[[1]], d[[2]]))
    #tmp <- as.factor(c(d$x, d$y))
    #d$x <- as.numeric(tmp)[seq_len(nrow(d))]
    #d$y <- as.numeric(tmp)[seq_len(nrow(d)) + nrow(d)]
    
    #nodes <- data.frame(label=levels(tmp))
    nodes <- data.frame(label=tmp)
    #edges <- d

    tg <- tidygraph::tbl_graph(nodes = nodes, edges = d) #edges = edges)
    return(tg)
}


#' @rdname mp_cal_cor-methods
#' @aliases mp_cal_cor,MPSE
#' @exportMethod mp_cal_cor
setMethod("mp_cal_cor", signature(.data = 'MPSE'), .internal_cal_cor)

#' @rdname mp_cal_cor-methods
#' @aliases mp_cal_cor,tbl_mpse
#' @exportMethod mp_cal_cor
setMethod("mp_cal_cor", signature(.data="tbl_mpse"), .internal_cal_cor)

#' @rdname mp_cal_cor-methods
#' @aliases mp_cal_cor,grouped_df_mpse
#' @exportMethod mp_cal_cor
setMethod("mp_cal_cor", signature(.data="grouped_df_mpse"), .internal_cal_cor)

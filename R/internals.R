.dist_to_df <- function(x){
    nms <- attr(x, 'Labels')
    res <- mapply(.generate.paris, x = rep(list(nms), length(nms)), 
           y = nms, SIMPLIFY = FALSE) %>% do.call('rbind', .) 
    res$d <- as.numeric(x)
    res <- rbind(data.frame(x=nms, y=nms, d=0), res)
    res %<>% dplyr::arrange(!!rlang::sym("x"))
    return(res)
}

.generate.paris <- function(x, y){
    x <- x[-seq_len(which(x == y))]
    if (length(x)>0){
        d <- data.frame(x=y, y=x)
        return(d)
    }
}

.df_to_dist <- function(x, diag = FALSE, upper = FALSE){
    nms <- x$x[!duplicated(x$x)]
    x <- x[x$x != x$y,]
    d <- x$d
    attr(d, 'Labels') <- nms
    attr(d, "Size") <- length(nms)
    attr(d, 'call') <- call('as.dist')
    attr(d, 'diag') <- diag
    attr(d, 'upper') <- upper
    attr(d, 'method') <- NULL
    class(d) <- 'dist'
    return(d)
}



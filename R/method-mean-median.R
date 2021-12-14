compare_mean <- function(x, ...){
    UseMethod("compare_mean")
}

compare_mean.default <- function(x, y, ...){
    x.mean <- mean(x, na.rm=TRUE)
    y.mean <- mean(y, na.rm=TRUE)
    diffmean <- x.mean - y.mean
    return (list(x.mean=x.mean, y.mean=y.mean, diffmean=diffmean))
}

compare_mean.formula <- function(x, data, subset, na.action, ...){
    if(missing(x)
       || length(x) !=3L
       || length(attr(terms(x[-2L]),
                           "term.labels"))!=1L){
        stop("'formula' missing or incorrect, please check it.")
    }
    m <- match.call(expand.dots = FALSE)
        names(m)[[2L]] <- "formula"
    if (is.matrix(eval(m$data, parent.frame()))){
        m$data <- as.data.frame(data)
    }
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if (nlevels(g) !=2L){
        stop("grouping factor must have exactly 2 levels")
    }
    dat <- setNames(split(mf[[response]], g),
                   c("x", "y"))
    fc <- do.call("compare_mean", c(dat, list(...)))
    fc
}

compare_median <- function(x, ...){
    UseMethod("compare_median")
}

compare_median.default <- compare_mean.default
#compare_median.default <- function(x, y, ...){
#    x.new <- x[x!=0]
#    y.new <- y[y!=0]
#    if (length(x.new) ==0 ){
#        x.median <- 0
#    }else{
#        x.median <- median(x.new[x.new != 0], na.rm=TRUE)/sum(x, na.rm=TRUE)
#    }
#    if (length(y.new) == 0){
#        y.median <- 0
#    }else{
#        y.median <- median(y.new[y.new != 0], na.rm=TRUE)/sum(y, na.rm=TRUE)
#    }
#    diffmedian <- x.median - y.median
#    return (list(x.median = x.median, y.median = y.median, diffmedian = diffmedian))
#}

compare_median.formula <- function(x, data, subset, na.action, ...){
    if(missing(x)
       || length(x) !=3L
       || length(attr(terms(x[-2L]),
                           "term.labels"))!=1L){
        stop("'formula' missing or incorrect, please check it.")
    }
    m <- match.call(expand.dots = FALSE)
        names(m)[[2L]] <- "formula"
    if (is.matrix(eval(m$data, parent.frame()))){
        m$data <- as.data.frame(data)
    }
    m[[1L]] <- quote(stats::model.frame)
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if (nlevels(g) !=2L){
        stop("grouping factor must have exactly 2 levels")
    }
    dat <- setNames(split(mf[[response]], g),
                   c("x", "y"))
    fc <- do.call("compare_median", c(dat, list(...)))
    fc
}


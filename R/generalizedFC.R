#' @title generalized fold change
#' @description calculate the mean difference in a set of predefined quantiles of the logarithmic
#' @param x numeric vector, numeric vector of data values
#' @param y numeric vector, numeric vector of data values
#' @param base a positive or complex number, the base with respect to which logarithms are computed.  
#' default is 10.
#' @param steps positive numeric, increment of the sequence, default is 0.05.
#' @param pseudo positive numeric, avoid the zero for logarithmic, default is 0.00001.
#' @param formula, example 'Ozone ~ Month', Ozone is a numeric variable giving the data values
#' and ‘Month’ a factor giving the corresponding groups.
#' @param data data.frame, an optional matrix or data frame,containing the variables in the formula.
#' @param subset, (similar: see 'wilcox.test')an optional vector specifying a subset of observations to be used.
#' @param na.action, a function which indicates what should happen when the data, contain 'NA's.  
#' Defaults to 'getOption("na.action")'.
#' @param ..., additional arguments.
#' @author ShuangbinXu
#' @export
## reference https://www.nature.com/articles/s41591-019-0406-6
generalizedFC <- function(x, y, ...){
	UseMethod("generalizedFC")
}

#' @method generalizedFC default
#' @rdname generalizedFC
#' @export
generalizedFC.default <- function(x, y, base=10, steps=0.05, pseudo=0.00001, ...){
	x.q <- quantile(log(x+pseudo, base), probs=seq(.05, .95, steps))
	y.q <- quantile(log(y+pseudo, base), probs=seq(.05, .95, steps))
	x.mean <- mean(x)
	x.med <- median(x)
	y.mean <- mean(y)
	y.med <- median(y)
	fc <- sum(x.q - y.q)/length(x.q)
	res <- list(x.mean=x.mean, x.median=x.median, y.mean=y.mean, y.median=y.med, gfc=fc)
	attr(res, "class") <- "gFC"
	return(res)
}

#' @method generalizedFC formula
#' @rdname generalizedFC
#' @export
## refrence wilcox.test
generalizedFC.formula <- function(formula, data, subset, na.action, ...){
	if(missing(formula) 
	   || length(formula) !=3L 
	   || length(attr(terms(formula[-2L]),
				   "term.labels"))!=1L){
		stop("'formula' missing or incorrect, please check it.")
	}
	m <- match.call(expand.dots = FALSE)
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
	fc <- do.call("generalizedFC", c(dat, list(...)))
	fc
}


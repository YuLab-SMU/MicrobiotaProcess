#' @title Split Large Vector or DataFrame
#' 
#' @description
#' Split large vector or dataframe to list class, which contian subset vectors
#' or dataframe of origin vector or dataframe.
#'
#' @param x vector class or data.frame class.
#' @param nums integer.
#' @param chunks integer. use chunks if nums is missing.
#' Note nums and chunks shouldn't concurrently be NULL, 
#' default is NULL. 
#' @param random bool, whether split randomly, default is FALSE,
#' if you want to split data randomly, you can set TRUE, and
#' if you want the results are reproducible, you should add 
#' seed before.
#' @export
#' @author ShuangbinXu
splitData <- function(x, 
			 nums, 
			 chunks=NULL, 
			 random=FALSE){
       if (is.null(nums)){
		nums <- trunc(x/chunks)
	}
	if (missing(nums) && is.null(chunks)){
		stop("The nums and chunks shouldn't concurrently be missing!")
	}
	if (!is.null(nums) && !is.null(chunks)){
		message("We would use nums to split. and the chunk numbers is the length(x) provided by nums.")
	}
       if (is.vector(x) && length(x)>0){
		ind <- seedind(length(x), nums, random=random)	    
	}
	if (is.vector(x) && length(x)==0){
		stop("The legth of vector should be larger than zero, if the class of x is vector.")
	}
	if (is.data.frame(x)){
		ind <- seedind(nrow(x), nums, random=random)
	}
	return(split(x, ind))
}

seedind <- function(l, n, random=FALSE, seed=TRUE){
    tmpind <- rep(1:(trunc(l/n)+1),n)
    if (random){
	if ((l%%n)==0){
		tmpind <- tmpind[!tmpind %in% max(tmpind)]
	}
	tmpind <- sample(sort(tmpind)[1:l])
    }else{
    	tmpind <- sort(tmpind)[1:l]
    }
}


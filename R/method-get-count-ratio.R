#' @title calculate the count or relative abundance of replicate element with a speficify column
#' 
#' @description
#' Caculate the count or relative abundance of replicate element with a speficify columns
#' 
#' @param data dataframe; a dataframe contained one character column and others is numeric,
#' if featurelist is NULL. Or a numeirc dataframe, if featurelist is non't NULL, all columns 
#' should be numeric.
#' @param featurelist dataframe; a dataframe contained one chatacter column, default is NULL.
#' @param ... additional parameters.
#' @return mean of data.frame by featurelist
#' @export 
#' @author Shuangbin Xu
#' @importFrom plyr ddply 
#' @importFrom plyr numcolwise
#' @examples
#' otudafile <- system.file("extdata", "otu_tax_table.txt", 
#'                       package="MicrobiotaProcess")
#' samplefile <- system.file("extdata", 
#'                  "sample_info.txt", package="MicrobiotaProcess")
#' otuda <- read.table(otudafile, sep="\t", header=TRUE, 
#'                     row.names=1, check.names=FALSE, 
#'                     skip=1, comment.char="")
#' sampleda <- read.table(samplefile, 
#'             sep="\t", header=TRUE, row.names=1)
#' taxdf <- otuda[!sapply(otuda, is.numeric)]
#' taxdf <- split_str_to_list(taxdf)
#' otuda <- otuda[sapply(otuda, is.numeric)]
#' phycount <- get_count(otuda, taxdf[,2,drop=FALSE])
#' phyratios <- get_ratio(otuda, taxdf[,2,drop=FALSE])
get_count <- function(data,
                      featurelist, ...){ 
    if (missing(featurelist) || is.null(featurelist)){
        data <- data
        nums <- !unlist(lapply(data, is.numeric))
        group <- names(data[,nums,drop=FALSE])
    }else{
        data <- merge(data, featurelist, by=0)
        rownames(data) <- as.vector(data$Row.names)
        data$Row.names <- NULL
        group <- colnames(featurelist)
    }
    data <- data.frame(plyr::ddply(data, group, plyr::numcolwise(sum), ...), 
    	     check.names=FALSE, stringsAsFactors=FALSE)   
    rownames(data) <- as.vector(data[[group]])
    data[[group]] <- NULL
    #if (!isTRUE(countmode)){
    #   	data <- data.frame(prop.table(as.matrix(data), 2), check.names=FALSE, stringsAsFactors=FALSE)
    #	data <- data*multiplenum
    #}
    #if (isTRUE(rownamekeep)){
    #	data <- data.frame(cbind(feature=rownames(data), data), check.names=FALSE, stringsAsFactors=FALSE)
    #}
    return (data)
}


#' @rdname get_count
#' @export
get_ratio <- function(data, featurelist, ...){
    data <- get_count(data=data, featurelist=featurelist, ...)
    data <- data.frame(prop.table(as.matrix(data),2), check.names=FALSE)
    return(data)
}

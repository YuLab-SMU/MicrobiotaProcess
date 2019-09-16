#' @title caculate the count or relative abundance of replicate element with a speficify column
#' 
#' @description
#' Caculate the count or relative abundance of replicate element with a speficify columns
#' 
#' @param data dataframe; a dataframe contained one character column and others is numeric,
#' if featurelist is NULL. Or a numeirc dataframe, if featurelist is non't NULL, all columns 
#' should be numeric.
#' @param featurelist dataframe; a dataframe contained one chatacter column, default is NULL.
#' @param countmode boolean; whether return the count results (TRUE),or relative abundance (FALSE).
#' @param multiplenum numeric; the multiple you want to increase,default is 1.
#' @param rownamekeep boolean; whether you return a dataframe contained the rownames,default is FALSE.
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
#' taxdf <- splitStrtoList(taxdf)
#' otuda <- otuda[sapply(otuda, is.numeric)]
#' phycount <- CountOrRatios(otuda, taxdf[,2,drop=FALSE])
#' phyratios <- CountOrRatios(otuda, taxdf[,2,drop=FALSE], 
#'                            countmode=FALSE)
CountOrRatios <- function(data, 
                          featurelist, 
                          countmode=TRUE, 
                          multiplenum=1, 
                          rownamekeep=FALSE){
    if (missing(featurelist) || is.null(featurelist)){
        data <- data
    }else{
        data <- merge(data, featurelist, by=0)
        rownames(data) <- as.vector(data$Row.names)
        data$Row.names <- NULL
    }
    nums <- !unlist(lapply(data, is.numeric))
    group <- names(data[,nums,drop=FALSE])
    data <- data.frame(plyr::ddply(data, group, plyr::numcolwise(sum)), 
    	     check.names=FALSE, stringsAsFactors=FALSE)   
    rownames(data) <- as.vector(data[[group]])
    data[[group]] <- NULL
    if (!isTRUE(countmode)){
       	data <- data.frame(prop.table(as.matrix(data), 2), check.names=FALSE, stringsAsFactors=FALSE)
    	data <- data*multiplenum
    }
    if (isTRUE(rownamekeep)){
    	data <- data.frame(cbind(feature=rownames(data), data), check.names=FALSE, stringsAsFactors=FALSE)
    }
    return (data)
}


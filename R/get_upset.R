#' @title generate the dataset for upset of UpSetR
#' @param obj object, phyloseq or data.frame, if it is data.frame, 
#' the shape of it should be row sample * columns features.
#' @param sampleda data.frame, if the obj is data.frame, the sampleda
#' should be provided.
#' @param factorNames character, the column names of factor in sampleda
#' @param threshold integer, default is 0.
#' @param ..., additional parameters.
#' @return a data.frame for the input of `upset` of `UpSetR`.
#' @author Shuangbin Xu
#' @export
#' @examples
#' data(test_otu_data)
#' upsetda <- get_upset(test_otu_data, factorNames="group")
#' otudafile <- system.file("extdata", "otu_tax_table.txt",
#'                          package="MicrobiotaProcess")
#' samplefile <- system.file("extdata","sample_info.txt", 
#'                          package="MicrobiotaProcess")
#' otuda <- read.table(otudafile, sep="\t", header=TRUE, 
#'                     row.names=1, check.names=FALSE,
#'                     skip=1, comment.char="")
#' sampleda <- read.table(samplefile,sep="\t", 
#'                        header=TRUE, row.names=1)
#' head(sampleda)
#' otuda <- otuda[sapply(otuda, is.numeric)]
#' otuda <- data.frame(t(otuda), check.names=FALSE)
#' head(otuda[1:5, 1:5])
#' upsetda2 <- get_upset(obj=otuda, sampleda=sampleda, 
#'                      factorNames="group")
#' #Then you can use `upset` of `UpSetR` to visualize the results.
#' #library(UpSetR)
#' #upset(upsetda, sets=c("B","D","M","N"), sets.bar.color = "#56B4E9",
#' #      order.by = "freq", empty.intersections = "on")
setGeneric("get_upset", function(obj, ...)standardGeneric("get_upset"))

#' @aliases get_upset,data.frame
#' @rdname get_upset
#' @importFrom stats na.omit
#' @export
setMethod("get_upset", "data.frame", function(obj, sampleda, factorNames, threshold=0, ...){
    flaglen <- length(na.omit(match(rownames(obj),rownames(sampleda))))
    sampleda <- sampleda[,match(factorNames, colnames(sampleda)),drop=FALSE]
    if (flaglen==0){
        stop("The sample names of obj and sampleda should be consistent!
              Please check the rownames of obj and rownames of sampleda!") 
    }
    if (flaglen!=0 & flaglen < nrow(obj)){
        message("There are some sample names are not consistent!")
    }
    dameta <- merge(obj, sampleda, by=0)
    rownames(dameta) <- as.vector(dameta$Row.names)
    dameta$Row.names <- NULL
    dameta <- count_or_ratios(dameta)
    daupset <- apply(dameta, 1, 
                     function(x){unlist(lapply(x, function(x){if(x>threshold){1}else{0}}))})
    daupset <- data.frame(daupset, check.names=FALSE)
    return(daupset)
})

#' @aliases get_upset,phyloseq
#' @rdname get_upset
#' @export
setMethod("get_upset", "phyloseq", function(obj,...){
    otuda <- checkotu(obj)
    sampledata <- checksample(obj)
    daupset <- get_upset(obj=otuda, sampleda=sampledata,...)
    return(daupset)
})


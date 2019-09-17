#' @title alpha index
#' 
#' @description
#' caculate the alpha index (Obseve,Chao1,Shannon,Simpson) of sample
#' with \code{\link[vegan]{diversity}}
#' @param obj object, data.frame of (nrow sample * ncol taxonomy(feature)) 
#' or phyloseq.
#' @param mindepth numeric, Subsample size for rarefying community.
#' @param ..., additional arguments.
#' @return data.frame contained alpha Index.
#' @author ShuangbinXu
#' @export
#' @examples
#' library(tidyverse)
#' library(vegan)
#' otudafile <- system.file("extdata", "otu_tax_table.txt", 
#'                         package="MicrobiotaProcess")
#' otuda <- read.table(otudafile, sep="\t", 
#'              header=TRUE, row.names=1, 
#'              check.names=FALSE, skip=1, comment.char="")
#' otuda <- otuda[sapply(otuda, is.numeric)] %>% t() %>% 
#'           data.frame(check.names=FALSE)
#' set.seed(1024)
#' alphatab <- alphaindex(otuda)
#' head(alphatab)
#' data(test_otu_data)
#' class(test_otu_data)
#' set.seed(1024)
#' alphatab2 <- alphaindex(test_otu_data)
#' head(alphatab2)
setGeneric("alphaindex", function(obj, ...){standardGeneric("alphaindex")})

#' @aliases alphaindex,matrix
#' @rdname alphaindex
#' @importFrom vegan rrarefy estimateR diversity specnumber
#' @export
setMethod("alphaindex", "matrix", function(obj, mindepth,...){
    if (!identical(all.equal(obj, round(obj)),TRUE)){
           stop("the data should be integer (counts)!")
    }
    if (missing(mindepth) || is.null(mindepth)){
           mindepth = min(rowSums(obj))
    }
    obj <- rrarefy(obj, mindepth)
    Chao <- estimateR(obj)
    Shannon <- diversity(obj)
    Simpson <- diversity(obj, index="simpson")
    J <- Shannon/log(specnumber(obj))
    alpha <- data.frame(Observe=Chao[1,],
                        Chao1=Chao[2,],
                        ACE=Chao[4,],
                        Shannon,
                        Simpson,
                        J)
    return(alpha)
})     

#' @aliases alphaindex,data.frame
#' @rdname alphaindex
#' @export
setMethod("alphaindex", "data.frame", function(obj, ...){
    obj <- obj[,colSums(obj)>0,drop=FALSE]
    obj <- as.matrix(obj)
    alpha <- alphaindex(obj, ...)
    return(alpha)
})

#' @aliases alphaindex,integer
#' @rdname alphaindex
#' @export
setMethod("alphaindex", "integer", function(obj, ...){
    obj <- obj[obj>0]
    obj <- as.matrix(obj)
    alpha <- alphaindex(obj, ...)
    return(alpha)
})

#' @aliases alphaindex,phyloseq
#' @rdname alphaindex
#' @export
setMethod("alphaindex", "phyloseq", function(obj, ...){
    otuda <- checkotu(obj)
    alpha <- alphaindex(otuda, ...)
    return(alpha)
})


#' @title alpha index
#' 
#' @description
#' calculate the alpha index (Obseve,Chao1,Shannon,Simpson) of sample
#' with \code{\link[vegan]{diversity}}
#' @param obj object, data.frame of (nrow sample * ncol taxonomy(feature)) 
#' or phyloseq.
#' @param mindepth numeric, Subsample size for rarefying community.
#' @param sampleda data.frame,sample information, row sample * column factors.
#' @param ... additional arguments.
#' @return data.frame contained alpha Index.
#' @author Shuangbin Xu
#' @rdname get_alphaindex
#' @export
#' @examples
#' otudafile <- system.file("extdata", "otu_tax_table.txt", 
#'                         package="MicrobiotaProcess")
#' otuda <- read.table(otudafile, sep="\t", 
#'              header=TRUE, row.names=1, 
#'              check.names=FALSE, skip=1, comment.char="")
#' otuda <- otuda[sapply(otuda, is.numeric)] %>% t() %>% 
#'           data.frame(check.names=FALSE)
#' set.seed(1024)
#' alphatab <- get_alphaindex(otuda)
#' head(as.data.frame(alphatab))
#' data(test_otu_data)
#' class(test_otu_data)
#' set.seed(1024)
#' alphatab2 <- get_alphaindex(test_otu_data)
#' head(as.data.frame(alphatab2))
setGeneric("get_alphaindex", function(obj, ...){standardGeneric("get_alphaindex")})

#' @aliases get_alphaindex,matrix
#' @rdname get_alphaindex
#' @importFrom vegan estimateR diversity specnumber
#' @export
setMethod("get_alphaindex", "matrix", function(obj, mindepth, sampleda,...){
    if (!identical(all.equal(obj, round(obj)),TRUE)){
           stop("the data should be integer (counts)!")
    }
    if (missing(mindepth) || is.null(mindepth)){
           mindepth <- min(rowSums(obj))
    }
    obj <- vegan::rrarefy(obj, mindepth)
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
    if (missing(sampleda)){
        sampleda <- NULL
    }
    alpha <- new("alphasample",
                 alpha=alpha,
		 sampleda=sampleda)
    return(alpha)
})     

#' @aliases get_alphaindex,data.frame
#' @rdname get_alphaindex
#' @export
setMethod("get_alphaindex", "data.frame", function(obj, ...){
    obj <- obj[,colSums(obj)>0,drop=FALSE]
    obj <- as.matrix(obj)
    alpha <- get_alphaindex(obj, ...)
    return(alpha)
})

#' @aliases get_alphaindex,integer
#' @rdname get_alphaindex
#' @export
setMethod("get_alphaindex", "integer", function(obj, ...){
    obj <- obj[obj>0]
    obj <- as.matrix(obj)
    alpha <- get_alphaindex(obj, ...)
    alpha <- alpha@alpha
    return(alpha)
})

#' @aliases get_alphaindex,numeric
#' @rdname get_alphaindex
#' @export
setMethod("get_alphaindex", "numeric",function(obj, ...){
    obj <- obj[obj>0]
    obj <- as.matrix(obj)
    alpha <- get_alphaindex(obj, ...)
    alpha <- alpha@alpha
    return(alpha)
})

#' @aliases get_alphaindex,phyloseq
#' @rdname get_alphaindex
#' @export
setMethod("get_alphaindex", "phyloseq", function(obj, ...){
    otuda <- checkotu(obj)
    sampleda <- checksample(obj)
    alpha <- get_alphaindex(obj=otuda, sampleda=sampleda,...)
    return(alpha)
})

#' @aliases get_alphaindex,tbl_ps
#' @rdname get_alphaindex
#' @export
setMethod("get_alphaindex", "tbl_ps", function(obj, ...){
    obj <- obj %>% as.phyloseq()
    alpha <- get_alphaindex(obj, ...)
    return(alpha)
})

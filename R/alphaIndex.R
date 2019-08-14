#' @title alpha index
#' 
#' @description
#' caculate the alpha index (Obseve,Chao1,Shannon,Simpson) of sample
#' with \code{\link[vegan]{diversity}}
#' @param data data.frame, (nrow sample * ncol taxonomy(feature))
#' @param mindepth numeric, Subsample size for rarefying community.
#' @return data.frame contained alpha Index.
#' @author ShuangbinXu
#' @importFrom vegan rrarefy estimateR diversity specnumber
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
#' otuda <- rrarefy(otuda, min(rowSums(otuda)))
#' alphatab <- alphaindex(otuda)
#' head(alphatab)
alphaindex <- function(data,mindepth=NULL){
    if (is.data.frame(data)){
           data <- data[,colSums(data)>0,drop=FALSE]
    }else{
           data <- data[data>0]}
    x <- as.matrix(data)
    if (!identical(all.equal(x, round(x)),TRUE)){
           stop("the data should be integer (counts)")
    }
    if (is.null(mindepth) || missing(mindepth)){
           mindepth = min(rowSums(data))
    }
    x <- rrarefy(data, mindepth)
    Chao <- estimateR(x)
    Shannon <- diversity(x)
    Simpson <- diversity(x, index="simpson")
    J <- Shannon/log(specnumber(data))
    alpha <- data.frame(Observe=Chao[1,],
                        Chao1=Chao[2,],
                        ACE=Chao[4,],
                        Shannon,
                        Simpson,
                        J)
    #attr(alpha, "indexs") <- alpha
    #attr(alpha, "class") <- "Alpha"
    return(alpha)
}     


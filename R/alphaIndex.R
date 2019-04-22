#' @title alpha index
#' 
#' @description
#' caculate the alpha index (Obseve,Chao1,Shannon,Simpson) of sample
#' with \code{\link[vegan]}
#' @param data data.frame, (nrow sample * ncol taxonomy(feature))
#' @param mindepth integer, Subsample size for rarefying community.
#' @author ShuangbinXu
#' @importFrom vegan rrarefy estimateR diversity specnumber
#' @export
alphaindex <- function(data,
                                      mindepth=NULL){
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
       return(alpha)
}     


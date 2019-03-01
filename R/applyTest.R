#' @title 
#' 
#' @description
#' @param
#' @param
#' @param
#' @param
#'
#' @export
#' @author Shuangbin Xu
#' @importFrom agricolae kruskal
MultiCompare <- function(data,
                         measure,
                         variable,
                         padj="fdr", ...){
    nums <- !unlist(lapply(data, is.numeric))
    if (length(num)==1){
	 variable <- names(data[,nums,drop=FALSE])
    }
    compare <- with(data, agricolae::kruskal(#eval(parse(text=measure)), eval(parse(text=variable)),
                                             data[[measure]],
                                             data[[variable]],
                                             p.adj=padj,
                                             group=FALSE, ...))
    comparisonDf <- compare$comparison
    comparisonDf <- cbind(measures=measure,
                          comparisonDf,
                          stringsAsFactors=FALSE)
    comparisonDf$compareGroup <- rownames(comparisonDf)
    return(comparisonDf)
}

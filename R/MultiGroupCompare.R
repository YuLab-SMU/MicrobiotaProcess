#' @title multiple comparison with Kruskal Wallis test
#' 
#' @description
#' Multiple comparision with Kruskal Wallis test. which based on the
#' \code{\link[agricolae]{kruskal}}
#' @param data data frame; n(row) sample * m(col) feature and a factor col.
#' @param measure character; must be the names of feature.
#' @param padj character; Method for adjusting p values, see p.adjust, default is `fdr`.
#' @param ...; Additional arguments passed to \code{\link[agricolae]{kruskal}}.
#'
#' @author Shuangbin Xu
#' @export
#' @importFrom agricolae kruskal
MultiGroupCompare <- function(data,
                         measure,
                         padj="fdr", ...){
	colnames(data) <- make.names(colnames(data))
    numstmp <- !unlist(lapply(data, is.numeric))
	variable <- names(data[,numstmp, drop=FALSE])
	measure <- make.names(measure)
    compare <- with(data, agricolae::kruskal(eval(parse(text=measure)),
											 eval(parse(text=variable)),
                                             p.adj=padj,
                                             group=FALSE, ...))
    comparisonDf <- compare$comparison
    comparisonDf <- cbind(measures=measure,
                          comparisonDf,
                          stringsAsFactors=FALSE)
    comparisonDf$compareGroup <- rownames(comparisonDf)
    return(comparisonDf)
}

#' @title multiple comparison for mulitple measures with Kruskal Wallis test   
#' 
#' @description
#' Multiple comparision for multiple measures with Kruskal Wallis test. which 
#' based on the \code{\link[MicrobitaProcess]{MultiGroupCompare}}
#' @param data data frame; n(row) sample * m(col) feature and a factor col.
#' @param measurelist vector character; must be the names of feature, 
#' default is NULL,which perform the test for the total features.
#' @param padjust character; Method for adjusting p values, see 
#' p.adjust, default is `fdr`.
#' @param ...; Additional arguments passed to 
#' \code{\link[MicrobitaProcess]{MultiGroupCompare}}.
#'
#' @export
#' @author Shuangbin Xu


mapplyMultiGroupCompare <- function(data, measurelist=NULL, padjust="fdr", ...){
	colnames(data) <- make.names(colnames(data))       
	nums <- unlist(lapply(data, is.numeric))   
	if (is.null(measurelist)){
		measurelist <- names(data[,nums,drop=FALSE])
	}
	multiMeasureTest <- mapply(MultiGroupCompare, 
		   measurelist,
		   MoreArgs=list(data=data,
						 padj=padjust,
						 ...),
		   SIMPLIFY=FALSE
		   )
	return(multiMeasureTest)
}

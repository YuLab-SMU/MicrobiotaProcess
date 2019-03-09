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
#' @export
#' @author Shuangbin Xu
#' @importFrom agricolae kruskal
MultiGroupCompare <- function(data,
                         measure,
                         padj="fdr", ...){
	colnames(data) <- make.names(colnames(data))
    nums <- !unlist(lapply(data, is.numeric))
	variable <- names(data[,nums,drop=FALSE])
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
#' @param measure vector character; must be the names of feature.
#' @param MultiGroupCompare character; Method for adjusting p values, see 
#' p.adjust, default is `fdr`.
#' @param ...; Additional arguments passed to 
#' \code{\link[MicrobitaProcess]{MultiGroupCompare}}.
#'
#' @export
#' @author Shuangbin Xu
#' @importFrom MicrobitaProcess MultiGroupCompare


mapplyMultiGroupCompare <- function(data, padjust="fdr",...){
	colnames(data) <- make.names(colnames(data))       
	nums <- unlist(lapply(data, is.numeric))   
	measurelist <- names(data[,nums,drop=FALSE])
	multiMeasureTest <- mapply(MultiGroupCompare, 
		   measurelist,
		   MoreArgs=list(data=data,
						 p.adj=padjust,
						 ...),
		   SIMPLIFY=FALSE
		   )
	return(multiMeasureTest)
}

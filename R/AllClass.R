#' @title prcomp class
#' @seealso \code{\link[stats]{prcomp}}
#' @name prcomp-class
#' @rdname prcomp-class
#' @keywords internal
setOldClass("prcomp")

#' @title pcoa class
#' @seealso \code{\link[ape]{pcoa}}
#' @name pcoa-class
#' @rdname pcoa-class
#' @keywords internal
setOldClass("pcoa")

#' @keywords internal
setClassUnion("prcompOrNull", c("prcomp", "pcoa", "NULL"))

#' @keywords internal
setClassUnion("dataframeOrNull", c("data.frame", "NULL"))

#' @title pcasample class
#' @name pcasample-class
#' @rdname pcasample-class
#' @exportClass pcasample
setClass("pcasample",
		 representation=representation(
					pca="prcompOrNull",
					#varcontrib="VarContrib",
					sampleda="dataframeOrNull"),
		 prototype=prototype(pca=NULL, sampleda=NULL))

#' @keywords internal
setClassUnion("matrixOrNull", c("matrix", "NULL"))

#' @keywords internal
setClassUnion("characterOrNull", c("character", "NULL"))

#' @title ordplotClass class
#' @name ordplotClass-class
#' @rdname ordplotClass-class
#' @exportClass ordplotClass
setClass("ordplotClass",
		 representation=representation(coord="matrixOrNull",
									   xlab="characterOrNull",
									   ylab="characterOrNull",
									   title="characterOrNull"),
		 prototype=prototype(coord=NULL, 
							 xlab=NULL, 
							 ylab=NULL, 
							 title=NULL))

#' @importClassesFrom phyloseq phylo
#' @keywords internal 
setClassUnion("phyloOrNULL", c("phylo", "NULL"))


#' @title clustplotClass class
#' @name clustplotClass-clcass
#' @rdname clustplotClass-class
#' @exportClass clustplotClass
setClass("clustplotClass",
		 representation=representation(hclustphylo="phyloOrNULL",
									   sampleda="dataframeOrNull",
									   distmethod="characterOrNull"
									   ),
		 prototype=prototype(hclustphylo=NULL,
							 sampleda=NULL,
							 distmethod=NULL))

#' @keywords internal
setClassUnion("listOrNull", c("list", "NULL"))

#' @title diffAnalysisClass class
#' @name diffAnalysisClass-class
#' @rdname diffAnalysisClass-class
#' @exportClass diffAnalysisClass
setClass("diffAnalysisClass",
		 representation=representation(originalD="dataframeOrNull",
									   sampleda="dataframeOrNull",
									   taxda="dataframeOrNull",
									   kwres="dataframeOrNull",
									   secondvars="listOrNull",
									   mlres="dataframeOrNull"
									   ),
		 prototype=prototype(originalD=NULL,
							 sampleda=NULL,
							 taxda=NULL,
							 kwres=NULL,
							 secondvars=NULL,
							 mlres=NULL))

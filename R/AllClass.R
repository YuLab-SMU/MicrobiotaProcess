# #' @name phyloseq-class
# #' @rdname phyloseq-class
# #' @importClassesFrom phyloseq phyloseq
#setOldClass("phyloseq")

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

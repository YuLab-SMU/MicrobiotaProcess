# #' @name phyloseq-class
# #' @rdname phyloseq-class
# #' @importClassesFrom phyloseq phyloseq
#setOldClass("phyloseq")

#' @title prcomp class
#' @seealso \code{\link[stats]{prcomp}}
#' @name prcomp-class
#' @rdname prcomp-class
#' @keywords internal
#prcomp <- structure(list(), class = "prcomp")
setOldClass("prcomp")

#' @title pcasample class
#' @name pcasample-class
#' @rdname pcasample-class
#' @exportClass pcasample
setClass("pcasample",
		 representation=representation(
					pca="prcomp",
					#varcontrib="VarContrib",
					sampleda="data.frame"),
		 prototype=prototype(pca=NULL, sampleda=NULL))

#' @title ordplotClass class
#' @name ordplotClass-class
#' @rdname ordplotClass-class
#' @exportClass ordplotClass
setClass("ordplotClass",
		 representation=representation(coord="matrix",
									   xlab="character",
									   ylab="character",
									   title="character"),
		 prototype=prototype(coord=NULL, 
							 xlab=NULL, 
							 ylab=NULL, 
							 title=NULL))

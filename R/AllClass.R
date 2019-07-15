#' @name phyloseq-class
#' @rdname phyloseq-class
#' @importClassesFrom phyloseq phyloseq
setOldClass("phyloseq")

#' @importClassesFrom stats prcomp
setOldClass("prcomp")

#' @name pcasample-class
#' @rdname pcasample-class
#' @exportClass pcasample
setClass("pcasample",
		 representation=representation(
					pca="prcomp",
					#varcontrib="VarContrib",
					sampleda="data.frame"),
		 prototype=prototype(pca=NULL, sampleda=NULL))



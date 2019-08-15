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
#' @docType class
#' @slot pca prcomp or pcoa object 
#' @slot sampleda associated sample information
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
#' @docType class
#' @slot coord matrix object contained the coordinate 
#' for ordination plot.
#' @slot xlab character object contained the text of xlab
#' for ordination plot.
#' @slot ylab character object contained the text of ylab
#' for ordination plot.
#' @slot title character object contained the text of title
#' for ordination plot.
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
#' @docType class
#' @slot hclustphylo phylo object (convert hclust to phylo).
#' @slot sampleda assocaited sample information.
#' @slot distmethod character the method of dist.
#' @name clustplotClass-class
#' @rdname clustplotClass-class
#' @exportClass clustplotClass
setClass("clustplotClass",
    representation=representation(hclustphylo="phyloOrNULL",
    sampleda="dataframeOrNull",
    distmethod="characterOrNull"),
    prototype=prototype(hclustphylo=NULL,
    sampleda=NULL,
    distmethod=NULL))

#' @keywords internal
setClassUnion("listOrNull", c("list", "NULL"))

#' @keywords internal
setClassUnion("numericOrNull", c("numeric", "NULL"))

#' @title diffAnalysisClass class
#' @docType class
#' @slot originalD original feature data.frame.
#' @slot sampleda associated sample information.
#' @slot taxda the data.frame contained taxonomy.
#' @slot kwres the results of first test, contained
#' feature names, pvalue and fdr.
#' @slot secondvars the results of second test, contained
#' features names, gfc (TRUE representation the relevant 
#' feantures is enriched in relevant factorNames), 
#' Freq(the number of TRUE or FALSE), factorNames.
#' @slot mlres the results of LDA or randomForest,
#' @slot classname character the factor names.
#' @slot normalization the argument of 
#' \code{\link[MicrobiotaProcess]{diffAnalysis}}
#' @name diffAnalysisClass-class
#' @rdname diffAnalysisClass-class
#' @exportClass diffAnalysisClass
setClass("diffAnalysisClass",
    representation=representation(originalD="dataframeOrNull",
    sampleda="dataframeOrNull",
    taxda="dataframeOrNull",
    kwres="dataframeOrNull",
    secondvars="listOrNull",
    mlres="dataframeOrNull",
    classname="characterOrNull",
    normalization="numericOrNull"),
    prototype=prototype(originalD=NULL,
    sampleda=NULL,
    taxda=NULL,
    kwres=NULL,
    secondvars=NULL,
    mlres=NULL,
    classname=NULL,
    normalization=NULL))

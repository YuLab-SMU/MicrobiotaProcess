#' @title prcomp class
#' @seealso \code{\link[stats]{prcomp}}
#' @name prcomp-class
#' @rdname prcomp-class
#' @keywords internal
setOldClass("prcomp")

#' @title tbl_mpse class
#' @name tbl_mpse-class
#' @rdname tbl_mpse-class
#' @keywords internal
#' @noRd
setOldClass("tbl_mpse")

#' @noRd
TreeSummarizedExperiment <- methods::getClassDef("TreeSummarizedExperiment", "TreeSummarizedExperiment") %>% 
                              suppressMessages()

#' @noRd
phyloseq <- methods::getClassDef("phyloseq", "phyloseq") %>% suppressMessages()

#' @title grouped_df_mpse class
#' @name grouped_df_mpse-class
#' @rdname grouped_df_mpse-class
#' @noRd
setOldClass("grouped_df_mpse")

#' @importClassesFrom tidytree treedata
setClassUnion("TREEDATA_OR_NULL", c("treedata", "NULL"))

#' @importClassesFrom Biostrings XStringSet
setClassUnion("XSTRINGSET_OR_NULL", c("XStringSet", "NULL"))

#' @title MPSE class
#' @docType class
#' @slot otutree A treedata object of tidytree package or NULL.
#' @slot taxatree A treedata object of tidytree package or NULL.
#' @slot refseq A XStringSet object of Biostrings package or NULL.
#' @slot ... Other slots from \code{\link[SummarizedExperiment:SummarizedExperiment]{SummarizedExperiment}}
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @exportClass MPSE
setClass("MPSE",
    contains = "SummarizedExperiment",
    slots    = c(
        otutree  = "TREEDATA_OR_NULL",
        taxatree = "TREEDATA_OR_NULL",
        refseq   = "XSTRINGSET_OR_NULL"
    ),
    prototype = list(
        otutree  = NULL,
        taxatree = NULL,
        refseq   = NULL    
    )
)

#' @title Construct a MPSE object
#' @param assays A 'list' or 'SimpleList' of matrix-like elements
#' All elements of the list must have the same dimensions, and must have
#' the names, eg. list(Abundance=xx1, RareAbundance=xx2).
#' @param colData An optional DataFrame describing the samples.
#' @param otutree A treedata object of tidytree package
#' @param taxatree A treedata object of tidytree package
#' @param refseq A XStingSet object of Biostrings package
#' @param ... additional parameters, see also the usage 
#' of \code{\link[SummarizedExperiment]{SummarizedExperiment}}.
#' @return MPSE object
#' @export
MPSE <- function(assays,
                 colData,
                 otutree = NULL, 
                 taxatree = NULL, 
                 refseq = NULL, 
                 ...){
    se <- SummarizedExperiment::SummarizedExperiment(assays=assays, colData=colData, ...)
    mpse <- new("MPSE",
                se,
                otutree = otutree,
                taxatree = taxatree,
                refseq = refseq
               )
}

.valid.MPSE <- function(object){
    if (!is.null(object@otutree)){
        ntip <- treeio::Ntip(object@otutree)
        if (nrow(object)!=ntip){
            rlang::abort(c("The number of tip labels of otutree is not equal the number of otu in assays.", 
                         "Please check the otutree or assays!"))
        }
        if (!all(object@otutree@phylo$tip.label %in% rownames(object))){
            rlang::abort(c("Some otu names of otutree are different with the otu names in assays.", 
                         "Please check the otutree or assays"))
        }
    }
    if (!is.null(object@taxatree)){
        ntip <- treeio::Ntip(object@taxatree)
        if (nrow(object) != ntip){
            rlang::abort(c("The number of tip labels of taxatree is not equal the number of otu in assays.", 
                           "Please check the taxatree or assays!"))
        }
        if (length(intersect(object@taxatree@phylo$tip.label, rownames(object)))!=ntip){
            rlang::abort(c("Some otu names of taxatree are different with the otu names in assays.", 
                           "Please check the taxatree or assays."))
        }
    }
    if (!is.null(object@refseq)){
       if (length(object@refseq) != nrow(object) || length(object@refseq)==0){
           rlang::abort(c("The number of reference sequences is not equal the numbers of the otu in assays.",
                        "Please check the refseq or assays"))
       }
       if (!all(names(object@refseq) %in% rownames(object))){
           rlang::abort(c("Some reference sequence names are different with the otu names in assays.", 
                        "Please check the refseq or assays"))
       }
    }
    return(NULL)
}

setValidity("MPSE", .valid.MPSE)

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

#' @title alphasample class
#' @docType class
#' @slot alpha data.frame contained alpha metrics of samples
#' @slot sampleda associated sample information
#' @name alphasample-class
#' @rdname alphasample-class
#' @exportClass alphasample
setClass("alphasample",
    representation=representation(
        alpha="dataframeOrNull",
	sampleda="dataframeOrNull"),
    prototype=prototype(alpha=NULL,sampleda=NULL)
)

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

# #' @importClassesFrom phyloseq phylo
# #' @keywords internal 
# #' setClassUnion("phyloOrNULL", c("phylo", "NULL"))

#' @keywords internal
setClassUnion("listOrNull", c("list", "NULL"))

#' @keywords internal
setClassUnion("numericOrNull", c("numeric", "NULL"))

#' @keywords internal
setClassUnion("callOrNull", c("call", "function", "NULL"))

#' @title diffAnalysisClass class
#' @docType class
#' @slot originalD original feature data.frame.
#' @slot sampleda associated sample information.
#' @slot taxda the data.frame contained taxonomy.
#' @slot result data.frame contained the results of 
#' first, second test and LDA or rf 
#' @slot kwres the results of first test, contained
#' feature names, pvalue and fdr.
#' @slot secondvars the results of second test, contained
#' features names, gfc (TRUE representation the relevant 
#' feantures is enriched in relevant factorNames), 
#' Freq(the number of TRUE or FALSE), factorNames.
#' @slot mlres the results of LDA or randomForest,
#' @slot someparams, some arguments will be used in other functions
#' \code{\link[MicrobiotaProcess]{diff_analysis}}
#' @name diffAnalysisClass-class
#' @rdname diffAnalysisClass-class
#' @exportClass diffAnalysisClass
setClass("diffAnalysisClass",
    representation=representation(originalD="dataframeOrNull",
    sampleda="dataframeOrNull",
    taxda="dataframeOrNull",
    result="dataframeOrNull",
    kwres="dataframeOrNull",
    secondvars="listOrNull",
    mlres="dataframeOrNull",
    someparams="listOrNull"),
    prototype=prototype(originalD=NULL,
    sampleda=NULL,
    taxda=NULL,
    kwres=NULL,
    secondvars=NULL,
    mlres=NULL,
    someparams=NULL))

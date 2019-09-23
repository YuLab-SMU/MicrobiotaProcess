#' @title A box or violin plot with significance test
#' @param obj object, alphasample or data.frame (row sample x column features).
#' @param sampleda data.frame, sample information if obj is data.frame, the 
#' sampleda should be provided.
#' @param factorNames character, the names of factor contained in sampleda.
#' @param indexNames character, the vector character, should be the names of 
#' features contained object.
#' @param geom character, "boxplot" or "violin", default is "boxplot".
#' @param factorLevels list, the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this.
#' @param compare logical, whether test the features among groups,default is TRUE.
#' @param testmethod character, the method of test, default is `wilcox.test`.
#' see also \code{\link[ggsignif]{geom_signif}}.
#' @param signifmap logical, whether the pvalue are directly written a annotaion
#' or asterisks are used instead, default is (pvalue) FALSE. see also
#' \code{\link[ggsignif]{geom_signif}}.
#' @param p_textsize numeric, the size of text of pvalue or asterisks, 
#' default is 2.
#' @param step_increase numeric, see also \code{\link[ggsignif]{geom_signif}},
#' default is 0.1.
#' @param facetnrow integer, the nrow of facet, default is 1.
#' @param ... additional arguments, see also \code{\link[ggsignif]{geom_signif}}.
#' @export
#' @examples
#' library(tidyverse)
#' otudafile <- system.file("extdata", "otu_tax_table.txt",
#'                          package="MicrobiotaProcess")
#' otuda <- read.table(otudafile, sep="\t", 
#'                     header=TRUE, row.names=1,
#'                     check.names=FALSE, skip=1, 
#'                     comment.char="")
#' samplefile <- system.file("extdata",
#'                           "sample_info.txt", 
#'                           package="MicrobiotaProcess")
#' sampleda <- read.table(samplefile,
#'                        sep="\t", header=TRUE, row.names=1)
#' otuda <- otuda[sapply(otuda, is.numeric)] %>% t() %>%
#'          data.frame(check.names=FALSE)
#' set.seed(1024)
#' alphaobj1 <- get_alphaindex(otuda, sampleda=sampleda)
#' p1 <- ggbox(alphaobj1, factorNames="group")
#' data(test_otu_data)
#' set.seed(1024)
#' alphaobj2 <- get_alphaindex(test_otu_data)
#' class(alphaobj2)
#' head(as.data.frame(alphaobj2))
#' p2 <- ggbox(alphaobj2, factorNames="group")
setGeneric("ggbox", function(obj, factorNames, ...){standardGeneric("ggbox")})

#' @aliases ggbox,data.frame
#' @rdname ggbox
#' @importFrom ggplot2 ggplot geom_boxplot geom_violin aes_string facet_wrap 
#' @importFrom ggsignif geom_signif
#' @importFrom reshape melt
#' @importFrom dplyr filter
#' @export
setMethod("ggbox", "data.frame", 
          function(obj, sampleda, factorNames, indexNames, geom="boxplot",
                   factorLevels=NULL, compare=TRUE, testmethod="wilcox.test",
                   signifmap=FALSE, p_textsize=2, step_increase=0.1, 
                   facetnrow=1,...){
    if (missing(sampleda) || is.null(sampleda)){
        stop("the sampleda should be provided!")
    }
    if (missing(factorNames) || is.null(factorNames)){
        stop("the factorNames should be provided!")
    }
    sampleda <- sampleda[,match(factorNames, colnames(sampleda)), drop=FALSE]
    obj <- merge(obj, sampleda, by=0)
    rownames(obj) <- obj$Row.names
    obj$Row.names <- NULL
    obj <- melt(obj, id.vars=c(factorNames), variable_name="feature")
    if (!is.null(factorLevels)){
        obj <- setfactorlevels(obj, factorLevels)
    }
    if (missing(indexNames)||is.null(indexNames)){
        indexNames <- unique(as.vector(obj$feature))[seq_len(6)]
    }	
    obj <- obj %>% filter(.data$feature %in% indexNames)
    comparelist <- get_comparelist(obj, factorNames)
    mapping <- aes_string(x=factorNames,y="value",fill=factorNames)
    p <- ggplot(data=obj,mapping)
    ifelse(geom=="boxplot",p <- p + geom_boxplot(outlier.size=0.5,outlier.shape=21),
           p <- p + geom_violin())
    if (compare){
        p <- p + geom_signif(comparisons = comparelist,
			     test = testmethod,
                             map_signif_level = signifmap,
                             textsize = p_textsize, 
                             step_increase = step_increase,
                             ...)
    }
    p <- p + facet_wrap(~feature, scales="free", nrow=facetnrow) +
         theme_bw() + xlab(NULL) + ylab(NULL)
    return(p)
})

#' @aliases ggbox,alphasample
#' @rdname ggbox
#' @export
setMethod("ggbox", "alphasample", function(obj,factorNames,...){
    alphatab <- obj@alpha
    sampleda <- obj@sampleda
    p <- ggbox(obj=alphatab,
               sampleda=sampleda,
               factorNames=factorNames,
               ...)
    return(p)
})

#' @keywords internal
get_comparelist <- function(data, class){
    groups <- getclasslevels(data, class)
    comparelist <- getcompareclass(groups)
    comparelist <- split(comparelist, slice.index(comparelist, 1))
    names(comparelist) <- NULL
    return(comparelist)
}

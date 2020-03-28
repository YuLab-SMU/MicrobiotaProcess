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
#' see also \code{\link[ggsignif]{stat_signif}}.
#' @param signifmap logical, whether the pvalue are directly written a annotaion
#' or asterisks are used instead, default is (pvalue) FALSE. see also
#' \code{\link[ggsignif]{stat_signif}}.
#' @param p_textsize numeric, the size of text of pvalue or asterisks, 
#' default is 2.
#' @param step_increase numeric, see also \code{\link[ggsignif]{stat_signif}},
#' default is 0.1.
#' @param boxwidth numeric, the width of boxplot when the geom is 'violin',
#' default is 0.2.
#' @param facetnrow integer, the nrow of facet, default is 1.
#' @param ... additional arguments, see also \code{\link[ggsignif]{stat_signif}}.
#' @return a 'ggplot' plot object, a box or violine plot.
#' @author Shuangbin Xu
#' @export
#' @examples
#' library(magrittr)
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
                   signifmap=FALSE, p_textsize=2, step_increase=0.1, boxwidth=0.2,
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
        indexNames <- unique(as.vector(obj$feature))
    }
    obj <- obj %>% filter(.data$feature %in% indexNames)
    comparelist <- get_comparelist(data=obj, classgroup=factorNames)
    mapping <- aes_string(x=factorNames,y="value",fill=factorNames)
    p <- ggplot(data=obj,mapping)
    ifelse(geom=="boxplot",p <- p + geom_boxplot(outlier.size=0.5,outlier.shape=21),
           p <- p + geom_violin(trim=FALSE)+
                geom_boxplot(outlier.size=0.5,outlier.shape=21, 
			     width=boxwidth,
                             fill="white", show.legend=FALSE))
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


#' @title boxplot for the result of diff_analysis
#' @param obj object, diffAnalysisClass class.
#' @param geom character, "boxplot" or "violin", default is "boxplot".
#' @param box_notch logical, see also `notch` of \code{\link[ggplot2]{geom_boxplot}},
#' default is TRUE. 
#' @param dodge_width numeric, the width of dodge of boxplot, default is 0.6.
#' @param box_width numeric, the width of boxplot, default is 0.05
#' @param addLDA logical, whether add the plot to visulize the result of LDA, 
#' default is TRUE. 
#' @param factorLevels list, the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this.
#' @param featurelist vector, the character vector, the sub feature of
#' originalD in diffAnalysisClass,default is NULL.
#' @param removeUnkown logical, whether remove the unknow taxonomy,
#' default is TRUE.
#' @param colorlist character, the color vector, default is NULL.
#' @param l_xlabtext character, the x axis text of left panel,
#' default is NULL.
#' @param ..., additional arguments.
#' @return a 'ggplot' plot object, a box or violine plot for the result
#' of diffAnalysisClass.
#' @author Shuangbin Xu
#' @export
#' @examples
#' data(kostic2012crc)
#' kostic2012crc
#' head(phyloseq::sample_data(kostic2012crc),3)
#' kostic2012crc <- phyloseq::rarefy_even_depth(kostic2012crc,
#'                  rngseed=1024)
#' table(phyloseq::sample_data(kostic2012crc)$DIAGNOSIS)
#' set.seed(1024)
#' diffres <- diff_analysis(kostic2012crc, classgroup="DIAGNOSIS",
#'                          mlfun="lda", filtermod="fdr",
#'                          firstcomfun = "kruskal.test",
#'                          firstalpha=0.05, strictmod=TRUE,
#'                          secondcomfun = "wilcox.test",
#'                          subclmin=3, subclwilc=TRUE,
#'                          secondalpha=0.01, ldascore=3)
#' library(ggplot2)
#' p <- ggdiffbox(diffres, box_notch=FALSE, l_xlabtext="relative abundance")
setGeneric("ggdiffbox", function(obj, ...){standardGeneric("ggdiffbox")})

#' @aliases ggdiffbox,diffAnalysisClass
#' @rdname ggdiffbox
#' @importFrom ggplot2 coord_flip
#' @importFrom patchwork plot_layout
#' @export
setMethod("ggdiffbox", "diffAnalysisClass", function(obj, geom="boxplot", 
          box_notch=TRUE, box_width=0.05, dodge_width=0.6,
          addLDA=TRUE, factorLevels=NULL, featurelist=NULL,
          removeUnkown=TRUE, colorlist=NULL, l_xlabtext=NULL,...){
    featureda <- obj@originalD
    classname <- getcall(obj, "classgroup")
    normalization <- getcall(obj, "normalization")
    if (!is.null(normalization)){
        featureda <- featureda / normalization
    }
    sampleda <- obj@sampleda
    tmpgroup <- unique(as.vector(sampleda[[classname]]))
    nodedfres <- tidyEffectSize(obj) 
    nodedfres <- set_newlevels(data=nodedfres, newlevels=tmpgroup, factorNames=classname)
    if (is.null(featurelist)){
        featurelist <- unique(as.vector(nodedfres$f))
    }else{
        featurelist <- featurelist
    }
    featurelist <- keep_know_taxa(featurelist, removeUnkown=removeUnkown)
    nodedfres <- nodedfres[nodedfres$f%in%featurelist,,drop=FALSE]
    nodedfres <- nodedfres[order(nodedfres[,2], decreasing=TRUE),,drop=FALSE]
    nodedfres$f <- factor(nodedfres$f, levels=as.vector(nodedfres$f))
    featurelist <- as.vector(nodedfres$f)
    featureda <- featureda[,match(featurelist,colnames(featureda)),drop=FALSE]
    if (is.null(colorlist)){colorlist <- getCols(length(tmpgroup))}
    if (is.null(names(colorlist))){names(colorlist) <- tmpgroup}
    if(is.null(getcall(obj,"standard_method"))){
        ifelse(is.null(l_xlabtext), xlabtext<-"abundance", xlabtext <- l_xlabtext)
    }else{xlabtext<-"abundance"}
    p <- plotdiffbox(obj=featureda, sampleda=sampleda, factorNames=classname, factorLevels=factorLevels,
                     featurelist=featurelist, geom=geom, box_notch=box_notch,
                     dodge_width=dodge_width, box_width=box_width) + coord_flip() +
         scale_fill_manual(values=colorlist) + ylab(xlabtext) 
    if (addLDA){
        ifelse ("LDAmean" %in% colnames(obj@mlres), 
                effectsizename<-"LDA", effectsizename <- "MeanDecreaseAccuracy")
        colorlist <- colorlist[unique(as.vector(nodedfres[[classname]]))]
        p2 <- ggeffectsize.data.frame(obj=nodedfres, factorName=classname,
                                      effectsizename=effectsizename,
                                      setFacet=FALSE, ...) +
              scale_color_manual(values=colorlist) +
              theme(legend.position="none", axis.text.y=element_blank(),
		    plot.margin =unit(c(2,2,2,0),"mm"),
                    axis.ticks.y=element_blank(),
                    panel.grid=element_blank())
        p <- p + p2 + plot_layout(guides = 'collect', widths = c(3, 2))
    }
    return(p)
})

#' @keywords intermal
set_newlevels <- function(data, newlevels, factorNames){
    oldlevels <- unique(as.vector(data[[factorNames]]))
    newlevels <- oldlevels[match(newlevels, oldlevels)]
    newlevels <- newlevels[!is.na(newlevels)]
    data[[factorNames]] <- factor(data[[factorNames]], levels=newlevels)
    return(data)
}


#' @importFrom ggplot2 position_dodge
#' @keywords internal
plotdiffbox <- function(obj, sampleda, factorNames, dodge_width=0.6, box_width=0.05, 
                        geom="boxplot", factorLevels=NULL, featurelist=NULL, box_notch=TRUE){
    if (missing(sampleda) || is.null(sampleda)){
        stop("the sampleda should be provided!")
    }
    if (missing(factorNames) || is.null(factorNames)){
        stop("the factorNames should be provided!")
    }
    sampleda <- sampleda[, match(factorNames,colnames(sampleda)),drop=FALSE]
    featureda <- merge(obj, sampleda, by=0) %>%
                 column_to_rownames(var="Row.names")
    featureda <- melt(featureda, id.vars=c(factorNames), variable_name="feature")
    if (!is.null(featurelist)){
        featureda$feature <- factor(featureda$feature, levels=featurelist)
    }
    if (!is.null(factorLevels)){featureda <- setfactorlevels(featureda, factorLevels)}
    boxmapping <- aes_string(x="feature", y="value", color=factorNames)
    violinmapping <- modifyList(boxmapping,aes_string(color=NULL, fill=factorNames))
    dodge <- position_dodge(width = dodge_width)
    p <- ggplot()
    if (geom=="boxplot"){
       p <- p + geom_boxplot(data=featureda, notch=box_notch, mapping=violinmapping,
                             width=10*box_width, outlier.color=NA)
    }
    if (geom=="violin"){
        p <- p + geom_violin(data=featureda,
                             mapping=violinmapping,
			     position = dodge,
                             trim=FALSE) +
             geom_boxplot(data=featureda,
                          mapping=boxmapping,
                          outlier.color=NA,
                          width=box_width,
                          position = dodge,
                          notch=box_notch,
                          show.legend=FALSE)+
             scale_color_manual(values=rep("black", length(unique(featureda[[factorNames]]))))
    }
    p <- p + theme_bw() + xlab(NULL) + 
            theme(axis.text.x = element_text(size=12), plot.margin =unit(c(2,0,2,2),"mm"),
                  panel.grid=element_blank())
    return (p)
}

#' @keywords internal
get_comparelist <- function(data, classgroup){
    groups <- getclasslevels(sampleda=data, classgroup=classgroup)
    comparelist <- getcompareclass(classlevels=groups)
    comparelist <- split(comparelist, slice.index(comparelist, 1))
    names(comparelist) <- NULL
    return(comparelist)
}

#' @keywords internal
keep_know_taxa <- function(featurelist, removeUnkown=TRUE){
    if (isTRUE(removeUnkown)){
        tmpflag <- grep("__un_",featurelist)
        if (length(tmpflag)>0){
            featurelist <- featurelist[-tmpflag]
        }
    }
    return(featurelist)
}

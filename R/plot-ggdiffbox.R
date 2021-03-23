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
#' @param removeUnknown logical, whether remove the unknown taxonomy,
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
#' # set factor levels
#' #p2 <- ggdiffbox(diffres, box_notch=FALSE, l_xlabtext="relative abundance", 
#' #                factorLevels=list(DIAGNOSIS=c("Tumor", "Healthy")))
setGeneric("ggdiffbox", function(obj, ...){standardGeneric("ggdiffbox")})

#' @aliases ggdiffbox,diffAnalysisClass
#' @rdname ggdiffbox
#' @importFrom ggplot2 coord_flip
#' @importFrom patchwork plot_layout
#' @export
setMethod("ggdiffbox", "diffAnalysisClass", function(obj, geom="boxplot", 
          box_notch=TRUE, box_width=0.05, dodge_width=0.6,
          addLDA=TRUE, factorLevels=NULL, featurelist=NULL,
          removeUnknown=TRUE, colorlist=NULL, l_xlabtext=NULL, ...){
    featureda <- obj@originalD
    classname <- extract_args(obj, "classgroup")
    normalization <- extract_args(obj, "normalization")
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
    params <- list(...)
    if (!is.null(params$removeUnkown) && inherits(params$removeUnkown, "logical")){
        message("The `removeUnkown` has been deprecated, Please use `removeUnknown` instead!")
        removeUnknown <- params$removeUnkown
    }
    featurelist <- keep_know_taxa(featurelist, removeUnknown=removeUnknown)
    nodedfres <- nodedfres[nodedfres$f%in%featurelist,,drop=FALSE]
    nodedfres <- nodedfres[order(nodedfres[,2], decreasing=TRUE),,drop=FALSE]
    nodedfres$f <- factor(nodedfres$f, levels=as.vector(nodedfres$f))
    featurelist <- as.vector(nodedfres$f)
    featureda <- featureda[,match(featurelist,colnames(featureda)),drop=FALSE]
    if (is.null(colorlist)){colorlist <- get_cols(length(tmpgroup))}
    if (is.null(names(colorlist))){names(colorlist) <- tmpgroup}
    if(is.null(extract_args(obj,"standard_method"))){
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
keep_know_taxa <- function(featurelist, removeUnknown=TRUE){
    if (isTRUE(removeUnknown)){
        tmpflag <- grep("__un_",featurelist)
        if (length(tmpflag)>0){
            featurelist <- featurelist[-tmpflag]
        }
    }
    return(featurelist)
}

#' @title visualization of effect size by the Linear Discriminant Analysis or randomForest
#' @param obj object, diffAnalysisClass see \code{\link[MicrobiotaProcess]{diff_analysis}},
#' or data.frame, contained effect size and the group information.
#' @param removeUnknown logical, whether do not show unknown taxonomy, default is TRUE.
#' @param factorName character, the column name contained group information in data.frame.
#' @param effectsizename character, the column name contained effect size information.
#' @param factorLevels list, the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this.
#' @param linecolor character, the color of horizontal error bars, default is grey50.
#' @param linewidth numeric, the width of horizontal error bars, default is 0.4.
#' @param lineheight numeric, the height of horizontal error bars, default is 0.2.
#' @param pointsize numeric, the size of points, default is 1.5.
#' @param setFacet logical, whether use facet to plot, default is TRUE.
#' @param ... additional arguments.
#' @return the figures of effect size show the LDA or MDA (MeanDecreaseAccuracy).
#' @author Shuangbin Xu
#' @export
#' @examples
#' \dontrun{
#' data(kostic2012crc)
#' kostic2012crc
#' head(phyloseq::sample_data(kostic2012crc),3)
#' kostic2012crc <- phyloseq::rarefy_even_depth(kostic2012crc,rngseed=1024)
#' table(phyloseq::sample_data(kostic2012crc)$DIAGNOSIS)
#' set.seed(1024)
#' diffres <- diff_analysis(kostic2012crc, classgroup="DIAGNOSIS",
#'                         mlfun="lda", filtermod="fdr",
#'                         firstcomfun = "kruskal.test",
#'                         firstalpha=0.05, strictmod=TRUE,
#'                         secondcomfun = "wilcox.test", 
#'                         subclmin=3, subclwilc=TRUE,
#'                         secondalpha=0.01, ldascore=3) 
#' library(ggplot2)
#' effectplot <- ggeffectsize(diffres) +
#'               scale_color_manual(values=c('#00AED7', 
#'                                           '#FD9347', 
#'                                           '#C1E168'))+
#'               theme_bw()+
#'               theme(strip.background=element_rect(fill=NA),
#'                     panel.spacing = unit(0.2, "mm"),
#'                     panel.grid=element_blank(),
#'                     strip.text.y=element_blank())
#' }
ggeffectsize <- function(obj,...){
    UseMethod("ggeffectsize")
}

#' @method ggeffectsize data.frame
#' @rdname ggeffectsize
#' @importFrom ggplot2 facet_grid ylab xlab scale_color_manual theme_bw element_text element_blank unit element_rect
#' @importFrom ggplot2 geom_errorbarh geom_point scale_x_continuous
#' @export 
ggeffectsize.data.frame <- function(obj, 
    factorName, 
    effectsizename,
    factorLevels=NULL,
    linecolor="grey50",
    linewidth=0.4,
    lineheight=0.2,
    pointsize=1.5,
    setFacet=TRUE,
    ...){
    if (effectsizename %in% "LDA"){
    	xlabtext <- bquote(paste(Log[10],"(",.("LDA"), ")"))
        xtext <- "LDAmean"
        xmintext <- "LDAlower"
        xmaxtext <- "LDAupper"
    }else{
    	xlabtext <- effectsizename
        xtext <- "MDAmean"
        xmintext <- "MDAlower"
        xmaxtext <- "MDAupper"
    }
    if (!is.null(factorLevels)){
    	obj <- setfactorlevels(obj,factorLevels)
    }
    #p <- ggplot(data=obj, 
    # 	        aes_(x=as.formula(paste0("~",effectsizename)),
    # 		     y=~f)) + 
    # 	geom_segment(aes_(xend=0, yend=~f), 
    #		     color="grey") + 
    p <- ggplot(data=obj, aes_(y=~f)) +
         geom_errorbarh(aes_(xmin=as.formula(paste0("~", xmintext)), 
                             xmax=as.formula(paste0("~", xmaxtext))),
                        color=linecolor, size=linewidth, height=lineheight) +
         geom_point(aes_(x=as.formula(paste0("~",xtext)),
                         color=as.formula(paste0("~",factorName))),
                    size=pointsize) 
    if (setFacet) {
        p <- p + facet_grid(as.formula(paste0(factorName," ~.")),
    			   scales = "free_y", space = "free_y")
    }
    p <- p + #scale_x_continuous(expand=c(0,0), 
        #limits=c(0, max(obj[[xmaxtext]])*1.1))+
    	ylab(NULL) + xlab(xlabtext) 
    #if (setColors){
    tmpn <- length(unique(as.vector(obj[[factorName]])))
    message("The color has been set automatically, you can reset it manually by adding scale_color_manual(values=yourcolors)")
    p <- p + scale_color_manual(values=get_cols(tmpn))
    #}
    p <- p + theme_bw()+
    	theme(axis.text.y = element_text(size=12),
	          axis.text.x = element_text(size=12),
      	      panel.grid=element_blank(),
    	      #panel.spacing = unit(0.2, "mm"),
    	      legend.title=element_text(size=10),
    	      legend.text=element_text(size=8),
    	      strip.background = element_rect(colour=NA, 
    				   fill="grey"))
    return(p)
}

#' @method ggeffectsize diffAnalysisClass
#' @rdname ggeffectsize
#' @export
ggeffectsize.diffAnalysisClass <- function(obj, removeUnknown=TRUE, setFacet=TRUE,...){
    #efres <- tidyEffectSize(obj)
    efres <- obj@result
    classname <- extract_args(obj, "classgroup")
    params <- list(...)
    if (!is.null(params$removeUnkown) && inherits(params$removeUnkown, "logical")){
        message("The `removeUnkown` has been deprecated, Please use `removeUnknown` instead!")
        removeUnknown <- params$removeUnkown
    }
    if (removeUnknown && length(grep("__un_",efres$f))){
    	efres <- efres[-grep("__un_",efres$f),,drop=FALSE]
    }
    if ("LDAmean" %in% colnames(obj@mlres)){
    	effectsizename <- "LDA"
    }else{
    	effectsizename <- "MeanDecreaseAccuracy" 
    }
    p <- ggeffectsize.data.frame(obj=efres, 
    			         factorName=classname, 
    			         effectsizename=effectsizename,
                                 setFacet=setFacet,
    			          ...)
    return (p)
}


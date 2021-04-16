#' @title plot the clade tree with highlight
#' 
#' @description plot results of different analysis or data.frame, 
#' contained hierarchical relationship or other classes,such like the 
#' tax_data of phyloseq.
#' @param obj object, diffAnalysisClass, the results of diff_analysis 
#' see also \code{\link[MicrobiotaProcess]{diff_analysis}}, or data.frame,
#' contained hierarchical relationship or other classes.
#' @param nodedf data.frame, contained the tax and the factor 
#' information and(or pvalue).
#' @param factorName character, the names of factor in nodedf.
#' @param removeUnknown logical, whether do not show unknown taxonomy, 
#' default is TRUE.
#' @param layout character, the layout of ggtree, but only "rectangular",
#' "roundrect", "ellipse", "radial", "slanted", "inward_circular" and 
#' "circular" in here, default is "radial".
#' @param linewd numeric, the size of segment of ggtree, default is 0.6.
#' @param skpointsize numeric, the point size of skeleton of tree, 
#' default is 0.8 .
#' @param alpha numeric, the alpha of clade, default is 0.4.
#' @param taxlevel positive integer, the full text of clade, default is 5.
#' @param cladetext numeric, the size of text of clade, default is 2.
#' @param factorLevels list, the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this.
#' @param setColors logical, whether set the color of clade, default is TRUE, 
#' or set FALSE,then use 'scale_fill_manual' setting.
#' @param xlim numeric, the x limits, only works for 'inward_circular' layout,
#' default is 12.
#' @param reduce logical, whether remove the unassigned taxonomy, which will 
#' remove the clade of unassigned taxonomy, but the result of `diff_analysis`
#' should remove the unknown taxonomy, default is FALSE. 
#' @param type character, the type of datasets, default is "species", 
#' if the dataset is not about species, such as dataset of kegg function, 
#' you should set it to "others".
#' @param ..., additional parameters.
#' @return figures of tax clade show the significant different feature.
#' @author Shuangbin Xu
#' @export
#' @examples
#' data(kostic2012crc)
#' kostic2012crc
#' head(phyloseq::sample_data(kostic2012crc),3)
#' kostic2012crc <- phyloseq::rarefy_even_depth(kostic2012crc,
#'                          rngseed=1024)
#' table(phyloseq::sample_data(kostic2012crc)$DIAGNOSIS)
#' set.seed(1024)
#' diffres <- diff_analysis(kostic2012crc, classgroup="DIAGNOSIS",
#'                         mlfun="lda", filtermod="fdr",
#'                         firstcomfun = "kruskal.test",
#'                         firstalpha=0.05, strictmod=TRUE,
#'                         secondcomfun = "wilcox.test",
#'                         subclmin=3, subclwilc=TRUE,
#'                         secondalpha=0.01, ldascore=3)
#' #library(ggplot2)
#' #diffcladeplot <- ggdiffclade(diffres,alpha=0.3, linewd=0.2, 
#' #                        skpointsize=0.4, 
#' #                        taxlevel=3,
#' #                        setColors=FALSE) +
#' #        scale_fill_manual(values=c('#00AED7', 
#' #                                   '#FD9347', 
#' #                                   '#C1E168'))
ggdiffclade <- function(obj,...){
    UseMethod("ggdiffclade")
}

#' @method ggdiffclade data.frame
#' @rdname ggdiffclade
#' @importFrom ggplot2 geom_point geom_rect geom_text scale_size_continuous scale_fill_manual 
#' @importFrom magrittr %<>%
#' @importFrom stats as.formula
#' @export
ggdiffclade.data.frame <- function(obj, nodedf, factorName, layout="radial", linewd=0.6, 
    skpointsize=0.8, alpha=0.4, taxlevel=6, cladetext=2, factorLevels=NULL, setColors=TRUE,
    xlim=12, reduce=FALSE, type="species", ...){
    params <- list(...)
    if (!is.null(params$size)){
        message("The `size` has been deprecated, Please use `linewd` instead!")
        linewd <- params$size
    }
    treedata <- convert_to_treedata(obj, type=type)
    layout %<>% match.arg(c("rectangular", "roundrect", "ellipse", "circular", "slanted", "radial", "inward_circular"))
    if (!is.null(factorLevels)){nodedf <- setfactorlevels(nodedf, factorLevels)}
    p <- treeskeleton(treedata, layout=layout, size=linewd, pointsize=skpointsize, xlim=xlim)
    nodedf <- get_node(treedata=p$data, nodedf=nodedf)
    if (reduce){
        df <- p$data[!grepl("__un_",p$data$label),]
        nodedf <- nodedf[!grepl("__un_", as.vector(nodedf[,1])),,drop=FALSE]
        p <- treeskeleton(treedata=df, layout=layout,size=linewd, pointsize=skpointsize, xlim=xlim) 
    }
    cladecoord <- get_cladedf(p, nodedf$node)
    cladecoord <- merge(cladecoord, nodedf, by.x="node", by.y="node")
    if (layout %in% c("rectangular","slanted", "roundrect", "ellipse")){
        #taxlevel <- 7
        labelannotcoord <- get_labeldf(p, nodedf$node, angle=90)
    }else{
        labelannotcoord <- get_labeldf(p, nodedf$node)
    }
    labelannotcoord <- get_annotlabel(labelannotcoord, classlevel=taxlevel)
    labelcoord <- labelannotcoord$labeldf
    annotcoord <- labelannotcoord$annotdf
    if ("pvalue" %in% colnames(nodedf)){
        nodedf[["-log10(pvalue)"]] <- -log10(nodedf$pvalue)
        pointmapping <- aes_(x=~x,y=~y,fill=as.formula(paste("~",factorName)), size=~-log10(pvalue))
    }else{
        pointmapping <- aes_(x=~x,y=~y, fill=as.formula(paste("~",factorName)))
    }
    p <- p + geom_rect(data=cladecoord,
                       aes_(fill=as.formula(paste("~",factorName)), 
                       xmin=~xmin, ymin=~ymin, xmax=~xmax, ymax=~ymax),
                       alpha=alpha,
    	               show.legend=FALSE) +
        geom_text(data=labelcoord, aes_(x=~x,y=~y, label=~label, angle=~angle), size=cladetext) +
        geom_point(data=cladecoord, pointmapping, shape=21)+
        geom_point(data=annotcoord, aes_(x=0,y=0, color=~label), size=0, stroke=0) +
        scale_size_continuous(range = c(1, 3))
    if (setColors){
        message("The color has been set automatically, you can reset it manually by adding scale_fill_manual(values=yourcolors)")
        tmpn <- length(unique(as.vector(nodedf[[match(factorName,colnames(nodedf))]])))
        p <- p + scale_fill_manual(values=get_cols(tmpn))
    }
    p <- p + theme_diffclade() + guides_diffclade() 
    return(p)
}


#' @method ggdiffclade diffAnalysisClass
#' @rdname ggdiffclade
#' @importFrom magrittr %>%
#' @export
ggdiffclade.diffAnalysisClass <- function(obj, removeUnknown=TRUE, ...){
    params <- list(...)
    if (!is.null(params$removeUnkown) && inherits(params$removeUnkown, "logical")){
        message("The `removeUnkown` has been deprecated, Please use `removeUnknown` instead!")
        removeUnknown <- params$removeUnkown
    }
    taxda <- obj@taxda
    #secondvars <- get_second_true_var(obj)
    #nodedfres <- merge(obj@kwres, secondvars, by.x="f", by.y="f") %>% 
    #			select(-c("gfc", "Freq"))
    #nodedfres <- nodedfres[as.vector(nodedfres$f)%in%as.vector(obj@mlres$f),]
    #nodedfres <- merge(nodedfres, obj@mlres$f, by.x="f", by.y="f")
    nodedfres <- as.data.frame(obj)
    classname <- extract_args(obj, "classgroup")
    if (removeUnknown){
        tmpflag <- grep("__un_",as.vector(nodedfres$f))
        if (length(tmpflag)>0){
   	        nodedfres <- nodedfres[-tmpflag,,drop=FALSE]
        }
    }
    p <- ggdiffclade(obj=taxda,
    			     nodedf=nodedfres,
    			     factorName=classname,
                     type=extract_args(obj, "type"),
    			     ...)
    return(p)
}
#' @importFrom dplyr rename
#' @importFrom rlang .data
#' @keywords internal
get_node <- function(treedata, nodedf){
    nodelist <- treedata[match(as.vector(nodedf[,1]),treedata$label),]$node
    if (!"node" %in% colnames(nodedf)){
    	nodedf$node <- nodelist
    }else{
    	nodedf %>% rename(nodetmp=.data$node)
    	nodedf$node <- nodelist
    }
    return(nodedf)
}


#' @importFrom ggtree ggtree
#' @importFrom ggtree geom_point
#' @keywords internal
treeskeleton <- function(treedata, layout, size, pointsize=1, xlim=12){
    if (layout=="inward_circular"){
        p <- ggtree(
                 treedata, 
                 layout=layout, 
                 size=size,
                 xlim=c(xlim,NA)
                 )
    
    }else{
        p <- ggtree(
                 treedata,
    	         layout=layout,
    	         size=size
                 )
    }
    p <- p + 
         geom_point(
             size=pointsize,
             shape=21, 
             fill="white",
             color="black"
         )
    return(p)
}

#' @importFrom ggplot2 guides guide_legend
#' @keywords internal
guides_diffclade <- function(...){
    guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5, order=1),
    	   size=guide_legend(keywidth = 0.5, keyheight = 0.5, order=2),
    	   color = guide_legend(keywidth = 0.1, ncol=1, keyheight = 0.6, order = 3),...)
}

#' @importFrom ggplot2 theme unit element_text element_rect margin
#' @keywords internal
theme_diffclade <- function(...){
    theme(legend.position="right",
    	  legend.margin=margin(0,0,0,0),
    	  legend.spacing.y = unit(0.02, "cm"),
    	  legend.box.spacing=unit(0.02,"cm"),
    	  legend.text = element_text(size=5),
    	  legend.title=element_text(size=6),
    	  legend.background=element_rect(fill=NA),
		  ...)
}

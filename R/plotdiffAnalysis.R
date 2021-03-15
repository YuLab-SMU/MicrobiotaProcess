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
#' "radial", "slanted", "inward_circular" and "circular" in here, 
#' default is circular.
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
ggdiffclade.data.frame <- function(obj, nodedf, factorName, layout="circular", linewd=0.6, 
    skpointsize=0.8, alpha=0.4, taxlevel=6, cladetext=2, factorLevels=NULL, setColors=TRUE,
    xlim=12, reduce=FALSE, type="species", ...){
    params <- list(...)
    if (!is.null(params$size)){
        message("The `size` has been deprecated, Please use `linewd` instead!")
        linewd <- params$size
    }
    treedata <- convert_to_treedata(obj, type=type)
    layout %<>% match.arg(c("rectangular", "circular", "slanted", "radial", "inward_circular"))
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
    if (layout %in% c("rectangular","slanted")){
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
    p <- p + geom_rect(data=cladecoord,aes_(fill=as.formula(paste("~",factorName)), 
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


#' @title significantly discriminative feature barplot
#' @param obj object, diffAnalysisClass see also 
#' \code{\link[MicrobiotaProcess]{diff_analysis}} or feMeanMedian class, 
#' see also \code{\link[MicrobiotaProcess]{get_mean_median}}.
#' @param filepath character, default is NULL, meaning current path. 
#' @param output character, the output dir name, default is "biomarker_barplot".
#' @param xtextsize numeric, the size of axis x label, default is 3.
#' @param removeUnknown logical, whether do not show unknown taxonomy, default is TRUE.
#' @param figwidth numeric, the width of figures, default is 6.
#' @param figheight numeric, the height of figures, default is 3.
#' @param ylabel character, the label of y, default is 'relative abundance'.
#' @param featurename character, the feature name, contained at the objet.
#' @param classgroup character, factor name.
#' @param subclass character, factor name. 
#' @param factorLevels list,  the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this. 
#' @param format character, the format of figure, default is pdf,
#' png, tiff also be supported.
#' @param coloslist vector, color vector, if the input is phyloseq, 
#' you should use this to adjust the color, not scale_color_manual.
#' @param dpi numeric, the dpi of output, default is 300.
#' @param ... additional arguments.
#' @return the figures of features show the distributions in samples.
#' @author Shuangbin Xu
#' @export
#' @examples
#' data(kostic2012crc)
#' kostic2012crc
#' head(phyloseq::sample_data(kostic2012crc),3)
#' kostic2012crc <- phyloseq::rarefy_even_depth(kostic2012crc,
#'                               rngseed=1024)
#' table(phyloseq::sample_data(kostic2012crc)$DIAGNOSIS)
#' #set.seed(1024)
#' #diffres <- diff_analysis(kostic2012crc, classgroup="DIAGNOSIS",
#' #                        mlfun="lda", filtermod="fdr",
#' #                        firstcomfun = "kruskal.test",
#' #                        firstalpha=0.05, strictmod=TRUE,
#' #                        secondcomfun = "wilcox.test",
#' #                        subclmin=3, subclwilc=TRUE,
#' #                        secondalpha=0.01, ldascore=3)
#' # not run in example
#' #ggdifftaxbar(diffres, output="biomarker_barplot")
ggdifftaxbar <- function(obj,...){
    UseMethod("ggdifftaxbar")
}

#' @keywords internal
setGeneric("ggdifftaxbar")

#' @aliases ggdifftaxbar,diffAnalysisClass
#' @rdname ggdifftaxbar
#' @importFrom ggplot2 ggsave
#' @export
setMethod("ggdifftaxbar","diffAnalysisClass",function(obj,
    filepath=NULL,
    output="biomarker_barplot",
    removeUnknown=TRUE,
    figwidth=6,
    figheight=3,
    ylabel="relative abundance",
    format="pdf",
    dpi=300,
    ...){
    featureda <- obj@originalD
    classname <- extract_args(obj, "classgroup")
    normalization <- extract_args(obj, "normalization")
    if (!is.null(normalization)){
        featureda <- featureda / normalization
    }
    sampleda <- obj@sampleda
    featureda <- merge(featureda, sampleda, by=0) %>%
        column_to_rownames(var="Row.names")
    nodedfres <- as.data.frame(obj)
    featurelist <- as.vector(nodedfres$f)
    params <- list(...)
    if (!is.null(params$removeUnkown) && inherits(params$removeUnkown, "logical")){
        message("The `removeUnkown` has been deprecated, Please use `removeUnknown` instead!")
        removeUnknown <- params$removeUnkown
    }
    if (removeUnknown && length(grep("__un_", featurelist))>0){
        featurelist <- featurelist[!grepl("__un_", featurelist)]
    }
    if (ncol(sampleda)>1){
        subclass <- colnames(sampleda)[-match(classname, colnames(sampleda))]
    }else{
        subclass <- classname
    }
    if(is.null(filepath)){filepath <- getwd()}
    filepath <- file.path(filepath, output)
    dir.create(filepath, showWarnings = FALSE)
    for (vars in featurelist){
        resdf <- get_mean_median(
                      datameta=featureda, 
                      feature=vars, 
                      subclass=subclass
                      )
        p <- ggdifftaxbar.featureMeanMedian(
                  obj=resdf,
                  featurename=vars,
                  classgroup=classname,
                  subclass=subclass,
                  ylabel=ylabel,
                  ...
                  )
        if (grepl("/", vars)){
            vars <- sub("/", "--", vars)
        }
        filename <- file.path(filepath, paste(vars, format, sep="."))
        ggsave(filename, p, device=format, width = figwidth, height=figheight, dpi=dpi, limitsize = FALSE)
    }
})

#' @method ggdifftaxbar featureMeanMedian
#' @rdname ggdifftaxbar
#' @importFrom ggplot2 aes geom_errorbar scale_linetype_manual unit theme_bw scale_y_continuous labs xlab facet_grid
#' @importFrom ggplot2 guide_legend element_text element_rect element_blank scale_fill_manual
#' @export
ggdifftaxbar.featureMeanMedian <- function(obj, featurename, classgroup, subclass, xtextsize=3,
    factorLevels=NULL, coloslist=NULL, ylabel="relative abundance", ...){
    data <- obj$singlefedf
    dastatistic <- obj$singlefestat
    if (!is.null(factorLevels)){
    	data <- setfactorlevels(data, factorLevels)
    	dastatistic <- setfactorlevels(dastatistic, factorLevels)
    }
    if (missing(subclass)){subclass <- classgroup}
    p <- ggplot(data, aes_(x=~sample, y=~RelativeAbundance, fill=as.formula(paste0("~",subclass))))+
    	geom_bar(stat="identity") +
    	geom_errorbar(data=dastatistic, aes_(x=~sample, ymax=~value, ymin=~value, linetype=~statistic), 
	              size=0.5, width=1, inherit.aes=FALSE)+
    	scale_linetype_manual(values=c("solid", "dotted"))+
    	facet_grid(as.formula(paste0("~",classgroup)), space="free_x", scales="free_x") + 
    	labs(title=featurename) + xlab(NULL)+ ylab(ylabel)+
    	scale_y_continuous(expand=c(0,0), limits=c(0,max(data$RelativeAbundance)*1.05))
    p <- p + theme_bw() + guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5, order=1),
    			   linetype=guide_legend(keywidth = 0.7, keyheight = 0.5, order=2))+
    	 theme(plot.title = element_text(face="bold",lineheight=25,hjust=0.5), legend.box.spacing=unit(0.02,"cm"),
    	       panel.grid=element_blank(), legend.text = element_text(size=6.5),  legend.title=element_text(size=7),
    	       legend.background=element_rect(fill=NA), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
    	       panel.spacing = unit(0.2, "mm"), strip.background = element_rect(colour=NA,fill="grey"))
    if (is.null(coloslist)){
        tmpn <- length(unique(as.vector(data[[match(subclass,colnames(data))]])))
        p <- p + scale_fill_manual(values=get_cols(tmpn))
    }else{
        p <- p + scale_fill_manual(values=coloslist)
    }
    return(p)
}


#' @title get the mean and median of specific feature.
#' @param datameta data.frame, nrow sample * ncol feature + factor.
#' @param feature character vector, the feature contained in datameta.
#' @param subclass character, factor name.
#' @importFrom dplyr rename group_by mutate
#' @importFrom magrittr %>%
#' @importFrom stats median
#' @return featureMeanMedian object, contained the abundance of feature, and the 
#' mean and median of feature by subclass.
#' @author Shuangbin Xu
#' @export
#' @examples
#' data(hmp_aerobiosis_small)
#' head(sampleda)
#' featureda <- merge(featureda, sampleda, by=0)
#' rownames(featureda) <- as.vector(featureda$Row.names)
#' featureda$Row.names <- NULL
#' feameamed <- get_mean_median(datameta=featureda, 
#'                     feature="p__Actinobacteria", 
#'                     subclass="body_site")
#' #not run in example
#' #fplot <- ggdifftaxbar(feameamed, featurename="p__Actinobacteria", 
#' #                     classgroup="oxygen_availability", subclass="body_site")
get_mean_median <- function(datameta, feature, subclass){
    RelativeAbundance <- NULL
    factornames <- colnames(datameta)[!unlist(vapply(datameta,is.numeric,logical(1)))]
    featuredatmp <- datameta %>% rownames_to_column(var="sample") %>%
    		select(c("sample", feature, factornames)) %>%
    		rename(RelativeAbundance=feature) %>% 
    		mutate(sample = factor(sample, levels=sample[order(eval(parse(text=subclass)), -RelativeAbundance)]))
    meantmp <- featuredatmp %>% group_by(eval(parse(text=subclass))) %>% 
    		mutate(value = mean(RelativeAbundance)) %>% mutate(statistic="mean")
    mediantmp <- featuredatmp %>% group_by(eval(parse(text=subclass))) %>% 
    		mutate(value = median(RelativeAbundance)) %>% mutate(statistic="median")
    festatic <- rbind(meantmp, mediantmp) %>% data.frame(check.names=FALSE)
    res <- structure(list(singlefedf=featuredatmp, singlefestat=festatic), 
    				 class="featureMeanMedian")
    return (res)
}

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
#' data(kostic2012crc)
#' kostic2012crc
#' head(phyloseq::sample_data(kostic2012crc),3)
#' kostic2012crc <- phyloseq::rarefy_even_depth(kostic2012crc,rngseed=1024)
#' table(phyloseq::sample_data(kostic2012crc)$DIAGNOSIS)
#' #set.seed(1024)
#' #diffres <- diff_analysis(kostic2012crc, classgroup="DIAGNOSIS",
#' #                        mlfun="lda", filtermod="fdr",
#' #                        firstcomfun = "kruskal.test",
#' #                        firstalpha=0.05, strictmod=TRUE,
#' #                        secondcomfun = "wilcox.test", 
#' #                        subclmin=3, subclwilc=TRUE,
#' #                        secondalpha=0.01, ldascore=3) 
#' #library(ggplot2)
#' #effectplot <- ggeffectsize(diffres) +
#' #              scale_color_manual(values=c('#00AED7', 
#' #                                          '#FD9347', 
#' #                                          '#C1E168'))+
#' #              theme_bw()+
#' #              theme(strip.background=element_rect(fill=NA),
#' #                    panel.spacing = unit(0.2, "mm"),
#' #                    panel.grid=element_blank(),
#' #                    strip.text.y=element_blank())
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
    efres <- tidyEffectSize(obj)
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

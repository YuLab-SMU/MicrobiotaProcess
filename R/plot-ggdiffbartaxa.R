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

#' @rdname ggdifftaxbar
#' @export
ggdiffbartaxa <- ggdifftaxbar

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

#' @title taxonomy barplot
#' @param obj phyloseq, phyloseq class or data.frame, 
#' (nrow sample * ncol feature (factor)) or the data.frame for geom_bar.
#' @param mapping set of aesthetic mapping of ggplot2, default is NULL,
#' if the data is the data.frame for geom_bar, the mapping should be set.
#' @param position character, default is `stack`. 
#' @param stat character, default is `identity`.
#' @param width numeric, the width of bar, default is 0.7.
#' @param topn integer, the top number of abundance taxonomy(feature).
#' @param count logical, whether show the relative abundance.  
#' @param sampleda data.frame, (nrow sample * ncol factor), the sample 
#' information, if the data doesn't contain the information.
#' @param factorLevels vector or list, the levels of the factors
#' (contained names e.g. list(group=c("B","A","C")) or
#' c(group=c("B","A","C"))), adjust the order of facet, default is NULL,
#' if you want to order the levels of factor, you can set this.
#' @param sampleLevels vector, adjust the order of x axis
#' e.g. c("sample2", "sample4", "sample3"), default is NULL.
#' @param facetNames character, default is NULL.
#' @param plotgroup logical, whether calculate the mean or median etc 
#' for each group, default is FALSE.
#' @param groupfun character, how to calculate for feature in each group,
#' the default is `mean`, this will plot the mean of feature in each group.
#' @param ... additional parameters, see \code{\link[ggplot2]{ggplot}}
#' @return barplot of tax
#' @author Shuangbin Xu
#' @export
#' @examples
#' library(ggplot2)
#' data(test_otu_data)
#' otubar <- ggbartax(test_otu_data) + 
#'          xlab(NULL) + ylab("relative abundance(%)")
ggbartax <- function(obj,...){
    UseMethod("ggbartax")
}

#' @rdname ggbartax
#' @export
ggbartaxa <- ggbartax

#' @method ggbartax phyloseq
#' @importFrom phyloseq otu_table taxa_are_rows
#' @rdname ggbartax
#' @export
ggbartax.phyloseq <- function(obj, ...){
    if (is.null(obj@otu_table)){
    	stop("The otu table is empty!")
    }else{
    	otudata <- get_otudata(obj)
    }
    if (!is.null(obj@sam_data)){
    	sampleda <- data.frame(obj@sam_data, check.names=FALSE)
    	p <- ggbartax(obj=otudata, sampleda=sampleda, ...)
    }else{
    	p <- ggbartax(obj=otudata,...)
    }
    return(p)	
}
#' @method ggbartax data.frame
#' @rdname ggbartax
#' @importFrom ggplot2 ggplot aes_ geom_bar scale_fill_manual facet_grid
#' @importFrom stats as.formula
#' @importFrom rlang as_name
#' @importFrom stats aggregate
#' @export
ggbartax.data.frame <- function(obj, mapping=NULL, position = "stack", stat="identity",
    width=0.7, topn=30, count=FALSE, sampleda=NULL, factorLevels=NULL, sampleLevels=NULL,
    facetNames=NULL, plotgroup=FALSE, groupfun=mean,...){
    if (is.null(mapping)){
        mapping <- aes_(x=~sample, y=~value, fill=~feature)
        obj <- mappingtaxda(data=obj, topn=topn, count=count, sampleda=sampleda, 
                            factorLevels=factorLevels, sampleLevels=sampleLevels)
        if(!is.null(facetNames) & plotgroup){
            formulatmp <- as.formula(paste0("value ~ ", 
                                     paste(c(as.name(facetNames), "feature"),collapse="+")))
            obj <- aggregate(formulatmp, obj, FUN=groupfun)
            mapping <- aes_string(x=facetNames, y="value", fill="feature")
        }
    }else{
        mapping <- mapping
    }
    p <- ggplot(data=obj,mapping=mapping,...) + 
         geom_bar(position = position,stat=stat, width=width) + 
         scale_y_continuous(expand=c(0,0))
    if ("fill" %in% names(mapping)){
        tmpn <- length(unique(as.vector(obj[,as_name(mapping[["fill"]])])))
    }else{
        tmpn <- nrow(obj)
    }
    message("The color has been set automatically, you can reset it 
            manually by adding scale_fill_manual(values=yourcolors)")
    p <- p + scale_fill_manual(values=get_cols(tmpn))
    if (!is.null(facetNames) & !plotgroup){
        tmpformula <- as.formula(paste0("~ ", facetNames))
        p <- p + facet_grid(tmpformula, scales="free_x", space="free_x")
    }
    p <- p + theme_taxbar()
    return(p)
}

#' @title generate the mapping data
#' @description
#' build the data.frame for `geom_bar` of `ggplot2`,
#' if the data.frame was built, you can visulize with `ggbartax`.
#' @param data data.frame, (nrow sample * ncol feature (factor))
#' @param topn integer, the top number of abundance taxonomy(feature).
#' @param count boolean, default is FALSE.
#' @param sampleda data.frame, (nrow sample * ncol factor) the 
#' sample information, if the data doesn't contain the information.
#' @param factorLevels vector, the levels of the factors (contained names e.g. 
#' list(group=c("B","A","C")) or c(group=c("B","A","C"))), adjust the order of
#' facet, default is NULL, if you want to order the levels of factor, you can set this.
#' @param sampleLevels vector, the levels of samples, adjust the order of x axis
#' e.g. c("sample2", "sample4"), default is NULL.
#' @param plotda boolean, default is TRUE, whether build the data.frame for
#' `geom_bar` of `ggplot2`.
#' @return the data.frame for ggbartax
#' @author Shuangbin Xu
#' @importFrom magrittr %>%
#' @keywords internal
#' @noRd
mappingtaxda <- function(data, topn=30, count=FALSE, sampleda=NULL, 
    factorLevels=NULL, sampleLevels=NULL, plotda=TRUE){
    tmpfeature <- colnames(data)[vapply(data,is.numeric,logical(1))]
    tmpfactor <- colnames(data)[!vapply(data,is.numeric,logical(1))]
    dat <- data[, tmpfeature,drop=FALSE] %>% t() %>% data.frame(check.names=FALSE)
    if(!count){dat <- apply(dat, 2, function(x){100*x/sum(x)}) %>% data.frame(check.names=FALSE)}
    dat$sums <- apply(dat, 1, sum)
    dat <- dat[order(dat$sums, decreasing = TRUE),,drop=FALSE]
    dat$sums <- NULL
    tmpsums <- matrix(colSums(dat),nrow=1) %>% data.frame()
    if (topn < nrow(dat)){
    	dat <- dat[seq_len(topn),,drop=FALSE]
    	if (!count){others <- 100 - (matrix(apply(dat,2,sum),nrow=1) %>% data.frame(check.names=FALSE))}
    	if (count){others <- tmpsums - (matrix(apply(dat,2,sum),nrow=1) %>% data.frame(check.names=FALSE))}
    	colnames(others) <- colnames(dat)
    	rownames(others) <- "Others"
    	dat <- rbind(dat, others)
    }
    featurelevels <- rownames(dat)
    if (plotda){
    	dat <- dat %>% 
               tibble::as_tibble(rownames="feature") %>% 
               tidyr::pivot_longer(!c("feature"), names_to="sample", values_to="value")
    	if (!is.null(sampleda)){
    		sampleda$sample <- rownames(sampleda)
    		dat <- merge(dat, sampleda)
    	}
    	if (is.null(sampleda) && length(tmpfactor)>0){
    		tmpsample <- data[,tmpfactor,drop=FALSE]
    		tmpsample$sample <- rownames(tmpsample)
    		dat <- merge(dat, tmpsample)
    	}
    	dat$feature <- factor(dat$feature, levels=featurelevels)
    }else{
    	if (!is.null(sampleda)){
    		dat <- dat %>% t() %>% data.frame(check.names=FALSE)
    		dat <- merge(dat, sampleda, by=0)
    	}
    	if (is.null(sampleda)&&length(tmpfactor)>0){
    		tmpsample <- data[,tmpfactor,drop=FALSE]
    		dat <- merge(dat, tmpsample, by=0)
    	}
    	colnames(dat)[1] <- "sample"		
    }
    if (is.numeric(dat$sample)){
       dat$sample <- as.character(dat$sample)
    }
    if (!is.null(factorLevels)){dat <- setfactorlevels(dat, factorLevels)}
    if (!is.null(sampleLevels)){dat$sample <- factor(dat$sample, levels=sampleLevels)}
    return(dat)
}


#setfactorlevels <- function(data, factorlist){
#	factornames <- intersect(colnames(dat), names(factorlist))
#	if (length(factornames)>0){
#		for(i in factornames){
#				   data[[i]] <- factor(data[[i]], 
#									 levels=as.vector(factorlist[[i]]))
#		}
#	}
#	return(data)
#}

#' @title theme_taxbar
#' @param axis.text.x element_text, x axis tick labels.
#' @param legend.position character, default is "bottom".
#' @param legend.box character, arrangement of legends, default is "horizontal".
#' @param legend.text element_text, legend labels text.
#' @param legend.title element_text, legend title text
#' @param strip.text.x element_text, strip text of x
#' @param strip.background element_rect, the background of x
#' @param ... additional parameters
#' @return updated ggplot object with new theme
#' @importFrom ggplot2 theme element_blank element_text unit element_rect theme_bw
#' @seealso \code{\link[ggplot2]{theme}}
#' @export
#' @examples
#' library(ggplot2)
#' data(test_otu_data)
#' otubar <- ggbartax(test_otu_data, settheme=FALSE) + 
#'     xlab(NULL) + ylab("relative abundance(%)") + 
#'     theme_taxbar()
theme_taxbar <- function(axis.text.x=element_text(angle = -45, hjust = 0, size=12),
						 #panel.grid = element_blank(),
						 legend.position = "bottom",
						 legend.box = "horizontal",
						 legend.text = element_text(size = 8),
						 legend.title=element_blank(),
						 #panel.spacing = unit(0.2, "mm"),
						 strip.text.x=element_text(size=12, face="bold"),
						 strip.background = element_rect(colour="white", fill="grey"),
						 ...
						 ){
    theme_bw()+
    theme(axis.text.x = axis.text.x,
	  axis.text.y = element_text(size=12),
          panel.grid = element_blank(),
          legend.position = legend.position, 
          legend.box = legend.box, 
          legend.text = legend.text, 
          legend.title = legend.title,
          plot.margin = unit(c(0.2,1,0.2,0.2),"cm"),
          panel.spacing = unit(0.2, "mm"),
          strip.text.x = strip.text.x,
          strip.background = strip.background,
		  ...)
}

###' @title legend guides for ggbartax
###' @description
###' the legend guides for `ggbartax`
###' @param keywidth numeric, the keywidth of `guide_legend`.
###' @param keyheight numeric, the keyheight of `guide_legend`.
###' @param ncol integer, the ncol of legend.
###' @param ..., additional parameter.
###' @return the guides of legend.
###' @author Shuangbin Xu
###' @importFrom ggplot2 guides guide_legend
###' @export
##
#taxbarguildes <- function(keywidth=0.3, keyheight=0.3, ncol=5, ...){
#	guides(fill=guide_legend(keywidth = keywidth, 
#							 keyheight = keyheight,
#							 ncol=ncol),...)
#}


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
#' @param size the column name for mapping the size of points, 
#' default is 'pvalue'. 
#' @param removeUnknown logical, whether do not show unknown taxonomy, 
#' default is TRUE.
#' @param layout character, the layout of ggtree, but only "rectangular",
#' "roundrect", "ellipse", "radial", "slanted", "inward_circular" and 
#' "circular" in here, default is "radial".
#' @param linewd numeric, the size of segment of ggtree, default is 0.6.
#' @param bg.tree.color character, the line color of tree, default is '#bed0d1'.
#' @param bg.point.color character, the color of margin of background node points 
#' of tree, default is '#bed0d1'.
#' @param bg.point.stroke numeric, the margin thickness of point of background nodes 
#' of tree, default is 0.2 . 
#' @param bg.point.fill character, the point fill (since point shape is 21) of background 
#' nodes of tree, default is 'white'.
#' @param skpointsize numeric, the point size of skeleton of tree, 
#' default is 2.
#' @param hilight.size numeric, the margin thickness of high light clade, default is 0.2.
#' @param alpha numeric, the alpha of clade, default is 0.4.
#' @param taxlevel positive integer, the full text of clade, default is 5.
#' @param cladetext numeric, the size of text of clade, default is 2.
#' @param tip.annot logcial whether to replace the differential tip labels with shorthand,
#' default is TRUE.
#' @param as.tiplab logical, whether to display the differential tip labels with 'geom_tiplab'
#' of 'ggtree', default is TRUE, if it is FALSE, it will use 'geom_text_repel' of 'ggrepel'.
#' @param factorLevels list, the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this.
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
#' \dontrun{
#' data(kostic2012crc)
#' kostic2012crc %<>% as.phyloseq()
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
#' library(ggplot2)
#' diffcladeplot <- ggdiffclade(diffres,alpha=0.3, linewd=0.2, 
#'                         skpointsize=0.4, 
#'                         taxlevel=5) +
#'                  scale_fill_diff_cladogram(
#'                         values=c('#00AED7', 
#'                                  '#FD9347'
#'                                  )
#'                  ) +
#'                  scale_size_continuous(range = c(1, 3))
#' }
ggdiffclade <- function(obj,...){
    UseMethod("ggdiffclade")
}

#' @method ggdiffclade data.frame
#' @rdname ggdiffclade
#' @importFrom ggplot2 geom_point geom_rect geom_text scale_size_continuous scale_fill_manual 
#' @importFrom magrittr %<>%
#' @importFrom stats as.formula
#' @importFrom rlang sym
#' @export
ggdiffclade.data.frame <- function(obj, nodedf, factorName, size, layout="radial", linewd=0.6, 
                                   bg.tree.color = '#bed0d1', bg.point.color = '#bed0d1', 
                                   bg.point.stroke = 0.2, bg.point.fill = 'white',
                                   skpointsize=2, hilight.size = .2, alpha=0.4, taxlevel=5, 
                                   cladetext=2.5, tip.annot = TRUE, as.tiplab = TRUE, 
                                   factorLevels=NULL, xlim=12, removeUnknown = FALSE, 
                                   reduce=FALSE, type="species", ...){
    params <- list(...)
    if (!is.null(params$size)){
        message("The `size` has been deprecated, Please use `linewd` instead!")
        linewd <- params$size
    }
    treedata <- convert_to_treedata(obj, type=type)
    layout %<>% match.arg(c("rectangular", "roundrect", "ellipse", "circular", "slanted", "radial", "inward_circular"))
    if (!is.null(factorLevels)){nodedf <- setfactorlevels(nodedf, factorLevels)}
    df <- ggtree::fortify(treedata, layout=layout)
    nodedf <- get_node(treedata=df, nodedf=nodedf)
    if (reduce){
        df <- df[!grepl("__un_", df$label),]
        nodedf <- nodedf[!grepl("__un_", as.vector(nodedf[,1])),,drop=FALSE]
    }

    df <- df %>% dplyr::left_join(nodedf[,-1,drop=FALSE], by='node')

    annot.index <- sort(unique(df$nodeDepth[df$nodeDepth > taxlevel]))
    if (!tip.annot){
        annot.index <- annot.index[-length(annot.index)]
    }

    df <- .generate_annot_df(df, annot.index=annot.index, .group=sym(factorName), removeUnknown) 

    offset.max <- max(dplyr::pull(df, .data$nodeDepth), na.rm = TRUE)
    
   if (!missing(size)){
      if (any(grepl(size, c("pvalue", "fdr"), ignore.case = TRUE))){
         tmp <- paste0("-log10(", size, ")")
         df[[tmp]] <- -log10(df[[size]])
         size <- tmp
      }
      pointmapping <- aes_string(fill=factorName, size=size)
   }else{
      pointmapping <- aes_(fill=as.formula(paste("~",factorName)))
   }

   p <- ggtree(df, layout = layout, size = linewd, color = bg.tree.color) +
        geom_point(
           data = ifelse(removeUnknown, td_filter(is.na(!!sym(factorName)) | grepl('__un_', .data$label)),
                    td_filter(is.na(!!sym(factorName)))
                  ),
           mapping = aes_(x = ~x, y = ~y),
           color = bg.point.color,
           fill = bg.point.fill,
           size = skpointsize,
           shape = 21,
           stroke = bg.point.stroke,
        ) +
        ggtree::geom_hilight(
           data = td_mutate(
                    extend = (offset.max + 1 - .data$nodeDepth) * 1,
                    .f = ifelse(removeUnknown, td_filter(!is.na(!!sym(factorName)) & !grepl('__un_', .data$label)),
                                td_filter(!is.na(!!sym(factorName))))
                  ),
           mapping = aes_string(node = "node", fill = factorName, extend = "extend"),
           size = hilight.size,
           alpha = alpha,
           show.legend = FALSE,
        )

   if (removeUnknown){
       flag.clade <- td_filter(!is.na(!!sym(factorName)) & !.data$isTip & !grepl('__un_', .data$label))(p$data)
   }else{
       flag.clade <- td_filter(!is.na(!!sym(factorName)) & !.data$isTip)(p$data)
   }

   if (nrow(flag.clade) > 0){
       p <- p +
        ggtree::geom_cladelab(
           data = td_mutate(
                    offset = (offset.max + 1 - .data$nodeDepth) * 0.88,
                    .f = ifelse(reduce, td_filter(!is.na(!!sym(factorName)) & !.data$isTip & !grepl('__un_', .data$label)),
                                td_filter(!is.na(!!sym(factorName)) & !.data$isTip))
                  ),
           mapping = aes_string(node = "node", label = "label", offset = "offset"),
           geom = 'shadowtext',
           angle = "auto",
           horizontal = FALSE,
           hjust = 0.5,
           fontsize = cladetext,
           barsize = 0,
           barcolor = NA,
           bg.colour = 'white'
        )
   }
   p <- p + 
       ggnewscale::new_scale_color() +
       geom_point(
           data = td_mutate(
                    Label = paste0(letters[seq_len(length(.data$label))], ":", .data$label),
                    .f = td_filter(!is.na(!!sym(factorName)) & .data$nodeDepth %in% annot.index),
                  ),
           mapping = aes_(x = 0, y = 0, color = ~.LABEL),
           size = 0, 
           stroke = 0
       ) + 
       geom_point( 
           data = ifelse(removeUnknown, td_filter(!grepl("__un_", .data$label) & !is.na(!!sym(factorName))), td_filter(!is.na(!!sym(factorName)))),
           mapping = pointmapping, 
           shape = 21, 
           stroke = .1,
       )
   if (!as.tiplab){
       p <- p +
            ggrepel::geom_text_repel(
               data = ifelse(removeUnknown, td_filter(!is.na(!!sym(factorName)) & .data$isTip & !grepl("__un_", .data$label)),
                             td_filter(!is.na(!!sym(factorName)) & .data$isTip)),
               mapping = aes_(x = ~x, y = ~y, label = ~label),
               size = cladetext,
               min.segment.length = 0, 
               bg.colour = 'white'
            )
   }else{
       p <- p +
            ggtree::geom_tiplab(
               data = ifelse(removeUnknown, td_filter(!is.na(!!sym(factorName)) & !grepl("__un_", .data$label)),
                        td_filter(!is.na(!!sym(factorName)))),
               offset = 0.78,
               size = cladetext,
               geom = 'shadowtext',
               color = 'black',
               bg.colour = 'white'
            ) 
   }
       
   p <- p + 
        theme(
          legend.key.width = unit(.3, 'cm'),
          legend.key.height = unit(.3, 'cm'),
          legend.text = element_text(size=6),
          legend.title = element_text(size=8),
          legend.margin = ggplot2::margin(-.25, 0, 0, 0, 'cm') 
        )
   p <- p + scale_fill_diff_cladogram()
   return (p)
}


#' @method ggdiffclade diffAnalysisClass
#' @rdname ggdiffclade
#' @importFrom magrittr %>%
#' @export
ggdiffclade.diffAnalysisClass <- function(obj, size, removeUnknown=TRUE, ...){
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
    if (missing(size)){
        size <- 'pvalue'
    }
    if (removeUnknown){
        tmpflag <- grep("__un_",as.vector(nodedfres$f))
        if (length(tmpflag)>0){
   	        nodedfres <- nodedfres[-tmpflag,,drop=FALSE]
        }
    }
    p <- ggdiffclade(obj=taxda,
    			     nodedf=nodedfres,
    			     factorName=classname,
                     size = size,
                     type=extract_args(obj, "type"),
    			     ...)
    return(p)
}
#' @importFrom dplyr rename
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

## #' @importFrom ggplot2 guides guide_legend
## #' @keywords internal
## guides_diffclade <- function(...){
##     guides(fill= guide_legend(keywidth = 0.5, keyheight = 0.5, order=1),
##     	   size=guide_legend(keywidth = 0.5, keyheight = 0.5, order=2),
##     	   color = guide_legend(keywidth = 0.1, ncol=1, keyheight = 0.6, order = 3),...)
## }

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

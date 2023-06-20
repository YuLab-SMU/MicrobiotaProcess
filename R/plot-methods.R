#' plotting the abundance of taxa via specified taxonomy class
#' @rdname mp_plot_abundance-methods
#' @param .data MPSE object or tbl_mpse object
#' @param .abundance the column name of abundance to be plotted.
#' @param .group the column name of group to be calculated and plotted,
#' default is NULL.
#' @param taxa.class name of taxonomy class, default is NULL, meaning the
#' Phylum class will be plotted.
#' @param topn integer the number of the top most abundant, default is 10.
#' @param relative logical whether calculate the relative abundance and plotted.
#' @param force logical whether calculate the relative abundance forcibly
#'  when the abundance is not be rarefied, default is FALSE.
#' @param plot.group logical whether plotting the abundance of specified taxa.class 
#' taxonomy with group not sample level, default is FALSE.
#' @param geom character which type plot, options is 'flowbar' 'bar' and 'heatmap', default is
#' 'flowbar'.
#' @param feature.dist character the method to calculate the distance between the features,
#' based on the '.abundance' of 'taxa.class', default is 'bray', options refer to
#' the 'distmethod' of [mp_cal_dist()] (except unifrac related).
#' @param feature.hclust character the agglomeration method for the features, default is
#' 'average', options are 'single', 'complete', 'average', 'ward.D', 'ward.D2', 'centroid'
#' 'median' and 'mcquitty'.
#' @param sample.dist character the method to calculate the distance between the samples
#' based on the '.abundance' of 'taxa.class', default is 'bray', options refer to the 
#' 'distmethod' of [mp_cal_dist()] (except unifrac related).
#' @param sample.hclust character the agglomeration method for the samples, default is
#' 'average', options are 'single', 'complete', 'average', 'ward.D', 'ward.D2', 'centroid'
#' 'median' and 'mcquitty'.
#' @param .sec.group the column name of second group to be plotted with nested facet,
#' default is NULL, this argument will be deprecated in the next version.
#' @param ... additional parameters, when the geom = "flowbar", it can specify the parameters of 
#' 'geom_stratum' of 'ggalluvial', when the geom = 'bar', it can specify the parameters of 
#' 'geom_bar' of 'ggplot2', when the geom = "heatmap", it can specify the parameter of
#' 'geom_tile' of 'ggplot2'.
#' @param rmun logical whether to remove the unknown taxa, such as "g__un_xxx",
#' default is FALSE (the unknown taxa class will be considered as 'Others').
#' @param rm.zero logical whether to display the zero abundance, which only work with geom='heatmap'
#' default is TRUE.
#' @author Shuangbin Xu
#' @export
#' @examples
#' \dontrun{
#' data(mouse.time.mpse)
#' mouse.time.mpse %<>%
#'   mp_rrarefy()
#' mouse.time.mpse
#' mouse.time.mpse %<>%
#'   mp_cal_abundance(.abundance=RareAbundance, action="add") %>%
#'   mp_cal_abundance(.abundance=RareAbundance, .group=time, action="add")
#' mouse.time.mpse
#' p1 <- mouse.time.mpse %>%
#'       mp_plot_abundance(.abundance=RelRareAbundanceBySample, 
#'                         .group=time, 
#'                         taxa.class="Phylum", 
#'                         topn=20)
#' p2 <- mouse.time.mpse %>%
#'       mp_plot_abundance(.abundance = Abundance,
#'                         taxa.class = Phylum,
#'                         topn = 20,
#'                         relative = FALSE,
#'                         force = TRUE
#'                        )
#' p3 <- mouse.time.mpse %>%
#'       mp_plot_abundance(.abundance = RareAbundance, 
#'                         .group = time,
#'                         taxa.class = Phylum, 
#'                         topn = 20,
#'                         relative = FALSE,
#'                         force = TRUE
#'                         )
#' p4 <- mouse.time.mpse %>%
#'       mp_plot_abundance(.abundance = RareAbundance,
#'                         .group = time,
#'                         taxa.class = Phylum,
#'                         topn = 20,
#'                         relative = FALSE,
#'                         force = TRUE,
#'                         plot.group = TRUE
#'                         )
#' }
setGeneric("mp_plot_abundance", 
           function(
              .data, 
              .abundance = NULL, 
              .group = NULL, 
              taxa.class = NULL, 
              topn = 10, 
              relative = TRUE,
              force = FALSE,
              plot.group = FALSE,
              geom = "flowbar",
              feature.dist = "bray",
              feature.hclust = "average",
              sample.dist = "bray",
              sample.hclust = "average",
              .sec.group = NULL,
              rmun = FALSE,
              rm.zero = TRUE,
              ...
           )
           standardGeneric("mp_plot_abundance")
)

.internal_plot_abundance <- function(.data, 
                                .abundance, 
                                .group = NULL, 
                                taxa.class = NULL, 
                                topn = 10,
                                relative = TRUE,
                                force = FALSE,
                                plot.group = FALSE,
                                geom = "flowbar",
                                feature.dist = "bray",
                                feature.hclust = "average",
                                sample.dist = "bray",
                                sample.hclust = "average",
                                .sec.group = NULL,
                                rmun = FALSE,
                                rm.zero = TRUE,
                                ...
                                ){
    .abundance <- rlang::enquo(.abundance)
    .group <- rlang::enquo(.group)
    .sec.group <- rlang::enquo(.sec.group)
    taxa.class <- rlang::enquo(taxa.class)
    geom %<>% match.arg(c("flowbar", "bar", "heatmap"))
    if (geom=="heatmap"){
        plot.group = FALSE
    }

    if (rlang::quo_is_null(taxa.class) || 
        (.data %>% mp_extract_tree() %>% is.null() %>% suppressMessages())){
        taxa.class <- rlang::sym("OTU")
    }

    if (rlang::quo_is_missing(.abundance)){
        .abundance <- rlang::sym("Abundance")
    }
    
    if (relative){
         if (force){
             abundance.nm <- paste0("Rel", rlang::as_name(.abundance), 'BySample')
         }else{
             abundance.nm <- "RelRareAbundanceBySample"
         }
         ylabs <- "Relative Abundance (%)"
     }else{
         if (force){
             abundance.nm <- rlang::as_name(.abundance)
         }else{
             abundance.nm <- "RareAbundance"
         }

         ylabs <- abundance.nm
     }

     if (!rlang::quo_is_null(.group) && plot.group){
         gp <- quo.vect_to_str.vect(.group)
         prefixBy <- paste0("By", paste0(gp, collapse="And"))
         if (!force && grepl('BySample$', abundance.nm)){
             abundance.nm %<>% gsub('BySample', "", .)
         }
         axis.x <- rlang::as_name(gp[1])
     }else{
         prefixBy <- ""
         axis.x <- "Sample"
     }
     abundance.nm <- paste0(abundance.nm, prefixBy)
     
     if (!(!relative && grepl('^Rel', abundance.nm) && !plot.group)){
         AbundBy <- abundance.nm %>% gsub("^Rel", "", .)
     }else{
         AbundBy <- abundance.nm
     }

     internal.cal <- FALSE
     if (!any(grepl(paste0("^", AbundBy), .data %>% mp_extract_feature() %>% colnames()))){
         internal.cal <- TRUE
         if (!rlang::quo_is_null(.group) && plot.group){
             .data %<>% mp_cal_abundance(.abundance=!!.abundance, .group=!!.group, force=force, relative=relative)
             gp <- quo.vect_to_str.vect(.group)
             prefixBy <- paste0("By", paste0(gp, collapse = "And"))
             if (relative){
                 if (force){
                     AbundBy <- paste0(rlang::as_name(.abundance), prefixBy)
                 }else{
                     AbundBy <- paste0("RareAbundance", prefixBy)
                 }
                 abundance.nm <- paste0('Rel', AbundBy)
             }else{
                 if (force){
                     AbundBy <- paste0(rlang::as_name(.abundance), prefixBy)
                 }else{
                     AbundBy <- paste0("RareAbundance", prefixBy)
                 }
                 abundance.nm <- AbundBy
             }
         }else{
             .data %<>% mp_cal_abundance(.abundance=!!.abundance, force=force, relative=relative)
             if (relative){
                 if (force){
                     AbundBy <- paste0(rlang::as_name(.abundance), 'BySample')
                 }else{
                     AbundBy <- "RareAbundanceBySample"
                 }
                 abundance.nm <- paste0('Rel', AbundBy)
             }else{
                 if (force){
                     AbundBy <- paste0(rlang::as_name(.abundance), 'BySample')
                     abundance.nm <- rlang::as_name(.abundance)
                     #if (rlang::as_name(.abundance) %in% c('Abundance', 'RareAbundance')){
                     #    abundance.nm <- 
                     #}
                 }else{
                     AbundBy <- "RareAbundanceBySample"
                     abundance.nm <- "RareAbundance"
                 }
             }
         }
     }
     
     if (!plot.group && !internal.cal){
         AbundBy %<>% gsub('BySample', "", .) %>% paste0('BySample')
         flag_prefix1 <- grepl('BySample$', abundance.nm)
         flag_prefix2 <- grepl("^Rel", abundance.nm)
         abundance.nm %<>% gsub('BySample', "", .) %>% gsub('Rel', "", .)
         if (flag_prefix1){
             abundance.nm %<>% paste0('BySample')
         }
         if (flag_prefix2){
             abundance.nm <- paste0('Rel', abundance.nm)
         }
     }

     
     tbl <- .data %>% 
            mp_extract_abundance(taxa.class=!!taxa.class, topn = topn, rmun = rmun) %>%
            tidyr::unnest(cols=AbundBy) %>% 
            dplyr::mutate(label=forcats::fct_rev(.data$label)) %>%
            dplyr::rename(!!taxa.class:="label") %>% 
            suppressMessages()

     if(geom %in% c("bar", 'flowbar')){
         if (geom == "flowbar"){
            p <- ggplot(data = tbl,
                        mapping = aes_string(
                           x = axis.x,
                           y = abundance.nm,
                           alluvium = rlang::as_name(taxa.class),
                           fill = rlang::as_name(taxa.class))
                 ) +
                 ggalluvial::geom_flow(stat="alluvium", lode.guidance = "frontback", color = "darkgray", ...) +
                 ggalluvial::geom_stratum(stat="alluvium", ...) 
         }else if (geom == "bar"){
            p <- ggplot(data = tbl,
                        mapping = aes_string(
                           x = axis.x,
                           y = abundance.nm,
                           fill = rlang::as_name(taxa.class)
                        )
                 ) +
                 geom_bar(stat = 'identity', ...)
         }
         p <- p + 
              ggplot2::labs(x=NULL, y=ylabs) +
              scale_fill_manual(
                  values = rev(get_cols(tbl %>% pull(!!taxa.class) %>% unique() %>% length())),
                  guide = guide_legend(reverse=TRUE)
              ) 
          
         if (!rlang::quo_is_null(.group)){
             gp <- quo.vect_to_str.vect(.group)
             if (!rlang::quo_is_null(.sec.group)){
                 warning_wrap("The .sec.group will be depcrecated, please use .group argument, which
                              supports multiple groups.")
                 gp2 <- quo.vect_to_str.vect(.sec.group)
                 if (!gp2 %in% gp){
                     gp <- append(gp, gp2, after=1)
                 }
             }
             if (plot.group && length(gp)>1 ){
                 gpformula <- as.formula(paste0(". ~ ", paste0(gp[-1], collapse="+")))
             }else if (!plot.group){
                 gpformula <- as.formula(paste0(". ~ ", paste0(gp, collapse="+")))
             }else{
                 gpformula <- NULL
             }
         }else{
             gpformula <- NULL
         }
         if (!is.null(gpformula)){
             p <- p + ggh4x::facet_nested(gpformula, scales="free_x", space="free")
         }   
         gb <- ggplot2::ggplot_build(p)
         expand.value = max(gb$layout$panel_scales_y[[1]]$range$range) * 0.5 / 100                  
         p <- p + 
              scale_y_continuous(expand = c(0, 0, 0, expand.value)) +
              theme_taxbar()

     }else if(geom=="heatmap"){
         lab.sty <- list(xlab(NULL),ylab(NULL))
         tbl %<>% dplyr::mutate_if(is.factor, as.character)
         leg.tl <- abundance.nm %>% gsub("BySample", "", .)
         p <- ggplot(data=tbl,
                     mapping = aes_string(
                         x = axis.x,
                         y = rlang::as_name(taxa.class),
                         fill = abundance.nm
                     )
               ) 
         if (rm.zero) {
             p <- p + ggplot2::geom_tile(
                   data = td_filter(!!rlang::sym(abundance.nm)!=0),
                   ...
               )
         }else{
             p <- p + ggplot2::geom_tile(...)
         }
         p <- p +
               theme(axis.text.x=element_text(angle=-45, hjust=0), 
                     panel.background=element_blank(), 
                     panel.border=element_rect(size=1, fill=NA)) + 
               ggplot2::scale_y_discrete(position="right") +
               ggplot2::guides(fill=ggplot2::guide_colourbar(title=leg.tl)) +
               lab.sty

         tbl2 <- tbl %>% 
                 as.tbl_mpse(.OTU = !!taxa.class, 
                             .Sample = !!rlang::sym(axis.x), 
                             .Abundance = !!rlang::sym(abundance.nm)
                 )

         sample.dist <- mp_cal_dist(tbl2, 
                                    .abundance = !!rlang::sym("Abundance"), 
                                    distmethod = sample.dist, 
                                    action = "get"
                        )

         feature.dist <- mp_cal_dist(tbl2, 
                                     .abundance = !!rlang::sym("Abundance"),
                                     distmethod = feature.dist, 
                                     action = "get", 
                                     cal.feature.dist = TRUE
                         )

         indexname <- leg.tl

         sample.hc <- hclust(sample.dist, method = sample.hclust) %>% ape::as.phylo()
         feature.hclust <- hclust(feature.dist, method = feature.hclust) %>% ape::as.phylo()

         if (!rlang::quo_is_null(.group)){
             gp <- quo.vect_to_str.vect(.group)
             if (!rlang::quo_is_null(.sec.group)){
                 warning_wrap("The .sec.group will be depcrecated, please use .group argument, which
                              supports multiple groups.")
                 gp2 <- quo.vect_to_str.vect(.sec.group)
                 gp <- append(gp, gp2, after=1)
             }
             gp %<>% unique()
             for (i in seq_len(length(gp))){
                 sampleda <- .data %>% mp_extract_sample()
                 f <- ggplot() +
                      ggplot2::geom_tile(
                         data = sampleda,
                         mapping = aes(x=!!rlang::sym("Sample"), y=gp[i], fill=!!rlang::sym(gp[i]))
                      ) +
                      ggplot2::scale_y_discrete(position="right", expand = c(0, 0), labels=gp[i]) +
                      theme(axis.text.x=element_blank(), 
                            axis.ticks.x=element_blank(), 
                            panel.background=element_blank()
                      ) +
                      lab.sty 
                 p %<>% insert_top(f, height = 0.04 + 0.0025 * (i - 1))
             }
             indexname <- c(indexname, gp)
         }
         p2 <- ggtree(feature.hclust, branch.length = "none", size = 0.8)
         p3 <- ggtree(sample.hc, branch.length = "none", size = 0.8, layout = "dendrogram")
         p %<>% insert_left(p2, width = 0.1)
         if (!rlang::quo_is_null(.group)){
             p %<>% insert_top(p3, height = 0.1 + 0.01 * length(gp))
         }else{
             p %<>% insert_top(p3, height = 0.1)
         }
         p$index <- c(indexname, paste0("tree", seq_len(2)))
         p %<>% add_class("aplot.heatmap")
     }
     
     return (p)
}

#' @rdname mp_plot_abundance-methods
#' @aliases mp_plot_abundance,MPSE
#' @export mp_plot_abundance
setMethod("mp_plot_abundance", signature(.data="MPSE"), .internal_plot_abundance)

#' @rdname mp_plot_abundance-methods
#' @aliases mp_plot_abundance,tbl_mpse
#' @export mp_plot_abundance
setMethod("mp_plot_abundance", signature(.data="tbl_mpse"), .internal_plot_abundance)

#' @rdname mp_plot_abundance-methods
#' @aliases mp_plot_abundance,grouped_df_mpse
#' @export mp_plot_abundance
setMethod("mp_plot_abundance", signature(.data="grouped_df_mpse"), .internal_plot_abundance)

#' Plotting the alpha diversity between samples or groups.
#' @rdname mp_plot_alpha-methods
#' @param .data MPSE or tbl_mpse object
#' @param .group the column name of sample group information
#' @param .alpha the column name of alpha index after run mp_cal_alpha or mp_cal_pd_metric.
#' @param test the name of the statistical test, default is 'wilcox.test'
#' @param comparisons A list of length-2 vectors. The entries in the vector are
#' either the names of 2 values on the x-axis or the 2 integers that 
#' correspond to the index of the columns of interest, default is NULL, meaning 
#' it will be calculated automatically with the names in the .group.
#' @param step_increase numeric vector with the increase in fraction of total
#' height for every additional comparison to minimize overlap, default is 0.05.
#' @param ... additional parameters, see also \code{\link[ggsignif]{geom_signif}}
#' @author Shuangbin Xu
#' @export
#' @examples
#' \dontrun{
#' data(mouse.time.mpse)
#' mpse <- mouse.time.mpse %>%
#'         mp_rrarefy() %>%
#'         mp_cal_alpha(.abundance=RareAbundance)
#' mpse
#' p <- mpse %>% 
#'      mp_plot_alpha(.group=time, .alpha=c(Observe, Shannon, Pielou))
#' p
#' }
setGeneric("mp_plot_alpha", 
  function(.data,
           .group,
           .alpha=c("Observe", "Shannon"),
           test = "wilcox.test",
           comparisons = NULL,
           step_increase = 0.05,
           ...
  )
  standardGeneric("mp_plot_alpha")
)

.internal_plot_alpha <- function(
    .data, 
    .group, 
    .alpha = c("Observe", "Shannon"), 
    test = "wilcox.test", 
    comparisons = NULL, 
    step_increase = 0.05, 
    ...){
    
    .group <- rlang::enquo(.group)
    .alpha <- rlang::enquo(.alpha)

    if (rlang::quo_is_missing(.group)){
        .group <- NULL
        #rlang::abort("The .group column name is required for the visualization of alpha diversity")
    }

    newlevels <- quo.vect_to_str.vect(.alpha)

    if ("J" %in% newlevels){
        warning_wrap("The 'J' exists in the .alpha parameter, it has been deprecated and it will not be supported 
                     in the next released version, please use 'Pielou' to replace it.", call. = FALSE)
        newlevels[newlevels == "J"] <- "Pielou"
    }

    tbl <- .data %>% 
        mp_extract_sample() %>%
        dplyr::select(!!rlang::sym("Sample"), !!newlevels, !!.group) %>%
        tidyr::pivot_longer(cols=-c(rlang::sym("Sample"), !!.group), 
                            names_to="Measure", values_to="Alpha")

    tbl$Measure <- factor(tbl$Measure, levels=newlevels)

    if (!is.null(.group)){
        gp <- quo.vect_to_str.vect(.group)
        if (is.null(comparisons)){
           comparisons <- tbl %>% 
                          pull(!!rlang::sym(gp[1])) %>% 
                          unique() %>% 
                          utils::combn(2) %>% 
                          apply(2, list) %>% 
                          unlist(recursive = FALSE)
        }
        mapping <- aes_string(x = gp[1], 
                              y = "Alpha", 
                              fill = gp[1])
    }else{
        gp <- NULL
        mapping <- aes_string(x = "Sample", 
                              y = "Alpha"
                   )
    }

    p <- ggplot(data=tbl, mapping = mapping)

    if (!is.null(gp)){
        if (is.numeric(tbl[[gp[1]]])){
           if (length(gp) > 1 ){
               mapping <- modifyList(mapping, aes_string(color=gp[2]))
               gp <- gp[-c(1, 2)]
           }else{
               mapping <- modifyList(mapping, aes_string(fill=NULL))
               gp <- gp[-1]
           }
           p <- ggplot(data = tbl, mapping = mapping)
           smoothparam <- modifyList(list(method="lm", se=TRUE), list(...))
           p <- p + 
                geom_point()
           p <- p +
                do.call(ggplot2::geom_smooth, smoothparam) 
        }else{
           p <- p +
              gghalves::geom_half_violin(color=NA, side="l", trim=FALSE) +
              gghalves::geom_half_point(side="r", shape=21, alpha=0.8) +
              ggplot2::geom_boxplot(aes_string(color=gp[1]),
                                    fill = NA,
                                    position=ggplot2::position_nudge(x=.22),
                                    size = 0.6,
                                    width = 0.2,
                                    outlier.shape = NA
                                    ) +
              ggsignif::geom_signif(comparisons=comparisons, test=test, step_increase=step_increase, ...) +
              ggplot2::scale_fill_manual(values=get_cols(tbl %>% pull(!!rlang::sym(gp[1])) %>% unique() %>% length())) +
              ggplot2::scale_color_manual(values=get_cols(tbl %>% pull(!!rlang::sym(gp[1])) %>% unique() %>% length()))
           gp <- gp[-1]
        }
        if (length(gp) > 0){
           gpformula <- as.formula(paste0("Measure ~ ", paste0(gp, collapse="+")))
           p <- p + ggh4x::facet_nested(gpformula, scales="free", space="free_x")
        }else{
           p <- p + ggplot2::facet_wrap(facets=ggplot2::vars(!!rlang::sym("Measure")), scales="free_y", nrow=1)
        }
    }else{
        p <- p +
             ggplot2::geom_col() +
             ggplot2::facet_wrap(facets=ggplot2::vars(!!rlang::sym("Measure")), scales="free_y")
    }

    p <- p + 
         labs(x=NULL, y="Alpha Index Value") + 
         theme_taxbar(
           legend.position="right", 
           strip.text.y = ggplot2::element_text(size = 12, face = "bold")
         )
    return (p)
}


#' @rdname mp_plot_alpha-methods
#' @aliases mp_plot_alpha,MPSE
#' @export mp_plot_alpha
setMethod("mp_plot_alpha", signature(.data="MPSE"), .internal_plot_alpha)

#' @rdname mp_plot_alpha-methods
#' @aliases mp_plot_alpha,tbl_mpse
#' @export mp_plot_alpha
setMethod("mp_plot_alpha", signature(.data="tbl_mpse"), .internal_plot_alpha)

#' @rdname mp_plot_alpha-methods
#' @aliases mp_plot_alpha,grouped_df_mpse
#' @export mp_plot_alpha
setMethod("mp_plot_alpha", signature(.data="grouped_df_mpse"), .internal_plot_alpha)

#' Plotting the different number of OTU between groups with Venn Diagram.
#' @rdname mp_plot_venn-methods
#' @param .data MPSE object or tbl_mpse object
#' @param .group the column names of group to be visualized
#' @param .venn the column names of result after run \code{mp_cal_venn}.
#' @param ... additional parameters, such as 'size', 'label_size', 'edge_size' etc, 
#' see also 'ggVennDiagram'.
#' @author Shuangbin Xu
#' @examples
#' \dontrun{
#' data(mouse.time.mpse)
#' mpse <- mouse.time.mpse %>%
#'         mp_rrarefy() %>%
#'         mp_cal_venn(.abundance=RareAbundance, .group=time, action="add")
#' mpse
#' p <- mpse %>% mp_plot_venn(.group=time, .venn=vennOftime) 
#' p
#' }
setGeneric("mp_plot_venn", function(.data, .group, .venn=NULL, ...) standardGeneric("mp_plot_venn"))

.internal_plot_venn <- function(.data, .group, .venn=NULL, ...){
    .group <- rlang::enquo(.group)
    .venn <- rlang::enquo(.venn)
    if (rlang::quo_is_null(.venn)){
        .venn <- rlang::sym(paste0("vennOf", rlang::as_name(.group)))
    }
    p <- .data %>% 
        mp_extract_sample() %>% 
        dplyr::select(!!.group, !!.venn) %>% 
        dplyr::distinct() %>%
        pull(var=!!.venn, name=!!.group) %>%
        ggVennDiagram::ggVennDiagram(., ...)
    return(p)
}

#' @rdname mp_plot_venn-methods
#' @aliases mp_plot_venn,MPSE
#' @export mp_plot_venn
setMethod("mp_plot_venn", signature(.data="MPSE"), .internal_plot_venn)

#' @rdname mp_plot_venn-methods
#' @aliases mp_plot_venn,tbl_mpse
#' @export mp_plot_venn 
setMethod("mp_plot_venn", signature(.data="tbl_mpse"), .internal_plot_venn)

#' @rdname mp_plot_venn-methods
#' @aliases mp_plot_venn,grouped_df_mpse
#' @export mp_plot_venn
setMethod("mp_plot_venn", signature(.data="grouped_df_mpse"), .internal_plot_venn)

#' Plotting the different number of OTU between group via UpSet plot
#' @rdname mp_plot_upset-methods
#' @param .data MPSE obejct or tbl_mpse object
#' @param .group the column name of group
#' @param .upset the column name of result after run \code{mp_cal_upset}
#' @param ... additional parameters, see also 'scale_x_upset' of 'ggupset'.
#' @export
#' @author Shuangbin Xu
#' @examples
#' \dontrun{
#' data(mouse.time.mpse)
#' mpse <- mouse.time.mpse %>%
#'         mp_rrarefy(.abundance=Abundance) %>%
#'         mp_cal_upset(.abundance=RareAbundance, .group=time) 
#' mpse
#' p <- mpse %>% mp_plot_upset(.group=time, .upset=ggupsetOftime)
#' p
#' }
setGeneric("mp_plot_upset", function(.data, .group, .upset=NULL, ...) standardGeneric("mp_plot_upset"))

.internal_plot_upset <- function(.data, .group, .upset=NULL, ...){
    .group <- rlang::enquo(.group)
    .upset <- rlang::enquo(.upset)

    if (rlang::quo_is_null(.upset)){
        .upset <- rlang::sym(paste0("ggupsetOf", rlang::as_name(.group)))
    }

    p <- .data %>%
         mp_extract_feature() %>%
         dplyr::select(!!rlang::sym("OTU"), !!.upset) %>%
         dplyr::filter(vapply(!!.upset, length, numeric(1))>0) %>%
         ggplot(mapping=aes(x=!!.upset)) +
         geom_bar() +
         ggupset::scale_x_upset(...) +
         ggupset::theme_combmatrix(combmatrix.label.extra_spacing=30) +
         ggplot2::labs(y=rlang::as_name(.group))

    return (p)
}

#' @rdname mp_plot_upset-methods
#' @aliases mp_plot_upset,MPSE
#' @export mp_plot_upset
setMethod("mp_plot_upset", signature(.data="MPSE"), .internal_plot_upset)

#' @rdname mp_plot_upset-methods
#' @aliases mp_plot_upset,tbl_mpse
#' @export mp_plot_upset
setMethod("mp_plot_upset", signature(.data="tbl_mpse"), .internal_plot_upset)

#' @rdname mp_plot_upset-methods
#' @aliases mp_plot_upset,grouped_df_mpse
#' @export mp_plot_upset
setMethod("mp_plot_upset", signature(.data="grouped_df_mpse"), .internal_plot_upset)

#' Rarefaction alpha index with MPSE
#' @rdname mp_plot_rarecurve-methods
#' @param .data MPSE object or tbl_mpse after it was performed \code{mp_cal_rarecurve} with \code{action='add'}
#' @param .rare the column names of 
#' @param .alpha the names of alpha index, which should be one or more of Observe, 
#' ACE, Chao1, default is Observe.
#' @param .group the column names of group, default is NULL, when it is provided, the
#' rarecurve lines will group and color with the \code{group}.
#' @param nrow integer Number of rows in \code{\link[ggplot2]{facet_wrap}}.
#' @param plot.group logical whether to combine the samples, default is FALSE,
#' when it is TRUE, the samples of same group will be represented by their group.
#' @param ... additional parameters, see also \code{\link[ggplot2]{geom_smooth}}.
#' @author Shuangbin Xu
#' @export
#' @examples
#' \dontrun{
#' data(mouse.time.mpse)
#' mpse <- mouse.time.mpse %>%
#'         mp_rrarefy()
#' mpse
#' mpse %<>% mp_cal_rarecurve(.abundance=RareAbundance, chunks=100, action="add") 
#' mpse
#' p1 <- mpse %>% mp_plot_rarecurve(.rare=RareAbundanceRarecurve, .alpha="Observe")
#' p2 <- mpse %>% mp_plot_rarecurve(.rare=RareAbundanceRarecurve, .alpha="Observe", .group=time)
#' p3 <- mpse %>% mp_plot_rarecurve(.rare=RareAbundanceRarecurve, .alpha="Observe", .group=time, plot.group=TRUE)
#' }
setGeneric("mp_plot_rarecurve", function(.data, 
                                         .rare, 
                                         .alpha=c('Observe', 'Chao1', 'ACE'),
                                         .group = NULL, 
                                         nrow = 1,
                                         plot.group = FALSE,
                                         ...
                                         )
    standardGeneric("mp_plot_rarecurve")
)

.internal_plot_rarecurve <- function(.data, .rare, .alpha = c('Observe', 'Chao1', 'ACE'), .group=NULL, nrow=1, plot.group=FALSE, ...){
    .rare <- rlang::enquo(.rare)
    .alpha <- rlang::enquo(.alpha)
    .group <- rlang::enquo(.group)
    params <- list(...)

    newlevels <- quo.vect_to_str.vect(.alpha)

    tbl <- .data %>% 
           mp_extract_sample() %>%
           dplyr::select(!!rlang::sym("Sample"), !!.rare) %>%
           tidyr::unnest(!!.rare) %>% 
           dplyr::filter(.data$Alpha %in% newlevels)

    tbl$Alpha <- factor(tbl$Alpha, newlevels)

    if (!rlang::quo_is_null(.group)){
        if (plot.group){
            maps <- aes_string(x="readsNums", 
                                  y="value", 
                                  color=rlang::as_name(.group), 
                                  fill=rlang::as_name(.group),
                                  group=rlang::as_name(.group)
                    )
            if (!'se' %in% names(params)){
                params$se = TRUE
            }
        }else{
            maps <- aes_string(x="readsNums", 
                                  y="value", 
                                  color=rlang::as_name(.group), 
                                  fill=rlang::as_name(.group), 
                                  group="Sample")
            if (!'se' %in% names(params)){
                params$se = FALSE
            }
        }
    }else{
        maps <- aes_(x=~readsNums, 
                     y=~value, 
                     color=~Sample, 
                     fill=~Sample)
        if (! 'se' %in% names(params)){
            params$se = FALSE
        }
    }

    if ("alpha" %in% names(params)){
        alpha <- params$alpha
        params$alpha <- NULL
    }else{
        alpha <- 0.5
    }

    p <- ggplot2::ggplot(data=tbl, mapping=maps) 
    smoothlayer <- do.call(geom_smooth, c(method="lm", formula=y~log(x), alpha=alpha, params))
    p <- p + smoothlayer + ggplot2::scale_y_continuous(expand=c(0, 0, 0, 0.8), oob=scales::squish)

    if (!rlang::quo_is_null(.group)){
         p <- p + ggplot2::stat_summary(fun.data = 'mean_se', geom = "ribbon", alpha = alpha)
    }
     
    p <- p + facet_wrap(facets=ggplot2::vars(!!rlang::sym("Alpha")), scales="free_y", nrow=nrow) 

    return (p)
}

#' @rdname mp_plot_rarecurve-methods
#' @aliases mp_plot_rarecurve,MPSE
#' @export mp_plot_rarecurve
setMethod("mp_plot_rarecurve", signature(.data="MPSE"), .internal_plot_rarecurve)

#' @rdname mp_plot_rarecurve-methods
#' @aliases mp_plot_rarecurve,tbl_mpse
#' @export mp_plot_rarecurve
setMethod("mp_plot_rarecurve", signature(.data="tbl_mpse"), .internal_plot_rarecurve)

#' @rdname mp_plot_rarecurve-methods
#' @aliases mp_plot_rarecurve,grouped_tbl_mpse
#' @export mp_plot_rarecurve
setMethod("mp_plot_rarecurve", signature(.data="grouped_df_mpse"), .internal_plot_rarecurve)

#' Plotting the distance between the samples with heatmap or boxplot.
#' @rdname mp_plot_dist-methods
#' @param .data the MPSE or tbl_mpse object after [mp_cal_dist()] is performed with action="add"
#' @param .distmethod the column names of distance of samples, it will generate after
#' [mp_cal_dist()] is performed.
#' @param .group the column names of group, default is NULL, when it is not provided 
#' the heatmap of distance between samples will be returned. If it is provided and
#' \code{group.test} is TURE, the comparisons boxplot of distance between the group
#' will be returned, but when \code{group.test} is FALSE, the heatmap of distance between
#' samples with group information will be returned.
#' @param group.test logical default is FALSE, see the \code{.group} argument.
#' @param hclustmethod character the method of \code{\link[stats]{hclust}}, default is
#' 'average' (= UPGMA).
#' @param test the name of the statistical test, default is 'wilcox.test'
#' @param comparisons A list of length-2 vectors. The entries in the vector are
#' either the names of 2 values on the x-axis or the 2 integers that
#' correspond to the index of the columns of interest, default is NULL, meaning
#' it will be calculated automatically with the names in the .group.
#' @param step_increase numeric vector with the increase in fraction of total
#' height for every additional comparison to minimize overlap, default is 0.1.
#' @param ... additional parameters, see also \code{\link[ggsignif]{geom_signif}}
#' @seealso [mp_cal_dist()] and [mp_extract_dist()]
#' @author Shuangbin Xu
#' @export
#' @examples
#' \dontrun{
#' data(mouse.time.mpse)
#' mouse.time.mpse %<>% mp_decostand(.abundance=Abundance)
#' mouse.time.mpse
#' mouse.time.mpse %<>% 
#'   mp_cal_dist(.abundance=hellinger, distmethod="bray")
#' mouse.time.mpse
#' p1 <- mouse.time.mpse %>% 
#'         mp_plot_dist(.distmethod=bray)
#' p2 <- mouse.time.mpse %>% 
#'         mp_plot_dist(.distmethod=bray, .group=time, group.test=TRUE)
#' p3 <- mouse.time.mpse %>% 
#'         mp_plot_dist(.distmethod=bray, .group=time)
#' }
setGeneric("mp_plot_dist", 
  function(.data, 
           .distmethod, 
           .group = NULL, 
           group.test = FALSE,
           hclustmethod = "average",
           test = "wilcox.test",
           comparisons = NULL,
           step_increase = 0.1,
           ...
  ) 
  standardGeneric("mp_plot_dist")
)

.internal_plot_dist <- function(
  .data, 
  .distmethod,
  .group = NULL,
  group.test = FALSE, 
  hclustmethod = "average",
  test = "wilcox.test",
  comparisons = NULL,
  step_increase = 0.1,
  ...
){
  .distmethod <- rlang::enquo(.distmethod)
  .group <- rlang::enquo(.group)
  params <- list(...)

  if (!rlang::quo_is_null(.group) && group.test){
     tbl <- .data %>%
            mp_extract_dist(distmethod=rlang::as_name(.distmethod), .group=!!.group)  
     if (is.null(comparisons)){
         comparisons <- tbl %>%
                        pull("GroupsComparison") %>%
                        unique() %>%
                        utils::combn(2) %>%
                        apply(2, list) %>%
                        unlist(recursive = FALSE)
     }
     mapps <- aes_string(x="GroupsComparison", y=rlang::as_name(.distmethod), color="GroupsComparison") 
  }else{
     dist.mat <- .data %>%
                 mp_extract_dist(distmethod=rlang::as_name(.distmethod))
     
     sample.hclust <- hclust(dist.mat, method=hclustmethod) %>% ape::as.phylo()

     tbl <- dist.mat %>% as.matrix() %>%
            tibble::as_tibble(rownames="Sample") %>%
            tidyr::pivot_longer(cols=-"Sample", 
                                names_to="col.samples", 
                                values_to=rlang::as_name(.distmethod))

     mapps <- aes_string(x="Sample", 
                         y="col.samples", 
                         color=rlang::as_name(.distmethod), 
                         size=rlang::as_name(.distmethod))

     if (!rlang::quo_is_null(.group)){
         sampleda <- .data %>% 
                     mp_extract_sample() %>% 
                     dplyr::select(!!rlang::sym("Sample"), !!.group)
         gp <- quo.vect_to_str.vect(.group) %>% unique()
         tbl %<>% dplyr::left_join(sampleda, by="Sample", suffix = c("", ".y"))
         gh.layers <- list()
         for (i in seq_len(length(gp))){
             gh.layers[[i]] <- ggplot2::geom_tile(data=sampleda, 
                                      mapping=aes(x=!!rlang::sym("Sample"), y=gp[i], fill=!!rlang::sym(gp[i])))
                                    #ggplot2::scale_y_discrete(position="right", expand = c(0, 0), labels=gp[i])
                               
         }
         #gh.layers <- unlist(gh.layers, recursive=FALSE)
     }
  }

  p <- ggplot(data=tbl, mapping=mapps)
  xylabs <- labs(x=NULL, y=NULL)

  if (!rlang::quo_is_null(.group) && group.test){
      if ("color" %in% names(params)){
          color <- params$color
          params$color <- NULL
      }else if( "colour" %in% names(params)){
          color <- params$colour
          params$color <- NULL
      }else{
          color <- "black"
      }
      p <- p + 
           labs(x=NULL)  
      signlayer <- do.call(geom_signif, 
                           c(list(comparisons=comparisons), 
                             test=test, 
                             step_increase=step_increase, 
                             color=color, 
                             params))
      p <- p + signlayer +
           ggplot2::geom_jitter(
               mapping = aes_string(fill="GroupsComparison"),
               position = ggplot2::position_jitter(seed=123), 
               shape = 21,
               color = "black",
               alpha = 0.8
           ) +
           ggplot2::geom_boxplot(fill=NA) +
           theme_taxbar(legend.position="none")
  }else{

      p <- p + ggplot2::geom_point() + xylabs +
           ggplot2::scale_color_viridis_c() +
           ggplot2::scale_y_discrete(position="right") + 
           theme(axis.text = element_text(size=8),
                 axis.text.x=element_text(angle=-45, hjust = 0),
           )
      indexname <- c(p$labels$colour)
      if (!rlang::quo_is_null(.group)){
          for (i in seq_len(length(gh.layers))){
              f <- ggplot() + gh.layers[i] + xylabs 
              f.left <- f + 
                        ggplot2::coord_flip(expand= FALSE) + 
                        theme(axis.text=element_blank(), axis.ticks=element_blank())

              f.top <- f +
                       ggplot2::scale_y_discrete(position="right", expand = c(0, 0), labels=gp[i]) + 
                       theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
              p %<>% insert_top(f.top, height = 0.04 + 0.0025 * (i - 1))
              p %<>% insert_left(f.left, width = 0.04 + 0.0025 * (i - 1))
          }
          indexname <- c(indexname, rep(gp, each = 2))
      }
      p2 <- ggtree(sample.hclust, branch.length="none", size=.8)
      p <- p %>% insert_left(p2, width=0.12)
      p <- p %>% insert_top(p2 + ggtree::layout_dendrogram(), height=0.12)
      p$index <- c(indexname, paste0("tree", seq_len(2)))
      p %<>% add_class("aplot.heatmap")
  }
  return(p)
}

#' @rdname mp_plot_dist-methods
#' @aliases mp_plot_dist,MPSE
#' @export mp_plot_dist
setMethod("mp_plot_dist", signature(.data="MPSE"), .internal_plot_dist)

#' @rdname mp_plot_dist-methods
#' @aliases mp_plot_dist,tbl_mpse
#' @export mp_plot_dist
setMethod("mp_plot_dist", signature(.data="tbl_mpse"), .internal_plot_dist)

#' @rdname mp_plot_dist-methods
#' @aliases mp_plot_dist,grouped_df_mpse
#' @export mp_plot_dist
setMethod("mp_plot_dist", signature(.data="grouped_df_mpse"), .internal_plot_dist)

#' Plotting the result of PCA, PCoA, CCA, RDA, NDMS or DCA
#' @rdname mp_plot_ord-methods
#' @param .data MPSE or tbl_mpse object, it is required.
#' @param .ord a name of ordination (required), options are PCA, PCoA, DCA, NMDS, RDA, CCA,
#' but the corresponding calculation methods (mp_cal_pca, mp_cal_pcoa, ...) 
#' should be done with action="add" before it.
#' @param .dim integer which dimensions will be displayed, it should be a vector (length=2) 
#' default is c(1, 2). if the length is one the default will also be displayed.
#' @param .group the column name of variable to be mapped to the color of points (fill character 
#' of \code{geom_star}) or one specified color code, default is NULL, meaning \code{fill}=NA,
#' the points are hollow.
#' @param .starshape the column name of variable to be mapped to the shapes of points (starshape 
#' character of \code{geom_star}) or one specified \code{starshape} of point of \code{ggstar},
#' default is NULL, meaning \code{starshape}=15 (circle point).
#' @param .size the column name of variable to be mapped to the size of points (size character of 
#' \code{geom_star}) or one specified \code{size} of point of \code{ggstar}, default is NULL,
#' meaning the \code{size}=1.5, the size of points.
#' @param .alpha the column name of variable to be mapped to the transparency of points (alpha
#' character of \code{geom_star}) or one specified \code{alpha} of point of \code{ggstar}.
#' default is NULL, meaning the \code{alpha}=1, the transparency of points.
#' @param .color the column name of variable to be mapped to the color of line of points (color 
#' character of \code{geom_star}) or one specified \code{starshape} of point of \code{ggstar},
#' default is NULL, meaning the \code{color} is 'black'.
#' @param starstroke numeric the width of edge of points, default is 0.5.
#' @param show.side logical whether display the side boxplot with the specified \code{.dim} 
#' dimensions, default is TRUE.
#' @param show.adonis logical whether display the result of \code{mp_adonis} with \code{action='all'},
#' default is FALSE.
#' @param ellipse logical, whether to plot ellipses, default is FALSE. (.group or .color variables 
#' according to the 'geom', the default geom is path, so .color can be mapped to the corresponding 
#' variable).
#' @param show.sample logical, whether display the sample names of points, default is FALSE.
#' @param show.envfit logical, whether display the result after run [mp_envfit()], default is FALSE.
#' @param p.adjust a character method of p.adjust \code{\link[stats]{p.adjust}}, default is NULL,
#' options are 'fdr', 'bonferroni', 'BH' etc.
#' @param filter.envfit logical or numeric, whether to remove the no significant environment factor after 
#' run [mp_envfit()], default is FALSE, meaning do not remove. If it is numeric, meaning the keep p.value
#' or the adjust p with \code{p.adjust} the factors smaller than the numeric, e.g when filter.envfit=0.05 or
#' (filter.envfit=TRUE), meaning the factors of p <= 0.05 will be displayed.
#' @param ... additional parameters, see also the \code{\link[ggplot2]{stat_ellipse}}.
#' @seealso [mp_cal_pca()], [mp_cal_pcoa], [mp_cal_nmds], [mp_cal_rda], 
#' [mp_cal_cca], [mp_envfit()] and [mp_extract_internal_attr()]
#' @export 
#' @examples
#' \dontrun{
#' library(vegan)
#' data(varespec, varechem)
#' mpse <- MPSE(assays=list(Abundance=t(varespec)), colData=varechem)
#' envformula <- paste("~", paste(colnames(varechem), collapse="+")) %>% as.formula
#' mpse %<>%
#' mp_cal_cca(.abundance=Abundance, .formula=envformula, action="add") %>%
#' mp_envfit(.ord=CCA, .env=colnames(varechem), permutations=9999, action="add")
#' mpse
#' p1 <- mpse %>% mp_plot_ord(.ord=CCA, .group=Al, .size=Mn)
#' p1
#' p2 <- mpse %>% mp_plot_ord(.ord=CCA, .group=Al, .size=Mn, show.sample=TRUE)
#' p2
#' p3 <- mpse %>% mp_plot_ord(.ord=CCA, .group="blue", .size=Mn, .alpha=0.8, show.sample=TRUE)
#' p3
#' p4 <- mpse %>% mp_plot_ord(.ord=CCA, .group=Al, .size=Mn, show.sample=TRUE, show.envfit=TRUE)
#' p4
#' }
setGeneric("mp_plot_ord", function(
    .data,
    .ord,
    .dim = c(1, 2),
    .group = NULL,
    .starshape = 15,
    .size  = 2,
    .alpha = 1,
    .color = "black",
    starstroke = 0.5,
    show.side = TRUE,
    show.adonis = FALSE,
    ellipse = FALSE,
    show.sample = FALSE,
    show.envfit = FALSE,
    p.adjust = NULL,
    filter.envfit = FALSE,
    ...    
)
  standardGeneric("mp_plot_ord")
)

.internal_plot_ord <- function(
    .data, 
    .ord,
    .dim = c(1, 2), 
    .group = NULL, 
    .starshape = 15,
    .size  = 2,
    .alpha = 1, 
    .color = "black",
    starstroke = 0.5,
    show.side = TRUE,
    show.adonis = FALSE,
    ellipse = FALSE,
    show.sample = FALSE,
    show.envfit = FALSE,
    p.adjust = NULL,
    filter.envfit = FALSE,
    ...){

    .ord <- rlang::enquo(.ord)
    if (rlang::quo_is_missing(.ord)){
        intnms <- .data %>% attr("internal_attr") %>% names() %>% toupper()
        .ord <- intersect(intnms, c("PCA", "PCOA", "NMDS", "RDA", "CCA", "DCA")) %>% 
                magrittr::extract(1) %>%
                rlang::sym()
    } 
    .ord %<>% rlang::as_name() %>% 
            gsub("\\s+", "", .) %>%
            toupper()

    .group <- rlang::enquo(.group)
    .starshape <- rlang::enquo(.starshape)
    .size <- rlang::enquo(.size)
    .alpha <- rlang::enquo(.alpha)
    .color <- rlang::enquo(.color)
    paramlist <- list(...)

    .ord %<>% match.arg(c("PCA", "PCOA", "NMDS", "RDA", "CCA", "DCA"))

    if (length(.dim)>2){
       message("The length of .dim can not be larger than 2.")
       return()
    }else if (length(.dim)==1){
       .dim <- c(1, 2) + 1
    }else{
       .dim <- .dim + 1
    }

    sampleda <- .data %>% mp_extract_sample()

    tbl <- .data %>% 
           mp_extract_internal_attr(name=!!rlang::sym(.ord)) %>%
           suppressMessages() %>%
           tidydr() %>% 
           dplyr::select(c(1, .dim)) %>%
           dplyr::rename(Sample="sites") %>%
           dplyr::inner_join(sampleda, by="Sample", suffix=c("", ".y"))

    maps <- aes(x=!!rlang::sym(colnames(tbl)[2]), 
                y=!!rlang::sym(colnames(tbl)[3]))
    params <- list(starstroke=starstroke)
    if (!rlang::quo_is_null(.group)){
       if (is.colors(.group)){
          params <- c(fill=rlang::as_name(.group), params)
       }else{
          maps <- modifyList(maps, aes(fill=!!.group))
       }
    }else{
       params <- c(fill=NA, params)
    }

    if (!rlang::quo_is_null(.color)){
       if (is.colors(.color)){
          params <- c(color=rlang::as_name(.color), params)
       }else{
          maps <- modifyList(maps, aes(color=!!.color))
       }
    }else{
       params <- c(color="black", params)
    }

    if (!rlang::quo_is_null(.starshape)){
       if (is.quo.numeric(.starshape) || is.starshape(.starshape)){
          starshape <- ifelse(is.quo.numeric(.starshape),
                              is.quo.numeric(.starshape, check=FALSE), 
                              rlang::as_name(.starshape))
          params <- c(starshape=starshape, params)
       }else{
          maps <- modifyList(maps, aes(starshape=!!.starshape))
       }
    }else{
       params <- c(starshape=15, params)
    }

    if (!rlang::quo_is_null(.size)){
       if (is.quo.numeric(.size)){
          params <- c(size = is.quo.numeric(.size, check=FALSE), params)
       }else{
          maps <- modifyList(maps, aes(size=!!.size)) 
       }
    }else{
       params <- c(size = 1.5, params)
    }

    if (!rlang::quo_is_null(.alpha)){
       if (is.quo.numeric(.alpha)){
          params <- c(alpha = is.quo.numeric(.alpha, check=FALSE), params)
       }else{
          maps <- modifyList(maps, aes(alpha=!!.alpha))
       }
    }else{
       params <- c(alpha = 1, params)
    }

    p <- ggplot(data=tbl, mapping=maps)
    ggstar <- "ggstar"
    require(ggstar, character.only=TRUE) %>% suppressMessages()
    if ("fill" %in% names(maps)){
        params <- c(params, color = "gray32")
    }
    point.layer <- do.call(geom_star, params)
    p <- p +
         ggplot2::geom_vline(xintercept=0, color="grey20", linetype=2) +
         ggplot2::geom_hline(yintercept=0, color="grey20", linetype=2) +
         point.layer +
         theme_bw() +
         theme(panel.grid=element_blank())

    labelparams <- list(mapping=aes_(label=~Sample), size=3, segment.size=0.25)
    labelparams <- extract_params(labelparams,
                                  paramlist,
                                  c(ggrepel::GeomTextRepel$default_aes %>% names(),
                                    ggrepel::geom_text_repel %>% 
                                        as.list() %>% 
                                        names() %>% 
                                        setdiff(., c("mapping", "data", "", "..."))
                                  )
                         )
    if (show.side){
       if ("fill" %in% names(maps) && !is.discrete(tbl, maps, "fill")){
           side.y <- do.call(ggside::geom_xsideboxplot, 
                             list(mapping=aes(y=!!.group), color = "black", orientation = "y", show.legend = FALSE)) %>%
                     suppressMessages()
           
           side.x <- do.call(ggside::geom_ysideboxplot,
                             list(mapping=aes(x=!!.group), color = "black",orientation = "x", show.legend = FALSE)) %>%
                     suppressMessages()

       }else if ("color" %in% names(maps) && !is.discrete(tbl, maps, "color")){
           side.y <- do.call(ggside::geom_xsideboxplot, 
                             list(mapping=aes(y=!!.color), fill=NA, orientation="y", show.legend = FALSE)) %>%
                     suppressMessages()

           side.x <- do.call(ggside::geom_ysideboxplot,
                             list(mapping=aes(x=!!.color), fill=NA, orientation="x", show.legend = FALSE)) %>%
                     suppressMessages()

       }else if ("starshape" %in% names(maps) && !is.discrete(tbl, maps, "starshape")){
           side.y <- do.call(ggside::geom_xsideboxplot, 
                             list(mapping=aes(y=!!.starshape), orientation="y", show.legend = FALSE)) %>%
                     suppressMessages()

           side.x <- do.call(ggside::geom_ysideboxplot, 
                             list(mapping=aes(x=!!.starshape), orientation="x", show.legend = FALSE)) %>%
                     suppressMessages()
       }else{
           side.y <- NULL
           side.x <- NULL
       }

       p <- suppressMessages(
            p +
            side.y + 
            ggside::scale_xsidey_discrete() + 
            side.x +
            ggside::scale_ysidex_discrete(guide = ggplot2::guide_axis(angle=-45)) + 
            theme(ggside.panel.scale = 0.25)
            )
    }

    if (show.sample){
        text.layer <- do.call(ggrepel::geom_text_repel, labelparams) 
        p <- p + text.layer
    }

    if (ellipse){
        if (length(paramlist)>0){
            ellipse.layer <- do.call(ggplot2::stat_ellipse, paramlist)
        }else{
            ellipse.layer <- ggplot2::stat_ellipse()
        }
        p <- p + ellipse.layer
    }

    if (show.envfit){
        envfit.res <- check.envfit(x=.data, ord=.ord)
        if (!is.null(envfit.res)){
            tbl2 <- envfit.res %>% mp_fortify() 
            if (!is.null(p.adjust)){
                tbl2 %<>% 
                    dplyr::mutate(pvals=stats::p.adjust(.data$pvals, method=p.adjust)) 
            }
            tbl2 %<>% 
                dplyr::mutate(label = 
                              dplyr::case_when(
                                  .data$pvals <= 0.05 & .data$pvals > 0.01 ~ paste0(.data$label, " *"),
                                  .data$pvals <= 0.01 & .data$pvals > 0.001 ~ paste0(.data$label, " **"),
                                  .data$pvals <= 0.001 ~ paste0(.data$label, " ***"),
                                  TRUE ~ .data$label
                              )
                )

            if (filter.envfit){
                tbl2 %<>% dplyr::filter(.data$pvals <= 0.05)
            }else if(is.numeric(filter.envfit)){
                tbl2 %<>%  dplyr::filter(.data$pvals <= filter.envfit)
            }

            if (nrow(tbl2)==0){
                message("No signature environment factor was found!")
                return(p) 
            }
            arrows.multi <- .internal.cal.arrows.multi(arrows=tbl2, data=tbl, dim=.dim)
            tbl2[, .dim] <- arrows.multi * tbl2[, .dim]
            nms <- colnames(tbl2)
            maps2 <- aes(x=0, y=0, 
                         xend=!!rlang::sym(nms[.dim[1]]), 
                         yend=!!rlang::sym(nms[.dim[2]]))
            maps3 <- aes(x=!!rlang::sym(nms[.dim[1]]), 
                         y=!!rlang::sym(nms[.dim[2]]), 
                         label=!!rlang::sym("label"))

            p <- p + 
                 ggplot2::geom_segment(
                    data = tbl2, 
                    mapping = maps2, 
                    arrow=ggplot2::arrow(length = unit(0.02, "npc")),
                    inherit.aes = FALSE
                 ) 
            labelparams <- modifyList(labelparams, list(data=tbl2, mapping=maps3, inherit.aes=FALSE))

            lab.layer <- do.call(ggrepel::geom_text_repel, labelparams) 
            p <- p + lab.layer
        }else{
            message("Please check whether the 'mp_envfit' has been performed with action='add'!")
            return(p)
        }
    }else{
       if (.ord %in% c("RDA", "CCA")){
           biplot.da <- tbl %>% attr("biplot")
           centroids.da <- tbl %>% attr("centroids")
           if (!is.null(biplot.da) && !is.null(centroids.da)){
               biplot.da %<>% dplyr::filter(!.data$biplot %in% centroids.da[,1])
           }
           if (!is.null(biplot.da) && nrow(biplot.da) > 0){
               nms <- colnames(biplot.da)
               arrows.multi <- .internal.cal.arrows.multi(arrows=biplot.da, data=tbl, dim=.dim)
               biplot.da[, .dim] <- arrows.multi * biplot.da[, .dim]
               maps2 <- aes(x = 0, 
                            y = 0,
                            xend = !!rlang::sym(nms[.dim[1]]),
                            yend = !!rlang::sym(nms[.dim[2]]))
               maps3 <- aes(x = !!rlang::sym(nms[.dim[1]]),
                            y = !!rlang::sym(nms[.dim[2]]),
                            label = !!rlang::sym(nms[1]))
               biplot.layer <- ggplot2::geom_segment(data=biplot.da, 
                                                     mapping = maps2, 
                                                     arrow=ggplot2::arrow(length = unit(0.02, "npc")), 
                                                     inherit.aes = FALSE)
               labelparams <- modifyList(labelparams, list(data=biplot.da, mapping=maps3, inherit.aes=FALSE))
               lab.layer <- do.call(ggrepel::geom_text_repel, labelparams)
               p <- p + biplot.layer + lab.layer
           }
           if (!is.null(centroids.da) && nrow(centroids.da) > 0){
               nms <- colnames(centroids.da) 
               maps3 <- aes(x = !!rlang::sym(nms[.dim[1]]), 
                            y = !!rlang::sym(nms[.dim[2]]), 
                            label = !!rlang::sym(nms[1])
                        )
               labelparams <- modifyList(labelparams, 
                                         list(data = centroids.da, 
                                              mapping = maps3, 
                                              hjust = 0.5, 
                                              vjust = 0.5, 
                                              nudge_x = 0, 
                                              nudge_y = 0,
                                              inherit.aes = FALSE))
               lab.layer <- do.call(ggplot2::geom_text, labelparams)
               p <- p + lab.layer
           }
       }
    }

    if (show.adonis){
        p <- .add_adonis_layer(plot = p, data=.data, show.side = show.side)
    }

    return(p)
}

.add_adonis_layer <- function(plot, data, show.side){
   adonis.res <- data %>% mp_extract_internal_attr(name='adonis') 
   if (is.null(adonis.res)){
       cli::cli_inform(c(
         "The {.cls {as.character(class(data))}} does not contain the result of {.arg adonis}. Please make sure ",
         "the {.fn mp_adonis} had been run with {.arg action = 'add'}."
         )
       ) 
       return(plot)
   }
   
   if (show.side){
       vjust <- 'outward'
       hjust <- 'outward'
   }else{
       vjust <- hjust <- 'inward'
   }
   
   label.tmp <- deparse(
                  bquote(
                    atop(
                     italic("Adonis")~":",
                    atop(
                     italic("R")^2~"="~.(round(adonis.res$R2[[1]], 4)),
                     italic("p")~"  ="~.(round(adonis.res$`Pr(>F)`[[1]], 4))
                    )
                  )
               ),
               width.cutoff = 200
             )   

   df.text <- data.frame(name = label.tmp, x = 1, y =1)
   plot <- plot +
           ggpp::geom_text_npc(
             data = df.text,
             mapping = aes(npcx = .data$x, npcy = .data$y, label = .data$name),
             hjust = hjust,
             vjust = vjust,
             parse = TRUE
           ) +
           ggplot2::coord_cartesian(clip = 'off')
   return(plot)
}

#' @rdname mp_plot_ord-methods
#' @aliases mp_plot_ord,MPSE
#' @export mp_plot_ord
setMethod("mp_plot_ord", signature(.data="MPSE"), .internal_plot_ord)

#' @rdname mp_plot_ord-methods
#' @aliases mp_plot_ord,tbl_mpse
#' @export mp_plot_ord
setMethod("mp_plot_ord", signature(.data="tbl_mpse"), .internal_plot_ord)

#' @rdname mp_plot_ord-methods
#' @aliases mp_plot_ord,grouped_df_mpse
#' @export mp_plot_ord
setMethod("mp_plot_ord", signature(.data="grouped_df_mpse"), .internal_plot_ord)

#' @title adjust the color of heatmap of mp_plot_dist
#' @param .data the plot of heatmap of mp_plot_dist
#' @param x the scale or theme
#' @param aes_var character the variable (column) name of color or size.
#' @export
set_scale_theme <- function(.data, x, aes_var){
    aes_var <- rlang::enquo(aes_var)
    if (inherits(.data, "aplot.heatmap")){
       index <- which(.data$index %in% rlang::as_name(aes_var))
       if (length(index)==0){
           return(.data)
       }else{
           for (i in index){
               .data$plotlist[[i]] <- .data$plotlist[[i]] + x
           }
       }
    }
    return(.data)
}


# this is refer to the ggvegan
.internal.cal.arrows.multi <- function(arrows, data, dim, at = c(0, 0), fill = 0.75){
    u <- c(range(data[, dim[1]], range(data[, dim[2]])))
    u <- u - rep(at, each = 2)
    r <- c(range(arrows[, dim[1]], na.rm = TRUE), range(arrows[, dim[2]],
        na.rm = TRUE))
    rev <- sign(diff(u))[-2]
    if (rev[1] < 0)
        u[1:2] <- u[2:1]
    if (rev[2] < 0)
        u[3:4] <- u[4:3]
    u <- u/r
    u <- u[is.finite(u) & u > 0]
    fill * min(u)
}


quo.vect_to_str.vect <- function(quoexpr){
   strvect <- rlang::quo_text(quoexpr) %>%
                 gsub("c\\(", "", .) %>%
                 gsub("\\)", "", .) %>%
                 base::strsplit(",") %>%
                 unlist() %>%
                 gsub("\\s+", "",.) %>%
                 gsub("\"", "", .) %>%
                 gsub("\'", "", .)
   return(strvect)
}

is.colors <- function(x) {
   x <- rlang::as_name(x)
   lapply(x, function(i) {
       tryCatch(is.matrix(grDevices::col2rgb(i)), 
                error = function(e) FALSE)
   }) %>% unlist()
}

is.quo.numeric <- function(quoexpr, check=TRUE){
   xx <- rlang::quo_text(quoexpr) %>% 
       gsub("\"", "", .) %>% 
       gsub("\'", "", .) %>%
       as.numeric() %>%
       suppressWarnings() 
   if (check){
       xx %<>% anyNA() %>%
             magrittr::not()
   }
   return(xx)
}

is.starshape <- function(quoexpr){
   x <- tryCatch(rlang::quo_text(quoexpr) %>% 
          magrittr::is_in(names(starshapes)),
          error = function(e) FALSE
        )
   return(x)
}

is.discrete <- function(data, mapping, aesname){
   data %>% pull(!!mapping[[aesname]]) %>% is.numeric()
}

check.envfit <- function(x, ord){
   if (!ord %in% c("NMDS", "CCA", "RDA", "DCA")){
       message("The mp_envfit only accept the result of CCA, RDA, NMDS or DCA!")
       return(NULL)
   }
   nm <- paste0(ord, "_ENVFIT")

   x <- tryCatch(x %>% 
                 mp_extract_internal_attr(name=!!rlang::sym(nm)) %>% 
                 suppressMessages(),
                 error = function(e) NULL
       )
   return (x)
}

extract_params <- function(originparam, inputparam, defaultparam){
    if (any(defaultparam %in% names(inputparam))){
        args <- intersect(defaultparam, names(inputparam))
        originparam <- modifyList(originparam, inputparam[names(inputparam) %in% args])
    }
    return (originparam)
}

starshapes <- getFromNamespace("starshape_table", "ggstar")
insert_top <- getFromNamespace("insert_top", "aplot")
insert_left <- getFromNamespace("insert_left", "aplot")


# #' set the theme of ggplot object with the STAMP style.
# #' @param colour character the color of theme stamp.
# #' @param axis character which grid of axis will be filled, default is 'y'.
# #' @param ... additional parameter, see also 'theme' of 'ggplot2'.
# #' @keywords internal
# theme_stamp <- function(colour=c('white', 'grey90'), axis = 'y', ...){
#     params <- list(...)
#     axis <- match.arg(axis, c('x', 'y'))
#     if ('color' %in% names(params)){
#         colour <- params$color
#         params$color <- NULL
#     }
#     if (length(colour)!=2){
#         message('The colour is not a vector contained two length.')
#         #colour <- c('white', 'grey90')
#     }
#     structure(
#       list(
#         colour = colour, 
#         axis = axis, 
#         params = params
#       ), 
#       class = 'theme_stamp'
#     )
# }
# 
# #' @method ggplot_add theme_stamp 
# #' @importFrom ggplot2 element_line geom_tile 
# ggplot_add.theme_stamp <- function(object, plot, object_name){
#     gb <- ggplot2::ggplot_build(plot)
#     axis <- paste0('panel_scales_', object$axis)
#     df <- data.frame(AXIS=gb$layout[[axis]][[1]]$get_labels())
#     len.ind <- length(object$colour)
#     axis.num <- nrow(df)
#     df$GROUP.GRID <- rep(object$colour, ceiling(axis.num/len.ind))[seq_len(axis.num)]
#     if (object$axis == 'y'){
#         grid.tile <- geom_tile(
#                        data = df,
#                        mapping = aes(x = 1, 
#                                      y = !!as.symbol("AXIS"), 
#                                      fill = I(!!as.symbol("GROUP.GRID")), 
#                                      height = 1,
#                                      width=Inf
#                                 ),
#                        inherit.aes = FALSE
#                      )
#     }else{
#         grid.tile <- geom_tile(
#                        data = df, 
#                        mapping = aes(x = !!as.symbol("AXIS"),
#                                      y = 1, 
#                                      fill = I(!!as.symbol("GROUP.GRID")), 
#                                      height = Inf,
#                                      width = 1
#                                  ), 
#                        inherit.aes = FALSE
#                      )
#     }
#     plot <- plot + ggnewscale::new_scale_fill() + grid.tile
#     plot$layers <- c(plot$layers[[length(plot$layers)]], plot$layers[-length(plot$layers)])
#     axis.keep <- paste0('axis.line.', setdiff(c('x', 'y'), object$axis))
#     default.theme <- list(element_blank(), element_line())
#     names(default.theme) <- c('panel.background', axis.keep)
#     if (axis.keep %in% names(object$params)){
#         object$params <- c(object$params, default.theme[[-2]])
#     }else{
#         object$params <- c(object$params, default.theme)
#     }
#     th <- do.call("theme", object$params)
#     plot <- plot + th
#     return(plot)
# }

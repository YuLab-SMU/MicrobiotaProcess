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
#' @param ... additional parameters, meaningless now.
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
                                ...
                                ){
    .abundance <- rlang::enquo(.abundance)
    .group <- rlang::enquo(.group)
    taxa.class <- rlang::enquo(taxa.class)

    if (rlang::quo_is_null(taxa.class) || 
        (.data %>% mp_extract_tree() %>% is.null() %>% suppressMessages())){
        taxa.class <- rlang::sym("OTU")
    }

    if (rlang::quo_is_missing(.abundance)){
        .abundance <- rlang::sym("Abundance")
    }
    
    if (relative){
         if (force){
             abundance.nm <- paste0("Rel", rlang::as_name(.abundance))
         }else{
             abundance.nm <- "RelRareAbundance"
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
         prefixBy <- paste0("By", rlang::as_name(.group))
         axis.x <- rlang::as_name(.group)
     }else{
         if (force){
             prefixBy <- ""
         }else{
             prefixBy <- "BySample"
         }
         axis.x <- "Sample"
     }
     abundance.nm <- paste0(abundance.nm, prefixBy)
     
     AbundBy <- abundance.nm %>% gsub("^Rel", "", .)
     if (!plot.group && force){
         AbundBy <- paste0(AbundBy, "BySample")
     }

     if (!any(grepl(paste0("^", AbundBy), .data %>% mp_extract_feature() %>% colnames()))){
         if (!rlang::quo_is_null(.group) && plot.group){
             .data %<>% mp_cal_abundance(.abundance=!!.abundance, .group=!!.group, force=force, relative=relative)
         }else{
             .data %<>% mp_cal_abundance(.abundance=!!.abundance, force=force, relative=relative)
         }
     }

     tbl <- .data %>% 
            mp_extract_abundance(taxa.class=!!taxa.class, topn = topn) %>%
            tidyr::unnest(cols=AbundBy) %>% 
            dplyr::rename(!!taxa.class:="label") %>%
            suppressMessages()
     
     if (!plot.group && ((force && !relative) || (!force && !relative))){
         abundance.nm %<>% gsub("BySample", "", .)
     }else if (!plot.group && (force && relative)){
         abundance.nm %<>% paste0("BySample")
     }

     p <- ggplot(data=tbl,
                 mapping=aes_string(x=axis.x,
                                    y=abundance.nm,
                                    alluvium=rlang::as_name(taxa.class),
                                    fill=rlang::as_name(taxa.class))
          ) +
          ggalluvial::geom_flow(stat="alluvium", lode.guidance = "frontback", color = "darkgray") +
          ggalluvial::geom_stratum(stat="alluvium") +
          ggplot2::labs(x=NULL, y=ylabs) +
          scale_fill_manual(values=get_cols(tbl %>% pull(!!taxa.class) %>% unique() %>% length()))

     if (!rlang::quo_is_null(.group) && !plot.group){
         p <- p + facet_grid(cols=ggplot2::vars(!!.group), scales="free_x", space="free")
     }

     p <- p + theme_taxbar()
     
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
#' @param .alpha the column name of alpha index after run mp_cal_alpha
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
#'      mp_plot_alpha(.group=time, .alpha=c(Observe, Shannon, J))
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

    tbl <- .data %>% 
        mp_extract_sample() %>%
        dplyr::select(!!rlang::sym("Sample"), !!.group, !!.alpha) %>%
        tidyr::pivot_longer(cols=-c(rlang::sym("Sample"), !!.group), 
                            names_to="Measure", values_to="Alpha")

    newlevels <- rlang::quo_text(.alpha) %>% 
        gsub("c\\(", "", .) %>% 
        gsub("\\)", "", .) %>% 
        base::strsplit(",") %>% 
        unlist() %>% 
        gsub("\\s+", "", .) %>%
        gsub("\"", "", .) %>%
        gsub("\'", "", .)

    tbl$Measure <- factor(tbl$Measure, levels=newlevels)

    if (!is.null(.group)){
        if (is.null(comparisons)){
           comparisons <- tbl %>% 
                          pull(!!.group) %>% 
                          unique() %>% 
                          utils::combn(2) %>% 
                          apply(2, list) %>% 
                          unlist(recursive = FALSE)
        }
        mapping <- aes_string(x = rlang::as_name(.group), 
                              y = "Alpha", 
                              fill = rlang::as_name(.group))
    }else{
        mapping <- aes_string(x = "Sample", 
                              y = "Alpha"
                   )
    }

    p <- ggplot(data=tbl, mapping = mapping)

    if (!is.null(.group)){
        p <- p +
           gghalves::geom_half_violin(color=NA, side="l", trim=FALSE) +
           ggplot2::geom_boxplot(aes_string(color=rlang::as_name(.group)), 
                                 fill=NA, 
                                 position=ggplot2::position_nudge(x=.22), 
                                 width=0.2) +
           gghalves::geom_half_point(side="r", shape=21) +
           ggsignif::geom_signif(comparisons=comparisons, test=test, ...) +
           ggplot2::facet_wrap(facets=ggplot2::vars(!!rlang::sym("Measure")), scales="free_y", nrow=1) +
           ggplot2::scale_fill_manual(values=get_cols(tbl %>% pull(!!.group) %>% unique() %>% length())) +
           ggplot2::scale_color_manual(values=get_cols(tbl %>% pull(!!.group) %>% unique() %>% length()))
    }else{
        p <- p +
             ggplot2::geom_col() +
             ggplot2::facet_wrap(facets=ggplot2::vars(!!rlang::sym("Measure")), scales="free_y")
    }

    p <- p + labs(x=NULL, y="Alpha Index Value") + theme_taxbar(legend.position="right")
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
#' @param ... additional parameters, meaningless now.
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
        ggVennDiagram::ggVennDiagram()
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
#' @param ... additional parameters, meaningless now.
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
         ggplot(mapping=aes(x=!!.upset)) +
         geom_bar() +
         ggupset::scale_x_upset() +
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
#' @param .alpha the names of alpha index, which should be one or more of Observe, Shannon, 
#' ACE, Chao1, Simpson, J, default is Observe.
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
                                         .alpha,
                                         .group = NULL, 
                                         nrow = 1,
                                         plot.group = FALSE,
                                         ...
                                         )
    standardGeneric("mp_plot_rarecurve")
)

.internal_plot_rarecurve <- function(.data, .rare, .alpha = "Observe", .group=NULL, nrow=1, plot.group=FALSE, ...){
    .rare <- rlang::enquo(.rare)
    .alpha <- rlang::enquo(.alpha)
    .group <- rlang::enquo(.group)
    params <- list(...)

    newlevels <- rlang::quo_text(.alpha) %>%
                 gsub("c\\(", "", .) %>%
                 gsub("\\)", "", .) %>%
                 base::strsplit(",") %>%
                 unlist() %>%
                 gsub("\\s+", "",.) %>%
                 gsub("\"", "", .) %>%
                 gsub("\'", "", .)

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
        }else{
            maps <- aes_string(x="readsNums", 
                                  y="value", 
                                  color=rlang::as_name(.group), 
                                  fill=rlang::as_name(.group), 
                                  group="Sample")
        }
    }else{
        maps <- aes_(x=~readsNums, 
                     y=~value, 
                     color=~Sample, 
                     fill=~Sample)
    }

    if ("alpha" %in% names(params)){
        alpha <- params$alpha
        params$alpha <- NULL
    }else{
        alpha <- 0.5
    }

    p <- ggplot2::ggplot() +
         geom_smooth(mapping=maps, data=tbl, method = "lm", formula = y~log(x), alpha=alpha, ...)

    if (!rlang::quo_is_null(.group)){
         p <- p + ggplot2::stat_summary(mapping=maps, data=tbl, fun.data = 'mean_se', geom = "ribbon", alpha = alpha)
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

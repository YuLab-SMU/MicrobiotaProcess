#' Fits an Environmental Vector or Factor onto an Ordination With MPSE or tbl_mpse Object
#' @rdname mp_envfit-methods
#' @param .data MPSE or tbl_mpse object
#' @param .ord a name of ordination, option it is DCA, NMDS, RDA, CCA.
#' @param .env the names of columns of sample group or environment information.
#' @param .dim integer The number of dimensions to be returned, default is 3. 
#' @param action character "add" joins the envfit result to internal attributes of the object, 
#' "only" return a non-redundant tibble with the envfit result. "get" return 'envfit' object can
#' be analyzed using the related vegan funtion. 
#' @param permutations the number of permutations required, default is 999.
#' @param seed a random seed to make the analysis reproducible, default is 123. 
#' @param ... additional parameters see also 'vegan::envfit'
#' @return update object according action
#' @export
#' @author Shuangbin Xu
#' @examples
#' library(vegan)
#' data(varespec, varechem)
#' mpse <- MPSE(assays=list(Abundance=t(varespec)), colData=varechem)
#' envformula <- paste("~", paste(colnames(varechem), collapse="+")) %>% as.formula
#' mpse %<>% 
#'        mp_cal_cca(.abundance=Abundance, .formula=envformula, action="add")
#' mpse2 <- mpse %>%
#'          mp_envfit(.ord=cca, 
#'                    .env=colnames(varechem), 
#'                    permutations=9999, 
#'                    action="add")
#' mpse2 %>% mp_plot_ord(.ord=cca, .group=Al, .size=Mn, show.shample=TRUE, show.envfit=TRUE)
#' \dontrun{
#' tbl <- mpse %>%
#'        mp_envfit(.ord=CCA, 
#'                  .env=colnames(varechem), 
#'                  permutations=9999, 
#'                  action="only")
#' tbl
#' library(ggplot2)
#' library(ggrepel)
#' x <- names(tbl)[grepl("^CCA1 ", names(tbl))] %>% as.symbol()
#' y <- names(tbl)[grepl("^CCA2 ", names(tbl))] %>% as.symbol()
#' p <- tbl %>%
#'      ggplot(aes(x=!!x, y=!!y)) + 
#'      geom_point(aes(color=Al, size=Mn)) + 
#'      geom_segment(data=dr_extract(
#'                             name="CCA_ENVFIT_tb", 
#'                             .f=td_filter(pvals<=0.05 & label!="Humdepth")
#'                        ), 
#'                   aes(x=0, y=0, xend=CCA1, yend=CCA2), 
#'                   arrow=arrow(length = unit(0.02, "npc"))
#'      ) + 
#'      geom_text_repel(data=dr_extract(
#'                               name="CCA_ENVFIT_tb", 
#'                               .f=td_filter(pvals<=0.05 & label!="Humdepth")
#'                           ), 
#'                   aes(x=CCA1, y=CCA2, label=label)
#'      ) +
#'      geom_vline(xintercept=0, color="grey20", linetype=2) +
#'      geom_hline(yintercept=0, color="grey20", linetype=2) +
#'      theme_bw() +
#'      theme(panel.grid=element_blank())
#' p
#' }
setGeneric("mp_envfit", function(.data, .ord, .env, .dim=3, action="only", permutations=999, seed=123, ...)standardGeneric("mp_envfit"))

.internal_cal_envfit <- function(.data, .ord, .env, .dim, action="only", permutations=999, seed=123, ...){

    .ord <- rlang::enquo(.ord) %>%
            rlang::as_name() %>%
            toupper()

    .env <- rlang::enquo(.env)

    .ord %<>% match.arg(c("NMDS", "RDA", "CCA", "DCA"))
    
    ordobj <- .data %>% 
              mp_extract_internal_attr(name=!!rlang::sym(.ord))

    if (is.null(ordobj)){
        ordfun <- switch(.ord,
                         NMDS = "mp_cal_nmds",
                         RDA  = "mp_cal_rda",
                         CCA  = "mp_cal_cca",
                         DCA  = "mp_cal_dca")
        rlang::abort(paste0("The ", rlang::as_name(.ord), " is not present in the object, please run ", ordfun, "before performing mp_envfit."))
    }

    tmpX <- vegan::scores(ordobj, display="site")
    envda <- .data %>%
             mp_extract_sample() %>%
             arrange(match(rownames(tmpX), !!as.symbol("Sample"))) %>%
             tibble::column_to_rownames(var="Sample") %>%
             select(!!.env)

    res <- withr::with_seed(seed, vegan::envfit(ordobj, envda, permutations=permutations, choices=seq_len(.dim), ...))

    if (action=="get"){
        return(res)
    }else if (action=="only"){
        da.ord <- .data %>% 
                 mp_extract_internal_attr(name=!!rlang::sym(.ord)) %>% 
                 tidydr()

        da <- .data %>% 
              mp_extract_sample() %>%
              add_total_attr(da.ord) %>%
              add_attr(res %>% mp_fortify(), name=paste0(.ord, "_ENVFIT_tb"))
        return(da)
    }else if (action=="add"){
        attrnames <- paste0(.ord, "_ENVFIT")
        message("The result of mp_envfit has been saved to the internal attribute of the object !")
        message(paste0("It can be extracted using this-object %>% mp_extract_internal_attr(name='",attrnames, "')"))
        .data %<>%
            add_internal_attr(object=res, name=attrnames)
        return(.data)
    }
}

#' @rdname mp_envfit-methods
#' @aliases mp_envfit,MPSE
#' @exportMethod mp_envfit
setMethod("mp_envfit", signature(.data="MPSE"), .internal_cal_envfit)

#' @rdname mp_envfit-methods
#' @aliases mp_envfit,tbl_mpse
#' @exportMethod mp_envfit
setMethod("mp_envfit", signature(.data="tbl_mpse"), .internal_cal_envfit)

#' @rdname mp_envfit-methods
#' @aliases mp_envfit,grouped_df_mpse
#' @exportMethod mp_envfit
setMethod("mp_envfit", signature(.data="grouped_df_mpse"), .internal_cal_envfit)

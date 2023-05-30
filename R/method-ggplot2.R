# #' Visualizing the MPSE or SummarizedExperiment object with ggplot2
# #' @param data MPSE or SummarizedExperiment object
# #' @param mapping Default list of aesthetic mappings to use for plot.
# #' If not specified, must be supplied in each layer added to the plot.
# #' @param .slot the slot name of the .data object, it can be one of 
# #' 'colData', 'rowData', 'otutree', 'taxatree', 'all', default is 'colData'.
# #' @param ... additional parameters, meaningless now.
# #' @importFrom ggplot2 ggplot aes
# #' @return ggplot-object
# #' @export
# #' @examples
# #' \dontrun{
# #' library(ggtree)
# #' library(ggplot2)
# #' data(mouse.time.mpse)
# #' mpse <- mouse.time.mpse %>%
# #'         mp_rrarefy() %>%
# #'         mp_cal_alpha(.abundance=RareAbundance)
# #' p1 <- ggmpse(
# #'         mpse, 
# #'         mapping = aes(x=time, y = Shannon), 
# #'         .slot=colData
# #'       ) +
# #'       geom_violin(aes(fill=time)) +
# #'       geom_boxplot() +
# #'       geom_jitter(aes(fill=time), color='black')
# #' mpse2 <- mouse.time.mpse %>%
# #'          mp_rrarefy() %>%
# #'          mp_diff_analysis(
# #'            .abundance = RareAbundance, 
# #'            .group = time,
# #'            first.test.alpha = 0.01
# #'          )
# #' p2 <- ggmpse(
# #'         mpse2,
# #'         mapping = aes(x=-log10(fdr), y=OTU, fill=Sign_time),
# #'         .slot = 'rowData'
# #'       ) +
# #'       geom_col(data=td_filter(!is.na(Sign_time)))
# #'}
# ggmpse <- function(data, mapping = aes(), .slot='colData', ...){
#     .slot <- rlang::enquo(.slot)
#     p <- ggplot(data = data, mapping = mapping, .slot = !!.slot, ...)
#     return(p)
# }

#' Differential internal and tip nodes (clades) analysis for MPSE or tbl_mpse object
#' @rdname mp_diff_clade-methods
#' @param .data MPSE or tbl_mpse object
#' @param .abundance the name of abundance to be calculated
#' @param .group the group name of the samples to be calculated.
#' @param .sec.group the second group name of the samples to be calculated.
#' @param action character, "add" joins the new information to the taxatree (if it exists) 
#' and otutree (if it exists) or \code{rowData} and return MPSE object,"only" return a 
#' non-redundant tibble with the result of different analysis. "get" return 'diffAnalysisClass' 
#' object.
#' @param force logical whether to calculate the relative abundance forcibly when the abundance
#' is not be rarefied, default is FALSE.
#' @param relative logical whether calculate the relative abundance, default is TRUE.
#' @param first.test.method the method for first test, option is "kruskal.test", "oneway.test", 
#' "lm", "glm", or "glm.nb", "kruskal_test", "oneway_test" of "coin" package. default is "kruskal.test".
#' @param first.test.alpha numeric the alpha value for the first test, default is 0.05.
#' @param p.adjust character the correction method, default is "fdr", see also p.adjust function
#' default is fdr.
#' @param filter.p character the method to filter pvalue, default is fdr, meanings the features 
#' that fdr <= .first.test.alpha will be kept, if it is set to pvalue, meanings the features that
#' pvalue <= .first.test.alpha will be kept.
#' @param strict logical whether to performed in one-against-one when .sec.group is provided, 
#' default is TRUE (strict).
#' @param fc.method character the method to check which group has more abundance for the 
#' significantly different features, default is "generalizedFC", options are \code{generalizedFC}, 
#' \code{compare_median}, \code{compare_mean}.
#' @param second.test.method the method for one-against-one (the second test), default is "wilcox.test" 
#' other option is one of 'wilcox_test' of 'coin'; 'glm'; 'glm.nb' of 'MASS'.
#' @param second.test.alpha numeric the alpha value for the second test, default is 0.05.
#' @param cl.min integer the minimum number of samples per group for performing test, default is 5.
#' @param cl.test logical whether to perform test (second test) between the groups (the number of sample
#' of the .group should be also larger that cl.min), default is TRUE.
#' @param subcl.min integer the minimum number of samples in each second groups for performing test, 
#' default is 3.
#' @param subcl.test logical whether to perform test for between the second groups (the .sec.group 
#' should be provided and the number sample of each .sec.group should be larger than subcl.min, and 
#' strict is TRUE), default is TRUE.
#' @param ml.method the method for calculating the effect size of features, option is 'lda' or 'rf'.
#' default is 'lda'.
#' @param normalization integer set a big number if to get more meaningful values for the LDA score, 
#' or you can set NULL for no normalization, default is 1000000.
#' @param ldascore numeric the threshold on the absolute value of the logarithmic LDA score, default is 2.
#' @param bootnums integer, set the number of bootstrap iteration for lda or rf, default is 30. 
#' @param sample.prop.boot numeric range from 0 to 1, the proportion of samples for calculating the effect
#' size of features, default is 0.7.
#' @param ci numeric, the confidence interval of effect size (LDA or MDA), default is 0.95.
#' @param seed a random seed to make the analysis reproducible, default is 123.
#' @param type character type="species" meaning the abundance matrix is from the species abundance, other 
#' option is "others", default is "species".
#' @param ... additional parameters
#' @return update object according to the action argument.
#' @author Shuangbin Xu
#' @export
#' @examples
#' \dontrun{
#'   suppressPackageStartupMessages(library(curatedMetagenomicData))
#'   xx <- curatedMetagenomicData('ZellerG_2014.relative_abundance', dryrun=F)
#'   xx[[1]] %>% as.mpse -> mpse
#'   mpse.agg.clade <- mpse %>%
#'     mp_aggregate_clade(
#'       .abundance = Abundance,
#'       force = TRUE,
#'       relative = FALSE,
#'       action = 'add' # other option is 'get' or 'only'.
#'     )
#'   mpse.agg.clade %>% mp_diff_clade(
#'       .abundance = Abundance,
#'       force = TRUE,
#'       relative = FALSE,
#'       .group = disease,
#'       fc.method = "compare_mean"
#'     ) %>%
#'   mp_extract_otutree() %>%
#'   dplyr::filter(!is.na(Sign_disease), keep.td = FALSE)
#' }
setGeneric("mp_diff_clade", function(.data, 
                                     .abundance, 
                                     .group, 
                                     .sec.group = NULL, 
                                     action = "add",
                                     force = FALSE,
                                     relative = TRUE,
                                     first.test.method = "kruskal.test",
                                     first.test.alpha = 0.05,
                                     p.adjust = "fdr",
                                     filter.p = "fdr",
                                     strict = TRUE,
                                     fc.method = "generalizedFC",
                                     second.test.method = "wilcox.test",
                                     second.test.alpha = 0.05,
                                     cl.min = 5,
                                     cl.test = TRUE,
                                     subcl.min = 3,
                                     subcl.test = TRUE,
                                     ml.method = "lda",
                                     normalization = 1000000,
                                     ldascore = 2,
                                     bootnums = 30, 
                                     sample.prop.boot = 0.7,
                                     ci = 0.95, 
                                     seed = 123,
                                     type = "species",
                                     ...
                                     )
     standardGeneric("mp_diff_clade")
)

.internal_mp_diff_clade <-  function(
              .data,
              .abundance,
              .group,
              .sec.group = NULL,
              action = "add",
              force = FALSE,
              relative = TRUE,
              first.test.method = "kruskal.test",
              first.test.alpha = 0.05,
              p.adjust = "fdr",
              filter.p = "fdr",
              strict = TRUE,
              fc.method = "generalizedFC",
              second.test.method = "wilcox.test",
              second.test.alpha = 0.05,
              cl.min = 5,
              cl.test = TRUE,
              subcl.min = 3,
              subcl.test = TRUE,
              ml.method = "lda",
              normalization = 1000000,
              ldascore = 2,
              bootnums = 30,
              sample.prop.boot = 0.7,
              ci = 0.95,
              seed = 123,
              type = "species",
              ...
              ){

     otu.tree <- .data %>% mp_extract_otutree() %>% suppressMessages()

     if (is.null(otu.tree)){
         warning_wrap("The MPSE did not contain otutree slot, which is required by mp_diff_clade.")
     }

     action %<>% match.arg(c("add", "get", "only"))
     .abundance <- rlang::enquo(.abundance)
     .group <- rlang::enquo(.group)
     .sec.group <- rlang::enquo(.sec.group)

     if (rlang::quo_is_missing(.group)){
         rlang::abort("The .group is required, please provide it!")
     }
     if (inherits(.data, "MPSE")){
         assaysvar <- .data %>% SummarizedExperiment::assayNames()
     }else{
         assaysvar <- .data %>% attr("assaysvar")
     }
     sampleda <- .data %>%
                 mp_extract_sample() %>%
                 select(!!as.symbol("Sample"), !!.group, !!.sec.group) %>%
                 tibble::column_to_rownames(var="Sample") %>%
                 dplyr::mutate(across(!!.group, as.factor))


     if (!rlang::quo_is_null(.sec.group)){
         sampleda %<>% 
             duplicatedtaxcheck() %>% 
             tibble::column_to_rownames(var="rowname") %>%
             dplyr::mutate(across(!!.sec.group, as.factor))
         .sec.group <- rlang::as_name(.sec.group)
     }else{
         .sec.group <- NULL
     }
     
     if (relative){
         if (force){
             abundance.nm <- paste0("Rel", rlang::as_name(.abundance), "BySample")
         }else{
             abundance.nm <- "RelRareAbundanceBySample"
         }
     }else{
         if (force){
             abundance.nm <- rlang::as_name(.abundance)
         }else{
             abundance.nm <- "RareAbundance"
         }
     }

     #AbundBy <- abundance.nm %>% gsub("^Rel", "", .)
     if (! abundance.nm %in% c(otu.tree %>% treeio::get.fields())){
         aggregate_fun <- getOption(x = "aggregate_fun", default = "mean")
         .data %<>% mp_aggregate_clade(
            .abundance = !!.abundance, 
            force = force, 
            relative = relative, 
            aggregate_fun = aggregate_fun,
            action = 'add'
         )
     }

     otu.tree <- .data %>% mp_extract_otutree()

     AbundBySample <- abundance.nm #%>% 
                      #gsub("^Rel", "", .) %>%
                      #gsub("BySample", "", .) %>%
                      #paste0("BySample")

     f_tb <- otu.tree %>%
             dplyr::select(c("label", AbundBySample)) %>% 
             suppressMessages() %>%
             tidyr::unnest(cols=AbundBySample) %>%
             dplyr::select(c("label", "Sample", abundance.nm)) %>% 
             distinct() %>%
             tidyr::pivot_wider(id_cols="Sample", names_from="label", values_from=abundance.nm) %>%
             tibble::column_to_rownames(var="Sample")
     #f_tb <- f_tb[, !apply(f_tb, 2, function(x)var(x)==0), drop = FALSE]
     vars <- f_tb %>% colnames()

     datameta <- merge(f_tb, sampleda, by=0) 

     first.res <- multi_compare(fun=first.test.method,
                               data=datameta,
                               feature=vars,
                               factorNames=rlang::as_name(.group)) %>%
                  lapply(get_pvalue)
     
     first.res <- do.call("rbind", first.res) %>% 
                  as.data.frame() %>% 
                  dplyr::rename(pvalue="V1") %>% 
                  dplyr::mutate(f=vars, fdr=p.adjust(.data$pvalue, method=p.adjust)) %>%
                  dplyr::select(c("f", "pvalue", "fdr"))
                   
     first.test.sig.vars <- first.res %>% 
                  dplyr::filter(!!as.symbol(filter.p) <= first.test.alpha & !is.na(!!as.symbol(filter.p))) %>% 
                  dplyr::pull(.data$f)

     if (!length(first.test.sig.vars)>0){
          message("There are not significantly discriminative features after internal first test !")
          return(NULL)
     }

     compareclass <- sampleda %>% 
                    get_classlevels(classgroup = rlang::as_name(.group)) %>% 
                    get_compareclass()
     
     if (!is.null(.sec.group) && strict){
         class2sub <- get_class2sub(sampleda = sampleda, 
                                    classgroup = rlang::as_name(.group), 
                                    subclass = .sec.group)
         comsubclass <- apply(compareclass,1,function(x)get_comparesubclass(x[1],x[2],class2sub))
         second.test.sig.vars <- diffsubclass(datasample = datameta, 
                                    features = first.test.sig.vars, 
                                    comsubclass = comsubclass,
                                    classgroup = rlang::as_name(.group),
                                    subclass = .sec.group,
                                    fcfun = fc.method, 
                                    secondcomfun = second.test.method, 
                                    submin = subcl.min, 
                                    subclwilc = subcl.test, 
                                    pfold = second.test.alpha, 
                                    ...)
     }else{
         second.test.sig.vars <- diffclass(datasample = datameta, 
                                 features = first.test.sig.vars, 
                                 comclass = compareclass, 
                                 classgroup = rlang::as_name(.group), 
                                 fcfun = fc.method,
                                 secondcomfun = second.test.method,
                                 classmin = cl.min,
                                 clwilc = cl.test,
                                 pfold = second.test.alpha, 
                                 ...)
     }

     if (length(second.test.sig.vars)==0){
         message("There are not significantly discriminative features after internal second test!")
         return(NULL)
     }

     leaveclasslevels <- unlist(lapply(names(second.test.sig.vars), 
                                       function(x){unlist(strsplit(x,"-vs-"))}))

     second.test.sig.vars <- get_consistentfeatures(diffsubclassfeature = second.test.sig.vars, 
                                          classgroup = rlang::as_name(.group),
                                          classlevels = leaveclasslevels) 
     second.test.sig.vars.vectors <- second.test.sig.vars %>% get_secondvarlist()
     if (!is.null(normalization)){
         normalization <- normalization / 100
         f_tb <- f_tb * normalization 
     }
     dameta <- merge(f_tb, sampleda, by=0) %>% 
               tibble::column_to_rownames(var = "Row.names") %>%
               select(c(second.test.sig.vars.vectors, rlang::as_name(.group))) %>%
               dplyr::group_split(!!.group)
     
     dameta <- withr::with_seed(seed, get_sampledflist(dameta, bootnums = bootnums, ratio = sample.prop.boot)) %>%
               remove_constant()

     if (ml.method=="lda"){
         ml.res <- LDAeffectsize(
                                dameta, 
                                compareclass, 
                                rlang::as_name(.group), 
                                bootnums = bootnums, 
                                LDA = ldascore,
                                ci = ci)
     }
     if (ml.method=="rf"){
         ml.res <- rfimportance(
                               dameta, 
                               rlang::as_name(.group), 
                               bootnums = bootnums, 
                               effsize = ldascore, 
                               ci = ci)
     }
 
     params <- list(mlfun = ml.method, firstcomfun = first.test.method, secondcomfun = second.test.method,
                    firstalpha = first.test.alpha, filtermod = filter.p, classgroup = rlang::as_name(.group),
                    normalization = normalization, subclass = .sec.group, strictmod = strict, type = type, 
                    fcfun = fc.method, clmin = cl.min, clwilc = cl.test, subclmin = subcl.min, subclwilc = subcl.test)

     result <- merge_total_res(kwres = first.res, secondvars = second.test.sig.vars, mlres = ml.res, params = params)

     if (action=="get"){
         res <- new("diffAnalysisClass", originalD = f_tb, sampleda = sampleda, taxda = NULL, result = result, kwres = first.res,
                   secondvars = second.test.sig.vars, mlres = ml.res, someparams = params)
         return(res)
     }
     result <- .combine_others(result, first.res)
     if (action=="only"){
         return(result)
     }else if (action=="add"){
         newgroup <- paste0("Sign_", rlang::as_name(.group))
         result %<>% dplyr::rename(label="f", !!newgroup := !!.group)
         otu.tree %<>% dplyr::left_join(result, by="label", suffix=c("", ".y"))
         otutree(.data) <- otu.tree
         return(.data)
     }
}

#' @rdname mp_diff_clade-methods
#' @aliases mp_diff_clade,MPSE
#' @exportMethod mp_diff_clade
setMethod("mp_diff_clade", signature(.data="MPSE"), .internal_mp_diff_clade)

#' @rdname mp_diff_clade-methods
#' @aliases mp_diff_clade,tbl_mpse
#' @exportMethod mp_diff_clade
setMethod("mp_diff_clade", signature(.data="tbl_mpse"), .internal_mp_diff_clade)

#' @rdname mp_diff_clade-methods
#' @aliases mp_diff_clade,grouped_df_mpse
#' @exportMethod mp_diff_clade
setMethod("mp_diff_clade", signature(.data="grouped_df_mpse"), .internal_mp_diff_clade)


## #' The visualization of result of mp_diff_analysis
## #' @rdname mp_plot_diff_res-methods
## #' @param .data MPSE or tbl_mpse after run mp_diff_analysis with \code{action="add"}
## #' @param layout the type of tree layout, should be one of "rectangular", "roundrect", "ellipse", 
## #' "circular", "slanted", "radial", "inward_circular".
## #' @param tree.type one of 'taxatree' and 'otutree', taxatree is the taxonomy class tree
## #' 'otutree' is the phylogenetic tree built with the representative sequences.
## #' @param .taxa.class character the name of taxonomy class level, default is NULL, meaning it will
## #' extract the phylum annotation automatically.
## #' @param tiplab.size numeric the size of tiplab, default is 2.
## #' @param offset.abun numeric the gap (width) (relative width to tree) between the tree and abundance 
## #' panel, default is 0.04.
## #' @param pwidth.abun numeric the panel width (relative width to tree) of abundance panel, 
## #' default is 0.3 .
## #' @param offset.effsize numeric the gap (width) (relative width to tree) between the tree and 
## #' effect size panel, default is 0.3 .
## #' @param pwidth.effsize numeric the panel width (relative width to tree) of effect size panel, 
## #' default is 0.5 .
## #' @param group.abun logical whether to display the relative abundance of group instead of sample,
## #' default is FALSE.
## #' @param tiplab.linetype numeric the type of line for adding line if 'tree.type' is 'otutree',
## #' default is 3 .
## #' @param ... additional parameters, meaningless now.
## #' @export
## setGeneric("mp_plot_diff_res", 
##                 function(
##                     .data, 
##                     layout = "radial", 
##                     tree.type = "taxatree", 
##                     .taxa.class = NULL, 
##                     tiplab.size = 2,
##                     offset.abun = 0.04,
##                     pwidth.abun = 0.8,
##                     offset.effsize = 0.3,
##                     pwidth.effsize = 0.5,
##                     group.abun = FALSE,
##                     tiplab.linetype = 3,
##                     ...) 
##                     standardGeneric("mp_plot_diff_res")
## )
## 
## 
## #' @importFrom ggplot2 geom_col
## #' @importFrom ggtreeExtra geom_fruit
## .internal_mp_plot_diff_res <- function(.data,
##                                        layout = "radial",
##                                        tree.type = "taxatree",
##                                        .taxa.class = NULL,
##                                        tiplab.size = 2,
##                                        offset.abun = 0.04,
##                                        pwidth.abun = 0.8,
##                                        offset.effsize = 0.3,
##                                        pwidth.effsize = 0.5,
##                                        group.abun = FALSE,
##                                        tiplab.linetype = 3,
##                                        ...
##                                       ){
##     .taxa.class <- rlang::enquo(.taxa.class)
##     layout %<>% match.arg(c("rectangular", "roundrect", "ellipse", "circular", 
##                             "slanted", "radial", "inward_circular"))
## 
##     tree.type %<>% match.arg(c("taxatree", "otutree"))
## 
##     if (tree.type == 'otutree'){
##         anno.tree <- .data %>% mp_extract_tree(type="otutree") %>% suppressMessages()
##         if (is.null(anno.tree)){
##             stop_wrap("The otutree slot of the MPSE class is empty, you can try to 
##                        select taxatree by setting tree.type to 'taxatree'.")
##         }else{
##             taxada <- .data %>% mp_extract_taxonomy() %>% suppressMessages()
##             if (!is.null(taxada)){
##                 anno.tree %<>% dplyr::left_join(taxada, by=c("label"="OTU"))
##                 if (rlang::quo_is_null(.taxa.class)){
##                     .taxa.class <- rlang::sym(colnames(taxada)[3])
##                 }
##             }
##         }
##     }
##     
##     if (tree.type == "taxatree"){
##         anno.tree <- .data %>% mp_extract_tree() %>% suppressMessages()
##         if (is.null(anno.tree)){
##             stop_wrap("The taxatree slot of the MPSE class is empty, you can try to 
##                        select otutree by setting tree.type to 'otutree'.")
##         }else{
##             if (rlang::quo_is_null(.taxa.class)){
##                 .taxa.class <- anno.tree %>% 
##                                tidytree::filter(!!rlang::sym("nodeDepth")==2, keep.td=FALSE) %>% 
##                                pull(!!rlang::sym("nodeClass")) %>% 
##                                unique()
##                 .taxa.class <- rlang::sym(.taxa.class)
##             }
##         }
##     }
## 
##     nsample <- .data %>% mp_extract_sample() %>% nrow()
##     field.da.nm <- tidytree::get.fields(anno.tree)
##     if ("LDAmean" %in% field.da.nm){
##         x.bar <- 'LDAmean'
##         x.bar.title <- "log10(LDA)"
##     }else{
##         x.bar <- "MDAmean"
##         x.bar.title <- "MDA"
##     }
##     
##     sign.field <- field.da.nm[grep("^Sign_", field.da.nm)][1]
## 
##     group.nm <- gsub("^Sign_", "", sign.field)
##     flag <- grepl(paste0("By", group.nm), field.da.nm)
## 
##     if (nsample > 50 || group.abun){
##         if (!any(flag)){
##             stop_wrap("The relative abundance of each group will be displayed, but the
##                        relative abundance of each group is not calculated, please run 
##                        the mp_cal_abundance specified group argument before !")
##         }
##         abun.col <- field.da.nm[flag]
##     }else{
##         abun.col <- field.da.nm[grepl("BySample", field.da.nm)]
##     }
##     abun.col <- abun.col[1]
##     x.abun.col <- anno.tree %>% 
##                   dplyr::select(!!rlang::sym(abun.col)) %>% 
##                   tidyr::unnest(!!rlang::sym(abun.col)) %>%
##                   colnames()
##     x.abun.col <- x.abun.col[grepl("^Rel", x.abun.col)]
## 
##     if (tree.type == "otutree"){
##         p1 <- ggtree(
##                 anno.tree,
##                 layout = layout,
##                 size = 0.3
##               )
##         if (!is.null(taxada)){
##             p1 <- p1 +
##                   geom_tiplab(
##                     align = TRUE,
##                     size = 0,
##                     linetype = tiplab.linetype
##                   ) +
##                   geom_fruit(
##                     data = td_filter(!!rlang::sym("isTip")),
##                     geom = geom_point,
##                     mapping = aes(colour = !!.taxa.class),
##                     size = 1.5,
##                     offset = 0
##                   )    
##         }
##     }
## 
##     if (tree.type == "taxatree"){
##         p1 <- suppressWarnings(
##                 ggtree(
##                   anno.tree,
##                   layout = layout,
##                   size = 0.3
##                   ) +
##                 geom_point(
##                   data = td_filter(!.data$isTip),
##                   fill = "white",
##                   size = 1,
##                   shape = 21
##                 )
##               )
##         
##         p1 <- suppressWarnings(
##                 p1 +
##                 ggtree::geom_hilight(
##                   data = td_filter(!!rlang::sym("nodeClass")==rlang::as_name(.taxa.class)),
##                   mapping = aes(
##                                 node = !!rlang::sym("node"), 
##                                 fill = !!rlang::sym("label")
##                             )
##                 )
##               )
##     }
## 
##     if (nsample > 50 || group.abun){
##          mapping <- aes(x=!!rlang::sym(group.nm), size = !!rlang::sym(x.abun.col), 
##                         fill = !!rlang::sym(group.nm), subset = !!rlang::sym(x.abun.col) > 0) 
##          n.pwidth <- .data %>%
##                       mp_extract_sample() %>%
##                       dplyr::pull(!!rlang::sym(group.nm)) %>%
##                       unique() %>% length()    
##     }else{
##          mapping <- aes(x=forcats::fct_reorder(!!rlang::sym("Sample"), !!rlang::sym(group.nm), .fun=min),
##                         size = !!rlang::sym(x.abun.col), fill = !!rlang::sym(group.nm), 
##                         fill = !!rlang::sym(group.nm), subset = !!rlang::sym(x.abun.col) > 0
##                     )
##          n.pwidth <- ncol(.data)
##     }
##     p2 <- suppressWarnings(
##             p1 + 
##             ggnewscale::new_scale_fill() +
##             geom_fruit(
##                data = td_unnest(!!rlang::sym(abun.col)),
##                geom = geom_star,
##                mapping = mapping,
##                starshape = 13,
##                starstroke = 0.05,
##                offset = offset.abun,
##                pwidth = pwidth.abun,
##                grid.params = list(linetype=2)
##             ) +  
##             scale_size_continuous(
##                name="Relative Abundance (%)",
##                range = c(.5, 3),
##                guide = guide_legend(override.aes = list(fill="black"))
##             )
##           )
## 
##     if (nsample > 50 || group.abun){
##         p3 <- suppressWarnings(
##                 p2 + 
##                 geom_tiplab(size=tiplab.size, offset = max(p2$data$xmaxtmp, na.rm=TRUE) - 0.98*max(p2$data$x, na.rm=TRUE), 
##                             align = TRUE, 
##                             linetype=NA)
##               )
##     }else{
##         p3 <- suppressWarnings(
##                 p2 + 
##                 geom_tiplab(size=tiplab.size, offset = max(p2$data$xmaxtmp, na.rm=TRUE) - 0.95*max(p2$data$x, na.rm=TRUE))
##               )
##     }
##     # display the LDA of significant OTU.
##     title.height <- 4.4e-06 * sum(p3$data$isTip) 
##     p4 <- suppressWarnings(
##             p3 +
##             ggnewscale::new_scale_fill() +
##             geom_fruit(
##                #data = td_filter(!is.na(!!rlang::sym(x.bar))),
##                geom = geom_col,
##                mapping = aes(
##                              x = !!rlang::sym(x.bar),
##                              fill = !!rlang::sym(sign.field)
##                              ),
##                orientation = "y",
##                offset = offset.effsize,
##                pwidth = pwidth.effsize,
##                axis.params = list(axis = "x",
##                                   title = x.bar.title,
##                                   title.height = title.height,
##                                   title.size = 2,
##                                   text.size = 1.8,
##                                   vjust = 1),
##                grid.params = list(linetype = 2)
##               )
##           )
## 
##     # display the significant (FDR) taxonomy after kruskal.test (default)
##     p5 <- suppressWarnings(
##             p4 +
##             ggnewscale::new_scale("size") +
##             geom_point(
##                data=td_filter(!is.na(!!rlang::sym("fdr"))),
##                mapping = aes(size = -log10(!!rlang::sym("fdr")),
##                               fill = !!rlang::sym(sign.field),
##                              ),
##                shape = 21
##             ) +
##             scale_size_continuous(range=c(1, 3)) #+
##             #scale_fill_manual(values=c("#1B9E77", "#D95F02"))
##           )
##     
##     p6 <- suppressWarnings(
##             p5 + theme(
##                legend.key.height = unit(0.3, "cm"),
##                legend.key.width = unit(0.3, "cm"),
##                legend.spacing.y = unit(0.02, "cm"),
##                legend.text = element_text(size = 7),
##                legend.title = element_text(size = 9),
##               )
##           )
## 
##     return (suppressWarnings(p6))
## }
## 
## #' @rdname mp_plot_diff_res-methods
## #' @aliases mp_plot_diff_res,MPSE
## #' @export mp_plot_diff_res
## setMethod("mp_plot_diff_res", signature(.data='MPSE'), .internal_mp_plot_diff_res)
## 
## #' @rdname mp_plot_diff_res-methods
## #' @aliases mp_plot_diff_res,tbl_mpse
## #' @export mp_plot_diff_res
## setMethod("mp_plot_diff_res", signature(.data="tbl_mpse"), .internal_mp_plot_diff_res)
## 
## #' @rdname mp_plot_diff_res-methods
## #' @aliases mp_plot_diff_res,grouped_df_mpse
## #' @export mp_plot_diff_res
## setMethod("mp_plot_diff_res", signature(.data="grouped_df_mpse"), .internal_mp_plot_diff_res)

.combine_others <- function(x, y){
    y %<>% dplyr::filter(!.data$f %in% as.vector(x$f))
    x %<>% dplyr::bind_rows(y)
    return(x)
}

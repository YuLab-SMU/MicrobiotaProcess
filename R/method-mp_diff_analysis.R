#' Differential expression analysis for MPSE or tbl_mpse object
#' @rdname mp_diff_analysis-methods
#' @param .data MPSE or tbl_mpse object
#' @param .abundance the name of abundance to be calculated
#' @param .group the group name of the samples to be calculated.
#' @param .sec.group the second group name of the samples to be calculated.
#' @param action character, "add" joins the new information to the taxatree (if it exists) 
#' or \code{rowData} and return MPSE object,"only" return a 
#' non-redundant tibble with the result of different analysis. "get" return 'diffAnalysisClass' 
#' object.
#' @param tip.level character the taxa level to be as tip level
#' @param force logical whether to calculate the relative abundance forcibly when the abundance
#' is not be rarefied, default is FALSE.
#' @param relative logical whether calculate the relative abundance.
#' @param taxa.class character if taxa class is not 'all', only the specified taxa class will
#' be identified, default is 'all'.
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
#' data(mouse.time.mpse)
#' mouse.time.mpse %<>%
#'   mp_rrarefy() 
#' mouse.time.mpse
#' mouse.time.mpse %<>%
#'   mp_diff_analysis(.abundance=RareAbundance, 
#'                    .group=time, 
#'                    first.test.alpha=0.01,
#'                    action="add") 
#' library(ggplot2)
#' p <- mouse.time.mpse %>% mp_plot_diff_res()
#' p <- p + 
#'      scale_fill_manual(
#'        aesthetics = "fill_new", # The fill aes was renamed to "fill_new" for the abundance dotplot layer
#'        values = c("skyblue", "orange")
#'      )  + 
#'      scale_fill_manual(
#'        values=c("skyblue", "orange") # The LDA barplot layer
#'      )
#' ### and the fill aes for hight light layer of tree was renamed to 'fill_new_new'
#' p <- p + 
#'      scale_fill_manual(
#'        aesthetics = "fill_new_new",
#'        values = c("#E41A1C", "#377EB8", "#4DAF4A", 
#'                   "#984EA3", "#FF7F00", "#FFFF33", 
#'                   "#A65628", "#F781BF", "#999999")
#'      )
#' p
#' \dontrun{
#'   ### visualizing the differential taxa with cladogram
#'   f <- mouse.time.mpse %>% 
#'        mp_plot_diff_cladogram(
#'          label.size = 2.5, 
#'          hilight.alpha = .3, 
#'          bg.tree.size = .5, 
#'          bg.point.size = 2, 
#'          bg.point.stroke = .25
#'        ) + 
#'        scale_fill_diff_cladogram(
#'          values = c('skyblue', 'orange')
#'        ) +
#'        scale_size_continuous(range = c(1, 4))
#'   f
#' }
setGeneric("mp_diff_analysis", function(.data, 
                                        .abundance, 
                                        .group, 
                                        .sec.group=NULL, 
                                        action="add",
                                        tip.level="OTU",
                                        force=FALSE,
                                        relative=TRUE,
                                        taxa.class = 'all',
                                        first.test.method="kruskal.test",
                                        first.test.alpha=0.05,
                                        p.adjust="fdr",
                                        filter.p="fdr",
                                        strict = TRUE,
                                        fc.method = "generalizedFC",
                                        second.test.method="wilcox.test",
                                        second.test.alpha=0.05,
                                        cl.min=5,
                                        cl.test=TRUE,
                                        subcl.min = 3,
                                        subcl.test=TRUE,
                                        ml.method="lda",
                                        normalization=1000000,
                                        ldascore=2,
                                        bootnums=30, 
                                        sample.prop.boot=0.7,
                                        ci=0.95, 
                                        seed=123,
                                        type="species",
                                        ...
                                        )
     standardGeneric("mp_diff_analysis")
)

#' @importFrom cli cli_inform make_ansi_style
.internal_mp_diff_analysis <-  function(
              .data,
              .abundance,
              .group,
              .sec.group=NULL,
              action="add",
              tip.level="OTU",
              force=FALSE,
              relative=TRUE,
              taxa.class = 'all',
              first.test.method="kruskal.test",
              first.test.alpha=0.05,
              p.adjust="fdr",
              filter.p="fdr",
              strict = TRUE,
              fc.method = "generalizedFC",
              second.test.method="wilcox.test",
              second.test.alpha=0.05,
              cl.min=5,
              cl.test=TRUE,
              subcl.min = 3,
              subcl.test=TRUE,
              ml.method="lda",
              normalization=1000000,
              ldascore=2,
              bootnums=30,
              sample.prop.boot=0.7,
              ci=0.95,
              seed=123,
              type="species",
              ...
              ){

     action %<>% match.arg(c("add", "get", "only"))
     .abundance <- rlang::enquo(.abundance)
     .group <- rlang::enquo(.group)
     .sec.group <- rlang::enquo(.sec.group)
     tip.level <- rlang::enquo(tip.level)
     taxa.class <- rlang::enquo(taxa.class)
     taxa.class <- rlang::as_name(taxa.class)

     if (rlang::quo_is_missing(.group)){
         rlang::abort("The .group is required, please provide it!")
     }
     if (inherits(.data, "MPSE")){
         assaysvar <- .data %>% SummarizedExperiment::assayNames()
     }else{
         assaysvar <- .data %>% attr("assaysvar")
     }
     .data %<>% mp_select_as_tip(tip.level = !!tip.level)
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
             #dplyr::mutate(!!.sec.group:=factor(!!.sec.group, levels=unique(as.vector(!!.sec.group))))
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

     AbundBy <- abundance.nm %>% gsub("^Rel", "", .)
     if (!any(grepl(paste0("^", AbundBy), .data %>% mp_extract_feature() %>% colnames()))){
         .data %<>% mp_cal_abundance(.abundance=!!.abundance, force=force, relative=relative)
     }

     taxatree <- .data %>% mp_extract_tree()

     if (is.null(taxatree)){
         f_tb <- .data %>% mp_extract_assays(.abundance=!!abundance.nm, byRow=FALSE)
     }else{
         AbundBySample <- abundance.nm %>% 
                          gsub("^Rel", "", .) %>%
                          gsub("BySample", "", .) %>%
                          paste0("BySample")
         f_tb <- taxatree %>%
                 as_tibble() %>%
                 dplyr::filter(.data$nodeClass!="Root")
         if (any(taxa.class %in% unique(f_tb$nodeClass))){
             f_tb %<>% dplyr::filter(.data$nodeClass %in% taxa.class)
         }
         f_tb %<>% dplyr::select(c("label", AbundBySample)) %>% 
                 tidyr::unnest(cols=AbundBySample) %>%
                 dplyr::select(c("label", "Sample", abundance.nm)) %>% 
                 distinct() %>%
                 tidyr::pivot_wider(id_cols="Sample", names_from="label", values_from=abundance.nm) %>%
                 tibble::column_to_rownames(var="Sample")
         
     }
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
                  #dplyr::mutate(f=vars) %>%
                  #tibble::column_to_rownames(var="f") %>%
                  dplyr::mutate(f=vars, fdr=p.adjust(.data$pvalue, method=p.adjust)) %>%
                  dplyr::select(c("f", "pvalue", "fdr"))
                   
     first.test.sig.vars <- first.res %>% 
                  dplyr::filter(!!as.symbol(filter.p) <= first.test.alpha & !is.na(!!as.symbol(filter.p))) %>% 
                  dplyr::pull(.data$f)
     
     msg1 <- function(x){
                 paste0("There are not significantly discriminative features after internal ", x, " test !")
             }
     msg2 <- paste0("The result returned was original input object {.cls ", class(.data)[1], "}.")
     msg3 <- "If you do not want to identify the differential features with such a conservative"
     msg4 <- " condition, you can set {.arg filter.p=\"pvalue\"} or {.arg first.test.alpha = .05} or {.arg strict = FALSE}"
     
     if (!length(first.test.sig.vars)>0){
         cli::cli_inform(cli::make_ansi_style("orange")(c(msg1('first'), msg2, msg3, msg4)))
         return(.data)
     }

     compareclass <- sampleda %>% 
                    get_classlevels(classgroup=rlang::as_name(.group)) %>% 
                    get_compareclass()
     
     if (!is.null(.sec.group) && strict){
         class2sub <- get_class2sub(sampleda=sampleda, 
                                    classgroup=rlang::as_name(.group), 
                                    subclass=.sec.group)
         comsubclass <- apply(compareclass,1,function(x)get_comparesubclass(x[1],x[2],class2sub))
         second.test.sig.vars <- diffsubclass(datasample=datameta, 
                                    features=first.test.sig.vars, 
                                    comsubclass=comsubclass,
                                    classgroup=rlang::as_name(.group),
                                    subclass=.sec.group,
                                    fcfun=fc.method, 
                                    secondcomfun=second.test.method, 
                                    submin=subcl.min, 
                                    subclwilc=subcl.test, 
                                    pfold=second.test.alpha, 
                                    ...)
     }else{
         second.test.sig.vars <- diffclass(datasample=datameta, 
                                 features=first.test.sig.vars, 
                                 comclass=compareclass, 
                                 classgroup=rlang::as_name(.group), 
                                 fcfun=fc.method,
                                 secondcomfun=second.test.method,
                                 classmin=cl.min,
                                 clwilc=cl.test,
                                 pfold=second.test.alpha, 
                                 ...)
     }

     if (!all(unlist(lapply(second.test.sig.vars,function(i)nrow(i) > 0)))){
         cli::cli_inform(cli::make_ansi_style("red")(c(msg1('second'), 
                                                       msg2, msg3, msg4, 
                                                       " or {.var second.test.alpha=.05}")))
         return(.data) 
     }

     leaveclasslevels <- unlist(lapply(names(second.test.sig.vars), 
                                       function(x){unlist(strsplit(x,"-vs-"))}))

     second.test.sig.vars <- get_consistentfeatures(diffsubclassfeature=second.test.sig.vars, 
                                          classgroup=rlang::as_name(.group),
                                          classlevels=leaveclasslevels) 
     second.test.sig.vars.vectors <- second.test.sig.vars %>% get_secondvarlist()
     if (!is.null(normalization)){
         normalization <- normalization / 100
         f_tb <- f_tb * normalization 
     }
     dameta <- merge(f_tb, sampleda, by=0) %>% 
               tibble::column_to_rownames(var="Row.names") %>%
               select(c(second.test.sig.vars.vectors, rlang::as_name(.group))) %>%
               dplyr::group_split(!!.group)
     
     dameta <- withr::with_seed(seed, get_sampledflist(dameta, bootnums=bootnums, ratio=sample.prop.boot)) %>%
               remove_constant()

     if (ml.method=="lda"){
         ml.res <- LDAeffectsize(
                                dameta, 
                                compareclass, 
                                rlang::as_name(.group), 
                                bootnums=bootnums, 
                                LDA=ldascore,
                                ci=ci)
     }
     if (ml.method=="rf"){
         ml.res <- rfimportance(
                               dameta, 
                               rlang::as_name(.group), 
                               bootnums=bootnums, 
                               effsize=ldascore, 
                               ci=ci)
     }
 
     params <- list(mlfun=ml.method, firstcomfun=first.test.method, secondcomfun=second.test.method,
                    firstalpha=first.test.alpha, filtermod=filter.p, classgroup=rlang::as_name(.group),
                    normalization=normalization, subclass=.sec.group, strictmod=strict, type=type, 
                    fcfun=fc.method, clmin=cl.min, clwilc=cl.test, subclmin=subcl.min, subclwilc=subcl.test)

     result <- merge_total_res(kwres=first.res, secondvars=second.test.sig.vars, mlres=ml.res, params=params)
      

     if (action=="get"){
         taxda <- .data %>% 
                  mp_extract_taxonomy() %>%
                  tibble::column_to_rownames(var="OTU")
         taxda[[rlang::as_name(tip.level)]] <- rownames(taxda)
         indx <- match(rlang::as_name(tip.level), colnames(taxda))
         taxda %<>% select(seq_len(indx))
         res <- new("diffAnalysisClass", originalD=f_tb, sampleda=sampleda, taxda=taxda, result=result, kwres=first.res,
                   secondvars=second.test.sig.vars, mlres=ml.res, someparams=params)
         return(res)
     }
     result <- .combine_others(result, first.res)
     if (action=="only"){
         return(result)
     }else if (action=="add"){
         newgroup <- paste0("Sign_", rlang::as_name(.group))
         result %<>% dplyr::rename(label="f", !!newgroup:=!!.group)
         if (is.null(taxatree)){
             if (!is.null(otutree(.data))){
                 otu.tree <- .data %>% mp_extract_otutree() %>%
                             treeio::full_join(result, by="label", suffix=c("", ".y"))
                 otutree(.data) <- otu.tree
             }else{             
                 otu_tb <- .data %>% 
                           mp_extract_feature() 
                 #result %<>% dplyr::select(c("label",setdiff(colnames(result), colnames(otu_tb))))
                 otu_tb %<>% 
                     dplyr::left_join(result %<>% dplyr::rename(OTU="label"), by="OTU", suffix=c("", ".y")) %>% 
                     tibble::column_to_rownames(var="OTU") %>%
                     S4Vectors::DataFrame(check.names=FALSE)
                 SummarizedExperiment::rowData(.data) <- otu_tb
             }
         }else{
             taxatree %<>% treeio::full_join(result, by="label", suffix=c("", ".y"))
             taxatree(.data) <- taxatree
         }
         return(.data)
     }
}

#' @rdname mp_diff_analysis-methods
#' @aliases mp_diff_analysis,MPSE
#' @exportMethod mp_diff_analysis
setMethod("mp_diff_analysis", signature(.data="MPSE"), .internal_mp_diff_analysis)

#' @rdname mp_diff_analysis-methods
#' @aliases mp_diff_analysis,tbl_mpse
#' @exportMethod mp_diff_analysis
setMethod("mp_diff_analysis", signature(.data="tbl_mpse"), .internal_mp_diff_analysis)

#' @rdname mp_diff_analysis-methods
#' @aliases mp_diff_analysis,grouped_df_mpse
#' @exportMethod mp_diff_analysis
setMethod("mp_diff_analysis", signature(.data="grouped_df_mpse"), .internal_mp_diff_analysis)


#' The visualization of result of mp_diff_analysis
#' @rdname mp_plot_diff_res-methods
#' @param .data MPSE or tbl_mpse after run mp_diff_analysis with \code{action="add"}
#' @param .group the column name for mapping the different color, default is the 
#' column name has 'Sign_' prefix, which contains the enriched group name, but the insignificant
#' should be NA.
#' @param layout the type of tree layout, should be one of "rectangular", "roundrect", "ellipse", 
#' "circular", "slanted", "radial", "inward_circular".
#' @param tree.type one of 'taxatree' and 'otutree', taxatree is the taxonomy class tree
#' 'otutree' is the phylogenetic tree built with the representative sequences.
#' @param barplot.x the column name of continuous value mapped to barplot, default is NULL, 
#' meaning the 'LDAmean' will be used internally.
#' @param point.size the column name of continuous value mapped to the size of point in the tree,
#' default is NULL, meaning the 'fdr' will be used internally.
#' @param sample.num integer when it is smaller than the sample number of '.data', the abundance 
#' of '.group' will replace the abundance of sample, default is 50.
#' @param .taxa.class character the name of taxonomy class level, default is NULL, meaning it will
#' extract the phylum annotation automatically.
#' @param tiplab.size numeric the size of tiplab, default is 2.
#' @param offset.abun numeric the gap (width) (relative width to tree) between the tree and abundance 
#' panel, default is 0.04.
#' @param pwidth.abun numeric the panel width (relative width to tree) of abundance panel, 
#' default is 0.3 .
#' @param offset.effsize numeric the gap (width) (relative width to tree) between the tree and 
#' effect size panel, default is 0.3 .
#' @param pwidth.effsize numeric the panel width (relative width to tree) of effect size panel, 
#' default is 0.5 .
#' @param group.abun logical whether to display the relative abundance of group instead of sample,
#' default is FALSE.
#' @param tiplab.linetype numeric the type of line for adding line if 'tree.type' is 'otutree',
#' default is 3 .
#' @param ... additional parameters, meaningless now.
#' @export
setGeneric("mp_plot_diff_res", 
                function(
                    .data,
                    .group, 
                    layout = "radial", 
                    tree.type = "taxatree", 
                    .taxa.class = NULL, 
                    barplot.x = NULL,
                    point.size = NULL,
                    sample.num = 50,
                    tiplab.size = 2,
                    offset.abun = 0.04,
                    pwidth.abun = 0.8,
                    offset.effsize = 0.3,
                    pwidth.effsize = 0.5,
                    group.abun = FALSE,
                    tiplab.linetype = 3,
                    ...) 
                    standardGeneric("mp_plot_diff_res")
)


#' @importFrom ggplot2 geom_col
#' @importFrom ggtreeExtra geom_fruit
.internal_mp_plot_diff_res <- function(.data,
                                       .group,
                                       layout = "radial",
                                       tree.type = "taxatree",
                                       .taxa.class = NULL,
                                       barplot.x = NULL,
                                       point.size = NULL,
                                       sample.num = 50,
                                       tiplab.size = 2,
                                       offset.abun = 0.04,
                                       pwidth.abun = 0.8,
                                       offset.effsize = 0.3,
                                       pwidth.effsize = 0.5,
                                       group.abun = FALSE,
                                       tiplab.linetype = 3,
                                       ...
                                      ){
    .taxa.class <- rlang::enquo(.taxa.class)
    .group <- rlang::enquo(.group)
    barplot.x <- rlang::enquo(barplot.x)
    point.size <- rlang::enquo(point.size)

    layout %<>% match.arg(c("rectangular", "roundrect", "ellipse", "circular", 
                            "slanted", "radial", "inward_circular"))

    tree.type %<>% match.arg(c("taxatree", "otutree"))

    if (tree.type == 'otutree'){
        anno.tree <- .data %>% mp_extract_tree(type="otutree") %>% suppressMessages()
        if (is.null(anno.tree)){
            stop_wrap("The otutree slot of the MPSE class is empty, you can try to 
                       select taxatree by setting tree.type to 'taxatree'.")
        }else{
            taxada <- .data %>% mp_extract_taxonomy() %>% suppressMessages()
            if (!is.null(taxada)){
                anno.tree %<>% dplyr::left_join(taxada, by=c("label"="OTU"))
                if (rlang::quo_is_null(.taxa.class)){
                    .taxa.class <- rlang::sym(colnames(taxada)[3])
                }
            }
        }
    }
    
    if (tree.type == "taxatree"){
        anno.tree <- .data %>% mp_extract_tree() %>% suppressMessages()
        if (is.null(anno.tree)){
            stop_wrap("The taxatree slot of the MPSE class is empty, you can try to 
                       select otutree by setting tree.type to 'otutree'.")
        }else{
            if (rlang::quo_is_null(.taxa.class)){
                .taxa.class <- anno.tree %>% 
                               tidytree::filter(!!rlang::sym("nodeDepth")==2, keep.td=FALSE) %>% 
                               pull(!!rlang::sym("nodeClass")) %>% 
                               unique()
                .taxa.class <- rlang::sym(.taxa.class)
            }
        }
    }

    tbl.f <- .data %>% mp_extract_feature()
    tbl.f %<>% dplyr::select(c(colnames(tbl.f)[1], 
                               setdiff(colnames(tbl.f)[-1], tidytree::get.fields(anno.tree)))
                            )

    anno.tree %<>% dplyr::left_join(tbl.f, by=c('label' = 'OTU'))

    nsample <- .data %>% mp_extract_sample() %>% nrow()
    field.da.nm <- tidytree::get.fields(anno.tree)

    if (any(grepl("LDAmean", field.da.nm)) && rlang::quo_is_null(barplot.x)){
        x.bar <- 'LDAmean'
        x.bar.title <- "log10(LDA)"
    }else if (any(grepl("MDAmean", field.da.nm)) && rlang::quo_is_null(barplot.x)){
        x.bar <- "MDAmean"
        x.bar.title <- "MDA"
    }else if (!rlang::quo_is_null(barplot.x)){
        x.bar <- x.bar.title <- gsub("^\"|\"$", "", rlang::as_label(barplot.x))
    }else{
        stop_wrap("Please provide barplot.x or verify the 'mp_diff_analysis' is done.")
    }
   
    if (rlang::quo_is_missing(.group)){
        if (any(grepl('^Sign_', field.da.nm))){
            sign.field <- field.da.nm[grepl('^Sign_', field.da.nm)][1]
            group.nm <- gsub('Sign_', "", sign.field)
        }else{
            stop_wrap('The .group name should be specified manually.')
        }
    }else{
        group.nm <- rlang::as_name(.group)
        if (!grepl('^Sign_', group.nm)){
            if (paste0('Sign_', group.nm) %in% field.da.nm){
                sign.field <- paste0('Sign_', group.nm)
            }else{
                #stop('Please check the mp_diff_analysis(..., action="add") has been run.')
                sign.field <- group.nm
            }
        }else{
            sign.field <- group.nm
            group.nm <- gsub('Sign_', "", group.nm)
        }
    }
    
    if (rlang::quo_is_null(point.size) && 'fdr' %in% field.da.nm){
        size.mapping <- '-log10(fdr)'
    }else if (!rlang::quo_is_null(point.size)){
        if (inherits(rlang::quo_get_expr(point.size), 'call')){
            size.mapping <- rlang::as_label(point.size)
        }else{
            size.mapping <- paste0("-log10(", gsub("^\"|\"$", "", rlang::as_label(point.size)), ")")
        }
    }else{
        stop_wrap("Please specify the 'point.size' argument manually !")
    }

    flag <- grepl(paste0("By", group.nm), field.da.nm)

    if (nsample > sample.num || group.abun){
        if (!any(flag)){
            stop_wrap("The relative abundance of each group will be displayed, but the
                       relative abundance of each group is not calculated, please run 
                       the mp_cal_abundance specified group argument before !")
        }
        abun.col <- field.da.nm[flag]
    }else{
        abun.col <- field.da.nm[grepl("BySample", field.da.nm)]
    }
    abun.col <- abun.col[1]
    x.abun.col <- anno.tree %>% 
                  dplyr::select(!!rlang::sym(abun.col)) %>%
                  tidyr::unnest(!!rlang::sym(abun.col)) %>%
                  colnames()
    if (any(grepl('BySample$',  abun.col))){
        if (any(grepl('^Rel', x.abun.col))){
            x.abun.col <- paste0('Rel', abun.col)
        }else{
            x.abun.col <- gsub('BySample$', '', abun.col)
        }
    }else{
        if (any(grepl("^Rel", x.abun.col))){
            x.abun.col <- paste0('Rel', abun.col)
        }else{
            x.abun.col <- abun.col
        }
    }
    #x.abun.col <- x.abun.col[grepl("^Rel", x.abun.col)]
    gplot.pck <- "ggplot2"
    require(gplot.pck, character.only=TRUE) %>% suppressMessages()
    if (tree.type == "otutree"){
        p1 <- ggtree(
                anno.tree,
                layout = layout,
                size = 0.3
              )
        if (!is.null(taxada)){
            p1 <- p1 +
                  geom_tiplab(
                    align = TRUE,
                    size = 0,
                    linetype = tiplab.linetype
                  ) +
                  geom_fruit(
                    data = td_filter(!!rlang::sym("isTip")),
                    geom = geom_point,
                    mapping = aes(colour = !!.taxa.class),
                    size = 1.5,
                    offset = 0
                  )    
        }
    }

    if (tree.type == "taxatree"){
        p1 <- suppressWarnings(
                ggtree(
                  anno.tree,
                  layout = layout,
                  size = 0.3
                  ) +
                geom_point(
                  data = td_filter(!.data$isTip),
                  fill = "white",
                  size = 1,
                  shape = 21
                )
              )
        
        p1 <- suppressWarnings(
                p1 +
                ggtree::geom_hilight(
                  data = td_filter(!!rlang::sym("nodeClass")==rlang::as_name(.taxa.class)),
                  mapping = aes(
                                node = !!rlang::sym("node"), 
                                fill = !!rlang::sym("label")
                            )
                )
              )
    }

    if (nsample > sample.num || group.abun){
         mapping <- aes(x=!!rlang::sym(group.nm), size = !!rlang::sym(x.abun.col), 
                        fill = !!rlang::sym(group.nm), subset = !!rlang::sym(x.abun.col) > 0) 
         n.pwidth <- .data %>%
                      mp_extract_sample() %>%
                      dplyr::pull(!!rlang::sym(group.nm)) %>%
                      unique() %>% length()    
    }else{
         mapping <- aes(x=forcats::fct_reorder(!!rlang::sym("Sample"), !!rlang::sym(group.nm), .fun=min),
                        size = !!rlang::sym(x.abun.col), 
                        fill = !!rlang::sym(group.nm), 
                        subset = !!rlang::sym(x.abun.col) > 0
                    )
         n.pwidth <- ncol(.data)
    }
    ggstar <- "ggstar"
    require(ggstar, character.only=TRUE) %>% suppressMessages()
    p2 <- suppressWarnings(
            p1 + 
            ggnewscale::new_scale_fill() +
            geom_fruit(
               data = td_unnest(!!rlang::sym(abun.col)),
               geom = geom_star,
               mapping = mapping,
               starshape = 13,
               starstroke = 0.05,
               offset = offset.abun,
               pwidth = pwidth.abun,
               grid.params = list(linetype=2)
            ) +  
            scale_size_continuous(
               name= ifelse(grepl('^Rel', abun.col), "Relative Abundance (%)", gsub("By.*", "", abun.col)),
               range = c(.5, 3),
               guide = guide_legend(override.aes = list(fill="black"))
            )
          )

    if (nsample > 50 || group.abun){
        p3 <- suppressWarnings(
                p2 + 
                geom_tiplab(size=tiplab.size, offset = max(p2$data$xmaxtmp, na.rm=TRUE) - 0.98*max(p2$data$x, na.rm=TRUE), 
                            align = TRUE, 
                            linetype=NA)
              )
    }else{
        p3 <- suppressWarnings(
                p2 + 
                geom_tiplab(size=tiplab.size, offset = max(p2$data$xmaxtmp, na.rm=TRUE) - 0.95*max(p2$data$x, na.rm=TRUE))
              )
    }
    # display the LDA of significant OTU.
    title.height <- 4.4e-06 * sum(p3$data$isTip) 
    p4 <- suppressWarnings(
            p3 +
            ggnewscale::new_scale_fill() +
            geom_fruit(
               #data = td_filter(!is.na(!!rlang::sym(x.bar))),
               data = td_filter(!is.na(!!rlang::sym(sign.field))),
               geom = "geom_col",
               mapping = aes(
                             x = !!rlang::sym(x.bar),
                             fill = !!rlang::sym(sign.field)
                             ),
               orientation = "y",
               offset = offset.effsize,
               pwidth = pwidth.effsize,
               axis.params = list(axis = "x",
                                  title = x.bar.title,
                                  title.height = title.height,
                                  title.size = 2,
                                  text.size = 1.8,
                                  vjust = 1),
               grid.params = list(linetype = 2)
              )
          )

    # display the significant (FDR) taxonomy after kruskal.test (default)
    p5 <- suppressWarnings(
            p4 +
            ggnewscale::new_scale("size") +
            geom_point(
               data=td_filter(!is.na(!!rlang::sym(sign.field))),
               mapping = aes_string(
                              size = size.mapping,
                              fill = sign.field,
                             ),
               shape = 21
            ) +
            scale_size_continuous(range=c(1, 3)) #+
            #scale_fill_manual(values=c("#1B9E77", "#D95F02"))
          )
    
    p6 <- suppressWarnings(
            p5 + theme(
               legend.key.height = unit(0.3, "cm"),
               legend.key.width = unit(0.3, "cm"),
               legend.spacing.y = unit(0.02, "cm"),
               legend.text = element_text(size = 7),
               legend.title = element_text(size = 9),
              )
          )

    return (suppressWarnings(p6))
}

#' @rdname mp_plot_diff_res-methods
#' @aliases mp_plot_diff_res,MPSE
#' @export mp_plot_diff_res
setMethod("mp_plot_diff_res", signature(.data='MPSE'), .internal_mp_plot_diff_res)

#' @rdname mp_plot_diff_res-methods
#' @aliases mp_plot_diff_res,tbl_mpse
#' @export mp_plot_diff_res
setMethod("mp_plot_diff_res", signature(.data="tbl_mpse"), .internal_mp_plot_diff_res)

#' @rdname mp_plot_diff_res-methods
#' @aliases mp_plot_diff_res,grouped_df_mpse
#' @export mp_plot_diff_res
setMethod("mp_plot_diff_res", signature(.data="grouped_df_mpse"), .internal_mp_plot_diff_res)


#' Visualizing the result of mp_diff_analysis with cladogram.
#' @param .data MPSE object or treedata which was from the taxatree slot 
#' after running the 'mp_diff_analysis'.
#' @param .group the column name for mapping the different color.
#' @param .size the column name for mapping the size of points, default is 'pvalue'.
#' @param taxa.class the taxonomy class name will be replaced shorthand, default is 
#' the one level above ‘OTU’.
#' @param removeUnknown logical, whether mask the unknown taxonomy information but differential species,
#' default is FALSE.
#' @param layout character, the layout of tree, default is 'radial', see also the 'layout' of 'ggtree'.
#' @param hilight.alpha numeric, the transparency of high light clade, default is 0.3.
#' @param hilight.size numeric, the margin thickness of high light clade, default is 0.2.
#' @param bg.tree.size numeric, the line size (width) of tree, default is 0.15.
#' @param bg.tree.color character, the line color of tree, default is '#bed0d1'.
#' @param bg.point.color character, the color of margin of background node points of tree, default is '#bed0d1'.
#' @param bg.point.fill character, the point fill (since point shape is 21) of background nodes of 
#' tree, default is 'white'.
#' @param bg.point.stroke numeric, the margin thickness of point of background nodes of tree,
#' default is 0.2 .
#' @param bg.point.size numeric, the point size of background nodes of tree, default is 2.
#' @param label.size numeric, the label size of differential taxa, default is 2.6.
#' @param tip.annot logcial whether to replace the differential tip labels with shorthand,
#' default is TRUE.
#' @param as.tiplab logical, whether to display the differential tip labels with 'geom_tiplab' 
#' of 'ggtree', default is TRUE, if it is FALSE, it will use 'geom_text_repel' of 'ggrepel'.
#' @param ... additional parameters, meaningless now.
#' @export
#' @details
#' The color scale of differential group can be designed by 'scale_fill_diff_cladogram'
#' @importFrom ggtree td_mutate td_filter
#' @examples
#' \dontrun{
#'   data(mouse.time.mpse)
#'   mouse.time.mpse %<>%
#'     mp_rrarefy()
#'   mouse.time.mpse
#'   mouse.time.mpse %<>%
#'     mp_diff_analysis(.abundance=RareAbundance,
#'                      .group=time,
#'                      first.test.alpha=0.01,
#'                      action="add")
#'   #' ### visualizing the differential taxa with cladogram
#'   library(ggplot2)
#'   f <- mouse.time.mpse %>%
#'        mp_plot_diff_cladogram(
#'          label.size = 2.5,
#'          hilight.alpha = .3,
#'          bg.tree.size = .5,
#'          bg.point.size = 2,
#'          bg.point.stroke = .25
#'        ) +
#'        scale_fill_diff_cladogram(
#'          values = c('skyblue', 'orange')
#'        ) +
#'        scale_size_continuous(range = c(1, 4))
#'   f
#' }
mp_plot_diff_cladogram <- function(
    .data, 
    .group, 
    .size = 'pvalue',
    taxa.class, 
    removeUnknown = FALSE,
    layout = 'radial', 
    hilight.alpha = .3,
    hilight.size = .2,
    bg.tree.size = .15, 
    bg.tree.color = '#bed0d1',
    bg.point.color = '#bed0d1',
    bg.point.fill = 'white',
    bg.point.stroke = .2,
    bg.point.size = 2, 
    label.size = 2.6,
    tip.annot = TRUE, 
    as.tiplab = TRUE,
    ...){
        .group <- rlang::enquo(.group)
        taxa.class <- rlang::enquo(taxa.class)
        .size <- rlang::enquo(.size)
        if (!(inherits(.data, 'MPSE') || inherits(.data, 'treedata') || inherits(.data, 'tbl_mpse') || inherits(.data, 'grouped_df_mpse'))){
            stop("The .data should be an MPSE or treedata object after running mp_diff_analysis ")
        }
        if (inherits(.data, 'MPSE') || inherits(.data, 'tbl_mpse') || inherits(.data, 'grouped_df_mpse')){
            tbl <- .data %>% mp_extract_feature()
            .data %<>% mp_extract_taxatree()
            if (is.null(.data)){
                stop_wrap('The .data object does not have the taxatree slot.')
            }else{
                tbl %<>% dplyr::select(c(colnames(tbl)[1], 
                                         setdiff(colnames(tbl)[-1], tidytree::get.fields(.data))))
                .data %<>% tidytree::left_join(tbl, by=c('label' = 'OTU'))
            }
        }
        .data %<>% as_tibble()
        nmda <- colnames(.data)
        if (rlang::quo_is_missing(.group)){
            if (any(grepl('^Sign_', nmda))){
                .group <- rlang::sym(nmda[grepl('^Sign_', nmda)][1])
            }else{
                stop_wrap('The .group name should be specified manually.')
            }
        }else{
            .group <- rlang::as_name(.group)
            if (!grepl('^Sign_', .group)){
                if (paste0('Sign_', .group) %in% nmda){
                    .group <- rlang::sym(paste0('Sign_', .group))
                }else{
                    .group <- rlang::sym(.group)
                }
            }else{
                .group <- rlang::sym(.group)
            }
        }
        if (!rlang::quo_is_missing(taxa.class)){
            taxa.class <- rlang::as_name(taxa.class)
        }else{
            if (any(is.na(.data$nodeDepth))){
                dat <- .data %>% dplyr::filter(!is.na(.data$nodeDepth)) 
            }else{
                dat <- .data
            }
            if (length(unique(dat$nodeClass))==1){
                taxa.class <- unique(dat$nodeClass)
            }else{
                taxa.class <- dat %>% 
	               dplyr::filter(.data$nodeDepth == (max(.data$nodeDepth) - 2)) %>%
                   dplyr::pull(.data$nodeClass) %>% unique()

            }
        }
        annot.index <- .get_annot_index(x = .data, taxa.class = taxa.class, tip.annot = tip.annot)
	    
        .data <- .generate_annot_df(x = .data, annot.index, .group, removeUnknown) %>% 
                  as.treedata()
        offset.max <- max(dplyr::pull(.data, .data$nodeDepth), na.rm = TRUE) 
        
        if (inherits(rlang::quo_get_expr(.size), "call")){
            mapping.point <- aes_string(size = rlang::as_label(.size), 
                                        fill = rlang::as_name(.group))
        }else{
            mapping.point <- aes_string(size = paste0("-log10(", gsub("^\"|\"$", "", rlang::as_name(.size)), ")"), 
                                        fill = rlang::as_name(.group))
        }

        p <- ggtree(.data, layout = layout, size = bg.tree.size, color = bg.tree.color) +
             geom_point(
                data = ifelse(removeUnknown, td_filter(is.na(!!.group) | grepl('__un_', .data$label)),
                         td_filter(is.na(!!.group))
                       ),
                mapping = aes_(x = ~x, y = ~y), 
                color = bg.point.color,
                fill = bg.point.fill, 
                size = bg.point.size, 
                shape = 21, 
                stroke = bg.point.stroke,
             ) +
             ggtree::geom_hilight(
                data = td_mutate(
                         extend = (offset.max + 1 - .data$nodeDepth) * 1,
                         .f = ifelse(removeUnknown, td_filter(!is.na(!!.group) & !grepl('__un_', .data$label)), 
                                     td_filter(!is.na(!!.group)))
                       ),
                mapping = aes_string(node = "node", fill = rlang::as_name(.group), color = rlang::as_name(.group), extend = "extend"),
                size = hilight.size,
                alpha = hilight.alpha,
                show.legend = FALSE,
             )
        
        if (removeUnknown){
            flag.clade <- td_filter(!is.na(!!.group) & !.data$isTip & !grepl('__un_', .data$label))(p$data)
        }else{
            flag.clade <- td_filter(!is.na(!!.group) & !.data$isTip)(p$data)
        }
        
        if (nrow(flag.clade) > 0){
            p <- p +
             ggtree::geom_cladelab(
                data = td_mutate(
                         offset = (offset.max + 1 - .data$nodeDepth) * 0.88,
                         .f = ifelse(removeUnknown, td_filter(!is.na(!!.group) & !.data$isTip & !grepl('__un_', .data$label)),
									 td_filter(!is.na(!!.group) & !.data$isTip))
                       ),
                mapping = aes_string(node = "node", label = "label", offset = "offset"),
                geom = 'shadowtext',
                angle = "auto",
                horizontal = FALSE,
                hjust = 0.5,
                fontsize = label.size,
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
                         .f = td_filter(!is.na(!!.group) & .data$nodeDepth %in% annot.index),
                ),
                mapping = aes_(x = 0, y = 0, color = ~.LABEL),
                size = 0,
                stroke = 0
             ) +
             geom_point(
               data = ifelse(removeUnknown, td_filter(!grepl("__un_", .data$label) & !is.na(!!.group)), td_filter(!is.na(!!.group))),
               mapping = mapping.point,
               shape = 21,
               stroke = .1,
               #show.legend = c(fill = FALSE)
             )
         
        if (!as.tiplab){
            p <- p + 
                 ggrepel::geom_text_repel(
                    data = ifelse(removeUnknown, td_filter(!is.na(!!.group) & .data$isTip & !grepl("__un_", .data$label)), 
                             td_filter(!is.na(!!.group) & .data$isTip)),
                    mapping = aes_(x = ~x, y = ~y, label = ~label),
                    size = label.size,
                    min.segment.length = 0,
                    bg.colour = 'white'
                 )
        }else{
	        p <- p +
                 ggtree::geom_tiplab(
                    data = ifelse(removeUnknown, td_filter(!is.na(!!.group) & !grepl("__un_", .data$label)),
                             td_filter(!is.na(!!.group))),
                    offset = 0.78,
                    size = label.size,
                    geom = 'shadowtext',
					color = 'black',
                    bg.colour = 'white'			
                 )
	    
	    }
        p <- p + theme(
          legend.key.width = unit(.3, 'cm'),
          legend.key.height = unit(.3, 'cm'),
          legend.text = element_text(size=6),
          legend.title = element_text(size=8),
          legend.margin = ggplot2::margin(-.25, 0, 0, 0, 'cm')
        )
        p <- p + scale_fill_diff_cladogram()
        return (p)
}

#' displaying the differential result contained abundance and LDA with 
#' boxplot (abundance) and error bar (LDA).
#' @rdname mp_plot_diff_boxplot-methods
#' @param .data MPSE or tbl_mpse after run mp_diff_analysis with 'action="add"'.
#' @param .group the column name for mapping the different color.
#' @param .size the column name for mapping the size of points or numeric, default is 2.
#' @param errorbar.xmin the column name for 'xmin' mapping of error barplot layer, default is NULL.
#' @param errorbar.xmax the column name for 'xmax' mapping of error barplot layer, default is NULL.
#' @param point.x the column name for 'x' mapping of point layer (right panel), default is NULL.
#' @param taxa.class the taxonomy class features will be displayed, default is 'all'.
#' @param group.abun logical whether plot the abundance in each group with bar plot,
#' default is FALSE.
#' @param removeUnknown logical whether mask the unknown taxonomy information but 
#' differential species, default is FALSE.
#' @param ... additional params, see also the 'geom_boxplot', 'geom_errorbarh' and 'geom_point'.
#' @importFrom ggfun theme_blinds
#' @export
#' @examples
#' data(mouse.time.mpse)
#' mouse.time.mpse %<>%
#'   mp_rrarefy()
#' mouse.time.mpse
#' mouse.time.mpse %<>%
#'   mp_diff_analysis(.abundance=RareAbundance,
#'                    .group=time,
#'                    first.test.alpha=0.01,
#'                    action="add")
#' library(ggplot2)
#' p1 <- mouse.time.mpse %>% 
#'         mp_plot_diff_boxplot(.group = time) %>%
#'         set_diff_boxplot_color(
#'           values = c("deepskyblue", "orange"),
#'           guide = guide_legend(title=NULL)
#'         )
#' p1
#' p2 <- mouse.time.mpse %>% 
#'         mp_plot_diff_boxplot(
#'           taxa.class = c(Genus, OTU),
#'           group.abun = TRUE, 
#'           removeUnknown = TRUE,
#'         ) %>%
#'         set_diff_boxplot_color(
#'           values = c("deepskyblue", "orange"),
#'           guide = guide_legend(title=NULL)
#'         )
#' p2 
setGeneric("mp_plot_diff_boxplot",
  function(
    .data,
    .group,
    .size = 2,
    errorbar.xmin = NULL,
    errorbar.xmax = NULL,
    point.x = NULL,
    taxa.class = 'all',
    group.abun = FALSE,
    removeUnknown = FALSE,
    ...
  )
  standardGeneric('mp_plot_diff_boxplot')
)

.internal_mp_plot_diff_boxplot <- function(.data, .group, .size = 2, 
                                           errorbar.xmin = NULL, 
                                           errorbar.xmax = NULL,
                                           point.x = NULL, 
                                           taxa.class = 'all', group.abun = FALSE, removeUnknown=FALSE, ...){
    taxa.class <- rlang::enquo(taxa.class)
    .group <- rlang::enquo(.group)
    .size <- rlang::enquo(.size)
    errorbar.xmin <- rlang::enquo(errorbar.xmin)
    errorbar.xmax <- rlang::enquo(errorbar.xmax)
    point.x <- rlang::enquo(point.x)

    params <- list(...)
    taxa.class <- quo.vect_to_str.vect(taxa.class)
    if (!is.null(suppressMessages(taxatree(.data)))){
        tbl <- .data %>% mp_extract_abundance(taxa.class = 'all')
        xx <- .data %>% mp_extract_feature()
        xx %<>% dplyr::select(setdiff(colnames(xx), colnames(tbl)))
        tbl %<>% left_join(xx, by = c('label'='OTU'))
    }else if (is.null(suppressMessages(taxatree(.data))) || taxa.class == 'OTU'){
        tbl <- .data %>% mp_extract_feature()
        if (!any(grepl("AbundanceBySample$", colnames(tbl)))){
            trash <- try(silent = TRUE,  
                         expr = {
                              .data <- suppressMessages(mp_cal_abundance(.data, .abundance = "Abundance"))  
                           }
                        )
            if (inherits(trash, "try-error")) {
                .data <- mp_cal_abundance(.data, .abundance = "Abundance", force = TRUE)
            }
            tbl <- .data %>% mp_extract_feature()
        }
        tbl %<>% dplyr::rename(label = 'OTU') 
    }
    if (!any(taxa.class %in% c('all', 'ALL', 'All')) && !is.null(suppressMessages(taxatree(.data)))){
       tbl %<>% dplyr::filter(.data$nodeClass %in% taxa.class)
    }
    if (removeUnknown){
        tbl %<>% dplyr::filter(!grepl('__un_', .data$label))
    }
    nmda <- colnames(tbl)
    nm.abun <- nmda[grepl('BySample', nmda)][1]
    tbl %<>% tidyr::unnest(nm.abun)
    nmda <- colnames(tbl)
    if (rlang::quo_is_missing(.group)){
        if (any(grepl('^Sign_', nmda))){
            sign.group <- nmda[grepl('^Sign_', nmda)][1]
            .group <- gsub('Sign_', "", sign.group)
        }else{
            stop_wrap('The .group name should be specified manually.')
        }
    }else{
        .group <- rlang::as_name(.group)
        if (!grepl('^Sign_', .group)){
            if (paste0('Sign_', .group) %in% nmda){
                sign.group <- paste0('Sign_', .group)
            }else{
                #stop('Please check the mp_diff_analysis(..., action="add") has been run.')
                sign.group <- .group
            }
        }else{
            sign.group <- .group
            .group <- gsub('Sign_', "", .group)
        }
    }
    tbl %<>% dplyr::filter(!is.na(!!rlang::sym(sign.group)))
    if (any(grepl('LDA', nmda)) && rlang::quo_is_null(errorbar.xmin) && 
        rlang::quo_is_null(errorbar.xmax) && rlang::quo_is_null(point.x)){
        xlabtext <- bquote(paste(Log[10],"(",.("LDA"), ")"))
        xtext <- "LDAmean"
        xmintext <- "LDAlower"
        xmaxtext <- "LDAupper"
    }else if ('MDA' %in% nmda && rlang::quo_is_null(errorbar.xmin) && 
              rlang::quo_is_null(errorbar.xmax) && rlang::quo_is_null(point.x)){
        xlabtext <- "MDA"
        xtext <- "MDAmean"
        xmintext <- "MDAlower"
        xmaxtext <- "MDAupper"
    }else if (!rlang::quo_is_null(errorbar.xmin) && !rlang::quo_is_null(errorbar.xmax) && !rlang::quo_is_null(point.x)){
        xlabtext <- gsub("^\"|\"$", "", rlang::as_label(point.x))
        xtext <- xlabtext
        xmintext <- gsub("^\"|\"$", "", rlang::as_label(errorbar.xmin))
        xmaxtext <- gsub("^\"|\"$", "", rlang::as_label(errorbar.xmax))

    }else if (!rlang::quo_is_null(point.x) && any(rlang::quo_is_null(errorbar.xmin) || rlang::quo_is_null(errorbar.xmax))){
        xlabtext <- gsub("^\"|\"$", "", rlang::as_label(point.x))
        xtext <- xlabtext
        xmintext <- NULL
        xmaxtext <- NULL
    }
    if (any(grepl('Rel.*BySample', nmda))){
        abunda <- nmda[grepl('Rel.*BySample', nmda)]
    }else{
        abunda <- gsub('BySample$', '', nm.abun) 
    }
    if (group.abun){
        mapping1 <- aes(x = !!rlang::sym(abunda), 
                        y = !!rlang::sym("label"),
                        group = !!rlang::sym(.group),
                        fill = !!rlang::sym(.group)
        )
        panel1.geom <- 'geom_bar'
        panel1.args <- .extract_args(geom = 'geom_bar')
        panel1.args <- params[names(params) %in% panel1.args]
        panel1.args <- panel1.args[!names(panel1.args) %in% c('fill', 'group')]
        panel1.args$fun <- mean
        panel1.args$stat <- "summary"
        panel1.args$position <- "dodge"
        panel1.args$orientation <- 'y'
    }else{
        mapping1 <- aes(
                  x = !!rlang::sym(abunda),
                  y = !!rlang::sym("label"),
                  fill = !!rlang::sym(.group)
                )
        panel1.geom <- 'geom_boxplot'
        panel1.args <- .extract_args(geom="geom_boxplot")
        panel1.args <- params[names(params) %in% panel1.args]
        panel1.args <- panel1.args <- panel1.args[!names(panel1.args) %in% c('fill')]
    }
    
    if (!is.null(xmintext)){
        mapping2 <- aes(
                  xmin = !!rlang::sym(xmintext),
                  xmax = !!rlang::sym(xmaxtext),
                  y = !!rlang::sym('label')
        )
        panel2.geom <- 'geom_errorbarh'
    }else{
        mapping2 <- aes(x = 0, 
                        xend = !!rlang::sym(xtext), 
                        y = !!rlang::sym('label'),
                        yend = !!rlang::sym('label')
        )
        panel2.geom <- 'geom_segment'
    }
    
    panel2.args <- .extract_args(geom = panel2.geom)
    panel2.args <- params[names(params) %in% panel2.args]
    if (panel2.geom == 'geom_errorbarh'){
        panel2.args$height <- .3
        if ('height' %in% names(params)){
            panel2.args$height <- params$height        
        }
    }

    mapping3 <- aes_string(
                    x = xtext, 
                    y = 'label', 
                    color = sign.group#, 
                    #size = paste0("-log10(", rlang::as_name(.size), ")")
                )

    point.args <- .extract_args(geom='geom_point')
    point.args <- params[names(params) %in% point.args]
    
    point.args <- point.args[!names(point.args) %in% c("color", 'colour', "shape")]
    
    panel1.args$mapping <- mapping1
    panel2.args$mapping <- mapping2
    point.args$mapping <- mapping3
    point.args$show.legend <- c(colour=FALSE)

    if (is.quo.numeric(.size)){
        point.args <- point.args[!names(point.args) %in% c("size")]
        point.args$size <- as.numeric(rlang::as_label(.size))
    }else{
        if (inherits(rlang::quo_get_expr(.size), 'call')){
            size.mapping <- aes(size=!!.size)
        }else{
            if(any(grepl("\\(*.\\)$", rlang::quo_get_expr(.size)))){
                warning_wrap('The .y argument might be a call, you should remove the \' or \\" quote symbols.')
                size.mapping <- aes_string(size=gsub("^\"|\"$", "", rlang::as_label(.size)))
            }else{
                size.mapping <- aes_string(size=paste0("-log10(", gsub("^\"|\"$", "", rlang::as_label(.size)), ")"))
            }
        }
        point.args$mapping <- modifyList(point.args$mapping, size.mapping)
    }
    
    p <- ggplot(tbl)

    panel1.geom <- do.call(panel1.geom, panel1.args)
    p1 <- p + 
          panel1.geom + 
          ylab(NULL) + 
          xlab('Abundance') +
          scale_x_continuous(expand = c(0, 0, 0, .2))
    p1 <- p1 + theme_blinds(colour = c("grey90", "white"), axis.ticks.y=element_blank(), 
                            panel.grid.major.x = element_blank(), 
                            panel.grid.minor.x = element_blank()
               )

    panel2.geom <- do.call(panel2.geom, panel2.args)
    point.geom <- do.call(geom_point, point.args)
    p2 <- p + 
          panel2.geom +
          point.geom +
          ylab(NULL) +
          xlab(xlabtext) +
          theme(
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank()
          )
    p2 <- p2 + theme_blinds(colour = c("grey90", "white"))
    fix.group <- levels(tbl[[.group]])
    if (is.null(fix.group)){
        fix.group <- unique(sort(tbl[[.group]]))
    }
    #fix.group.color <- p1$scales$get_scales('fill')$palette.cache
    fix.group.color <- scales::hue_pal()(length(fix.group))
    names(fix.group.color) <- fix.group
    fix.group.color <- fix.group.color[match(unique(as.character(tbl[[sign.group]])), names(fix.group.color))]
    p2 <- p2 + scale_color_manual(values=fix.group.color)

    p <- p1 %>% aplot::insert_right(p2, width = .8)
    return(p)
}

.extract_args <- function(geom){
    xx <- do.call(geom, list())
    c(names(xx$geom_params), names(xx$geom$default_aes))
}

#' @rdname mp_plot_diff_boxplot-methods
#' @aliases mp_plot_diff_boxplot,MPSE
#' @export mp_plot_diff_boxplot
setMethod("mp_plot_diff_boxplot", signature(.data='MPSE'), .internal_mp_plot_diff_boxplot)

#' @rdname mp_plot_diff_boxplot-methods
#' @aliases mp_plot_diff_boxplot,tbl_mpse
#' @export mp_plot_diff_boxplot
setMethod("mp_plot_diff_boxplot", signature(.data="tbl_mpse"), .internal_mp_plot_diff_boxplot)

#' @rdname mp_plot_diff_boxplot-methods
#' @aliases mp_plot_diff_boxplot,grouped_df_mpse
#' @export mp_plot_diff_boxplot
setMethod("mp_plot_diff_boxplot", signature(.data="grouped_df_mpse"), .internal_mp_plot_diff_boxplot)

#' displaying the differential result contained abundance and LDA with
#' manhattan plot.
#' @rdname mp_plot_diff_manhattan-methods
#' @param .data MPSE or tbl_mpse after run 'mp_diff_analysis' with 'action="add"'.
#' @param .group the column name for mapping the different color.
#' @param .y the column name for mapping the y axis, default is 'fdr'.
#' @param .size the column name for mapping the size of points or numeric, default is 2.
#' @param taxa.class the taxonomy class features will be displayed, default is 'OTU'.
#' @param anno.taxa.class the taxonomy class to annotate the sign taxa with color, default is
#' 'Phylum' if 'taxatree' is not empty.
#' @param removeUnknown logical whether mask the unknown taxonomy information but
#' differential species, default is FALSE.
#' @param ... additional params, see also the 'geom_text_repel' and 'geom_point'.
#' @export
#' @examples
#' data(mouse.time.mpse)
#' mouse.time.mpse %<>%
#'   mp_rrarefy()
#' mouse.time.mpse
#' mouse.time.mpse %<>%
#'   mp_diff_analysis(.abundance=RareAbundance,
#'                    .group=time,
#'                    first.test.alpha=0.01,
#'                    action="add")
#' p <- mouse.time.mpse %>% 
#'        mp_plot_diff_manhattan(
#'            .group = Sign_time, 
#'            .y = fdr,
#'            .size = 2,
#'            taxa.class = OTU,
#'            anno.taxa.class = Phylum,
#'        )
setGeneric('mp_plot_diff_manhattan', 
   function(
     .data, 
     .group,
     .y='fdr',
     .size=2,
     taxa.class='OTU', 
     anno.taxa.class=NULL, 
     removeUnknown=FALSE, ...)
       standardGeneric('mp_plot_diff_manhattan')
)

#' @importFrom ggplot2 guide_axis labs element_line geom_hline coord_cartesian scale_x_discrete
.internal_mp_plot_diff_manhattan <- function(.data,
                                             .group,
                                             .y = 'fdr',
                                             .size = 2,
                                             taxa.class = 'OTU',
                                             anno.taxa.class = NULL,
                                             removeUnknown=FALSE, ...){
    taxa.class <- rlang::enquo(taxa.class)
    anno.taxa.class <- rlang::enquo(anno.taxa.class)
    .group <- rlang::enquo(.group)
    .size <- rlang::enquo(.size)
    .y <- rlang::enquo(.y)

    params <- list(...)
    taxa.class <- quo.vect_to_str.vect(taxa.class)
    taxa.tree <- suppressMessages(taxatree(.data))
    if (!is.null(taxa.tree)){
        tbl <- taxa.tree %>% tibble::as_tibble() %>% dplyr::filter(.data$label!='r__root')
        if (rlang::quo_is_null(anno.taxa.class)){
            anno.taxa.class <- tbl %>% dplyr::filter(.data$nodeDepth==2) %>%
                dplyr::pull(.data$nodeClass) %>% unique()
        }else{
            anno.taxa.class <- rlang::as_name(anno.taxa.class)
        }
        depth.value <- tbl %>% dplyr::filter(.data$nodeClass == anno.taxa.class) %>%
            dplyr::pull(.data$nodeDepth) %>% unique()
        tbl %<>% dplyr::filter(.data$nodeDepth >= depth.value)
        if (!any(taxa.class %in% c('all', 'ALL', 'All'))){
            tbl <- tbl %>% dplyr::filter(.data$nodeClass %in% taxa.class)
        }else{
            taxa.class <- tbl %>% dplyr::pull(.data$nodeClass) %>% unique()
        }
        taxa.da <- .data %>% mp_extract_taxonomy()
        taxa.da <- lapply(taxa.class, function(i){
            df <- taxa.da[,colnames(taxa.da) %in% c(i, anno.taxa.class),drop=FALSE]
            if (i == anno.taxa.class){
                df <- df %>% dplyr::mutate(.ANNO.TAXA=!!rlang::sym(anno.taxa.class))
            }else{
                df <- df %>% dplyr::rename(.ANNO.TAXA=!!rlang::sym(anno.taxa.class))
            }
            df %<>% dplyr::rename(label=!!rlang::sym(i)) %>% dplyr::distinct()
            return(df)
        }) %>% dplyr::bind_rows()
        tbl %<>% dplyr::left_join(taxa.da, by='label')
        if (removeUnknown){
            tbl %<>% dplyr::filter(!grepl('__un_', .data$label))
        }
    }else{
        tbl <- .data %>% mp_extract_feature() %>%
                dplyr::rename(label = 'OTU')
        taxa.da <- NULL
    }
    nmda <- colnames(tbl)
    if (rlang::quo_is_missing(.group)){
         if (any(grepl('^Sign_', nmda))){
             sign.group <- nmda[grepl('^Sign_', nmda)][1]
             .group <- gsub('Sign_', "", sign.group)
         }else{
             stop_wrap('The .group name should be specified manually.')
         }
     }else{
         .group <- rlang::as_name(.group)
         if (!grepl('^Sign_', .group)){
             if (paste0('Sign_', .group) %in% nmda){
                 sign.group <- paste0('Sign_', .group)
             }else{
                 #stop('Please check the mp_diff_analysis(..., action="add") has been run.')
                 sign.group <- .group
             }
         }else{
             sign.group <- .group
             .group <- gsub('Sign_', "", .group)
         }
     }
     tbl %<>% dplyr::filter(!is.na(!!.y))
     tbl[[sign.group]][is.na(tbl[[sign.group]])] <- 'NoSign'
     if (!is.null(taxa.da)){
         tbl %<>% dplyr::group_split(.data$.ANNO.TAXA) %>% dplyr::bind_rows()
         label.order <- tbl %>% dplyr::pull("label")
         tbl %<>% dplyr::mutate(label=factor(.data$label, levels=label.order))
         taxa.da <- tbl %>% dplyr::group_by(.data$.ANNO.TAXA) %>% dplyr::slice(ceiling(dplyr::n()/2)) %>%
                    dplyr::mutate_at('label', as.character)
     }
     mapping1 <- aes(x=!!rlang::sym('label'), shape=!!rlang::sym(sign.group)) 
     if (inherits(rlang::quo_get_expr(.y), 'call')){
         y.aes <- aes(y=!!.y)
     }else{
         if(any(grepl("\\(*.\\)$", rlang::quo_get_expr(.y)))){
             warning_wrap('The .y argument might be a call, you should remove the \' or \\" quote symbols.')
             y.aes <- aes_string(y=gsub("^\"|\"$", "", rlang::as_label(.y)))
         }else{
             y.aes <- aes_string(y=paste0("-log10(", gsub("^\"|\"$", "", rlang::as_label(.y)), ")"))
         }
     }
     mapping1 <- modifyList(mapping1, y.aes)
     if (!is.null(taxa.da)){
         mapping1 <- modifyList(mapping1, aes_string(color=".ANNO.TAXA"))
     }

     p <- ggplot(data=tbl)
     point.args <- params[names(params) %in% setdiff(.extract_args('geom_point'), c('size', 'color'))]
     point.args$mapping <- mapping1
     if (is.quo.numeric(.size)){
         point.args$size <- as.numeric(rlang::as_label(.size))
     }else{
         if (inherits(rlang::quo_get_expr(.size), 'call')){
             size.mapping <- aes(size=!!.size)
         }else{
             if(any(grepl("\\(*.\\)$", rlang::quo_get_expr(.size)))){
                 warning_wrap('The .y argument might be a call, you should remove the \' or \\" quote symbols.')
                 #size.mapping <- aes_string(size=gsub("^\"|\"$", "", rlang::as_label(.size)))
             }
             #size.mapping <- aes_string(size=paste0("-log10(", gsub("^\"|\"$", "", rlang::as_label(.size)), ")"))
             size.mapping <- aes_string(size=gsub("^\"|\"$", "", rlang::as_label(.size)))
         }
         point.args$mapping <- modifyList(point.args$mapping, size.mapping)
     }
     point.layer <- do.call('geom_point', point.args)
     if (!is.null(taxa.da)){
         axis.x.breaks <- taxa.da %>% dplyr::pull('label')
         axis.x.labels <- taxa.da %>% dplyr::pull('.ANNO.TAXA')
         rect.args <- params[names(params) %in% setdiff(.extract_args('geom_rect'), c('fill', 'color', 'alpha'))]
         rect.args$data <- taxa.da
         rect.args$mapping <- aes_string(fill=".ANNO.TAXA")
         rect.args$xmin <- -Inf
         rect.args$xmax <- Inf
         rect.args$ymin <- -Inf
         rect.args$ymax <- Inf
         rect.args$color <- NA
         rect.args$alpha <- .3
         rect.layer <- do.call('geom_rect', rect.args)
         p <- p +
              rect.layer +
              scale_fill_manual(values=rep(c('white', 'grey50'), ceiling(nrow(taxa.da)/2)), guide='none')
     }
     if (utils::packageVersion("ggplot2") >= '3.4.0'){
         p <- p + geom_hline(yintercept=-log10(0.05), color='red', linewidth = .2, linetype=2)
         line.theme <- element_line(linewidth = .2)
     }else{
         p <- p + geom_hline(yintercept=-log10(0.05), color='red', size = .2, linetype=2)
         line.theme <- element_line(size = .2)
     }
     p <- p + point.layer
     text.args <- list(min.segment.length = 0,
                       max.overlaps = 100,
                       color = 'black',
                       bg.color = 'white',
                       size = 2.5)
     text.args <- modifyList(text.args, params[names(params) %in% .extract_args('geom_text_repel')])
     text.args$data <- ggtree::td_filter(!!rlang::sym(sign.group)!='NoSign')
     text.args$mapping <- modifyList(mapping1, aes_string(label='label', shape=NULL))
     text.layer <- do.call('geom_text_repel', text.args)
     p <- p + text.layer

     if (!is.null(taxa.da)){
         axis.x.breaks <- taxa.da %>% dplyr::pull('label')
         axis.x.labels <- taxa.da %>% dplyr::pull('.ANNO.TAXA')
         p <- p +
              labs(color=anno.taxa.class) +
              facet_grid(. ~ .ANNO.TAXA, scales='free_x', space='free_x') +
              scale_x_discrete(breaks=axis.x.breaks, labels=axis.x.labels, guide=guide_axis(angle=-45))
     }
     p <- p +
          labs(x=NULL) +
          theme(
            strip.text.x = element_blank(),
            strip.background = element_blank(),
            panel.grid = element_blank(),
            panel.spacing.x = unit(.15, 'cm'),
            panel.background = element_blank(),
            axis.line = line.theme
          ) +
          coord_cartesian(clip = 'off')
     return(p)
}     
     
#' @rdname mp_plot_diff_manhattan-methods
#' @aliases mp_plot_diff_manhattan,MPSE
#' @export mp_plot_diff_manhattan
setMethod("mp_plot_diff_manhattan", signature(.data='MPSE'), .internal_mp_plot_diff_manhattan)

#' @rdname mp_plot_diff_manhattan-methods
#' @aliases mp_plot_diff_manhattan,tbl_mpse
#' @export mp_plot_diff_manhattan
setMethod("mp_plot_diff_manhattan", signature(.data="tbl_mpse"), .internal_mp_plot_diff_manhattan)

#' @rdname mp_plot_diff_manhattan-methods
#' @aliases mp_plot_diff_manhattan,grouped_df_mpse
#' @export mp_plot_diff_manhattan
setMethod("mp_plot_diff_manhattan", signature(.data="grouped_df_mpse"), .internal_mp_plot_diff_manhattan)

.get_annot_index <- function(x, taxa.class, tip.annot = TRUE){
    start.annot <- x %>%                 
       dplyr::filter(.data$nodeClass == taxa.class) %>%
       dplyr::pull(.data$nodeDepth) %>%      
       max(na.rm = TRUE)
    end.annot <- x %>% dplyr::pull(.data$nodeDepth) %>% max(na.rm = TRUE)
    if (tip.annot){                          
        annot.index <- seq(start.annot, end.annot)
    }else{                                   
        annot.index <- seq(start.annot, end.annot - 1)
    }
    return(annot.index)	
}

.generate_annot_df <- function(x, annot.index, .group, removeUnknown){
    abbre <- c(letters, toupper(letters))
    x1 <- x %>% dplyr::filter(.data$nodeDepth %in% annot.index & !is.na(!!.group))
    if (removeUnknown){
        x1 %<>% dplyr::filter(!grepl('__un_', .data$label))
    }
	if (nrow(x1) <= 52){
        x1[['.LABEL']] <- paste0(abbre[seq_len(nrow(x1))], ": ", x1$label)
        x1$label <- abbre[seq_len(nrow(x1))]	
    }else{
        x1[['.LABEL']] <- paste0(c(abbre, seq_len(nrow(x1) - 52))[seq_len(nrow(x1))], ": ", x1$label)
        x1$label <- c(abbre, seq_len(nrow(x1) - 52))[seq_len(nrow(x1))]
    }
    x2 <- x %>% dplyr::filter(!.data$node %in% x1$node)
	x <- dplyr::bind_rows(x1, x2) %>% dplyr::arrange(.data$node) %>% add_class("tbl_tree") 
	return (x)
}

#' Create the scale of mp_plot_diff_cladogram.
#' @param values a set of aesthetic values (different group (default)) to map data values to.
#' @param breaks One of 'NULL' for no breaks, ‘waiver()’ for the default breaks, A character vector of breaks.
#' @param na.value The aesthetic value to use for missing (‘NA’) values.
#' @param ... see also 'discrete_scale' of 'ggplot2'.
#' @importFrom ggplot2 waiver
#' @export
scale_fill_diff_cladogram  <- function(
    values, 
    breaks = waiver(), 
    na.value = 'grey50', 
    ...){
    
    if (missing(values)){
        values <- NULL
    }
    params <- list(...)
    scl <- list(
              values = values,
              breaks = breaks,
              na.value = na.value,
			  params = params
           )
    structure(scl, class = 'ScaleDiffClade')
    
}

#' set the color scale of plot generated by mp_plot_diff_boxplot
#' @param .data the aplot object generated by mp_plot_diff_boxplot.
#' @param values the color vector, required.
#' @param ... additional parameters, see also the 'scale_fill_manual' of 'ggplot2'
#' @export
set_diff_boxplot_color <- function(
    .data,
    values,
    ...){

    .data[[1]] <- .data[[1]] + scale_fill_manual(values = values, aesthetics='fill', ...)
    tmp.group <- .data[[1]]$data %>% dplyr::pull(!!.data[[1]]$layers[[1]]$mapping$fill) 
    fix.group <- unique(sort(tmp.group))
    if (!is.null(levels(tmp.group))){
        fix.group <- levels(tmp.group)
    }
    fix.group.color <- values
    if (is.null(names(fix.group.color))){
        names(fix.group.color) <- fix.group 
    }
    tmp.group2 <- .data[[2]]$data %>% dplyr::pull(!!.data[[2]]$layers[[2]]$mapping$colour)
    fix.group.color <- fix.group.color[match(unique(as.character(tmp.group2)), names(fix.group.color))]
    .data[[2]] <- .data[[2]] + scale_color_manual(values = fix.group.color)
    return(.data)
}


#' @importFrom ggplot2 ggplot_add
#' @method ggplot_add ScaleDiffClade
#' @export
ggplot_add.ScaleDiffClade <- function(object, plot, object_name){
    index <- which(plot$scales$find('fill'))
    ss <- plot$scales$scales[[index]]
    original.fill <- FALSE
    if (is.null(object$values)){
        object$values <- ss$palette.cache
        original.fill <- TRUE
    }
	params <- object$params
    object$params <- NULL
    object <- c(object, params)
    if (!'guide' %in% names(object)){
        object$guide <- ggplot2::guide_legend(order = 1, override.aes = list(alpha=1, size=3, shape=22))
    }else{
        if (object$guide$order %in% c(0, 2)){
            object$guide$order <- 1
        }
    }
    new.fill.scale <- do.call('scale_fill_manual', object)
    object$aesthetics <- 'colour_new'
    object$guide <- 'none'
    original.color <- do.call(scale_color_manual, object)
    color.label <- .build_color_values(plot, values=object$values)
    object$values <- color.label
    object$aesthetics <- 'colour'
    object$labels <- waiver()
    object$breaks <- names(color.label)
    object$guide <- ggplot2::guide_legend(
                       ncol = ifelse(length(color.label) > 30, 2, 1), 
                       override.aes = list(size = 3),
                       order = 2
                    )
    new.color.scale <- do.call('scale_color_manual', c(list(name = NULL), object))
    if (original.fill){
        object <- new.color.scale
    }else{
        object <- list(new.fill.scale, original.color, new.color.scale) 
    }
    ggplot_add(object, plot, object_name)
}

.build_color_values <- function(plot, values){
    idx <- which(unlist(lapply(plot$layers, function(x) any(grepl('GeomHilightRect', class(x[['geom']]))))))
    group <- plot$layers[[idx]]$mapping$fill    
    group.cols <- plot$data %>%
        dplyr::filter(!is.na(!!group)) %>%	
        dplyr::select(!!group) %>% 
        dplyr::distinct() %>%
        dplyr::arrange(!!group) %>%
        dplyr::mutate(colors=values[seq_len(nrow(.))])
    group.cols <- plot$data %>% 
        dplyr::filter(!is.na(.data$.LABEL)) %>% 
        dplyr::select(!!group, !!rlang::sym('.LABEL')) %>%
        dplyr::left_join(group.cols, by=rlang::as_name(group)) %>%
        dplyr::pull(.data$colors, name=.data$.LABEL)
    return(group.cols)
}


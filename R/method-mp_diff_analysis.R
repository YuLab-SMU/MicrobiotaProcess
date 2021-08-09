#' Differential expression analysis for MPSE or tbl_mpse object
#' @rdname mp_diff_analysis-methods
#' @param .data MPSE or tbl_mpse object
#' @param .abundance the name of abundance to be calculated
#' @param .group the group name of the samples to be calculated.
#' @param .sec.group the second group name of the samples to be calculated.
#' @param action character, "add" joins the new information to the taxatree and return,
#' "only" return a non-redundant tibble with the result of different analysis. "get" return
#' 'diffAnalysisClass' object.
#' @param tip.level character the taxa level to be as tip level
#' @param force logical whether to calculate the relative abundance forcibly when the abundance
#' is not be rarefied, default is FALSE.
#' @param relative logical whether calculate the relative abundance.
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
#' significantly different features, default is "generalizedFC".
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
#' @param seed a random seed to make the adonis analysis reproducible, default is 123.
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
#' mouse.time.mpse %>%
#'   mp_diff_analysis(.abundance=RareAbundance, 
#'                    .group=time, 
#'                    first.test.alpha=0.01,
#'                    action="get") %>%
#' ggdiffclade(linewd=0.1)
setGeneric("mp_diff_analysis", function(.data, 
                                        .abundance, 
                                        .group, 
                                        .sec.group=NULL, 
                                        action="add",
                                        tip.level="OTU",
                                        force=FALSE,
                                        relative=TRUE,
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

.internal_mp_diff_analysis <-  function(
              .data,
              .abundance,
              .group,
              .sec.group=NULL,
              action="add",
              tip.level="OTU",
              force=FALSE,
              relative=TRUE,
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
                 tibble::column_to_rownames(var="Sample")
                 #dplyr::mutate(across(!!.group, ~if(!is.factor(.x))as.factor(.x)))


     if (!rlang::quo_is_null(.sec.group)){
         sampleda %<>% 
             duplicatedtaxcheck() %>% 
             tibble::column_to_rownames(var="rowname")
         .sec.group <- rlang::as_name(.sec.group)
     }else{
         .sec.group <- NULL
     }
     
     if (relative){ 
         if(force){
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

     if (!any(grepl(abundance.nm, assaysvar))){
         .data %<>% mp_cal_abundance(.abundance=!!.abundance, force=force, relative=relative)
     }

     taxatree <- .data %>% mp_extract_tree(tip.level=tip.level)

     if (is.null(taxatree)){
         f_tb <- .data %>% mp_extract_assays(.abundance=!!abundance.nm, byRow=FALSE)
     }else{
         AbundBySample <- abundance.nm %>% 
                          gsub("^Rel", "", .) %>%
                          gsub("BySample", "", .) %>%
                          paste0("BySample")
         f_tb <- taxatree %>%
                 as_tibble() %>%
                 dplyr::filter(.data$nodeClass!="Root") %>%
                 dplyr::select(c("label", AbundBySample)) %>% 
                 tidyr::unnest(cols=AbundBySample) %>%
                 dplyr::select(c("label", "Sample", abundance.nm)) %>% 
                 distinct() %>%
                 tidyr::pivot_wider(id_cols="Sample", names_from="label", values_from=abundance.nm) %>%
                 tibble::column_to_rownames(var="Sample")
     }

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

     if (!length(first.test.sig.vars)>0){
          message("There are not significantly discriminative features after internal first test !")
          return(NULL)
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

     if (length(second.test.sig.vars)==0){
         message("There are not significantly discriminative features after internal second test!")
         return(NULL)
     }

     leaveclasslevels <- unlist(lapply(names(second.test.sig.vars), 
                                       function(x){unlist(strsplit(x,"-vs-"))}))

     second.test.sig.vars <- get_consistentfeatures(diffsubclassfeature=second.test.sig.vars, 
                                          classgroup=rlang::as_name(.group),
                                          classlevels=leaveclasslevels) 
     second.test.sig.vars.vectors <- second.test.sig.vars %>% get_secondvarlist()
     if (!is.null(normalization)){
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
         taxda$OTU <- rownames(taxda)
         indx <- match(tip.level, colnames(taxda))
         taxda %<>% select(seq_len(indx))
         res <- new("diffAnalysisClass", originalD=f_tb, sampleda=sampleda, taxda=taxda, result=result, kwres=first.res,
                   secondvars=second.test.sig.vars, mlres=ml.res, someparams=params)
         return(res)
     }else if (action=="only"){
         return(result)
     }else if (action=="add"){
         if (is.null(taxatree)){
             otu_tb <- .data %>% 
                       mp_extract_feature() 
             result %<>% dplyr::rename(label="f")
             #result %<>% dplyr::select(c("label",setdiff(colnames(result), colnames(otu_tb))))
             otu_tb %<>% dplyr::left_join(result, by="label")
             return(otu_tb)
         }else{
             newgroup <- paste0("Sign_", rlang::as_name(.group))
             result %<>% 
                 dplyr::rename(label="f", !!newgroup:=!!.group) 
             #result %<>% dplyr::select(setdiff(colnames(result), 
             #                          c(colnames(taxatree@data), 
             #                            colnames(taxatree@extraInfo)))
             #                         )
             taxatree %<>% treeio::full_join(result, by="label")
             return(taxatree)
         }
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

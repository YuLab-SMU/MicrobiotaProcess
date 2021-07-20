#setGeneric("mp_diff_analysis", function(.data, 
#                                        .abundance, 
#                                        .group, 
#                                        .sec.group=NULL, 
#                                        .tiplevel=OTU, 
#                                        .first.test.method=kruskal.test,
#                                        .first.test.alpha=0.05,
#                                        .p.adjust=fdr,
#                                        .strict = TRUE,
#                                        .second.test.method=wilcox.test,
#                                        .second.test.alpha=0.05,
#                                        .cl.min=5
#                                        .cl.test=TRUE,
#                                        .subcl.min = 3,
#                                        .subcl.test=TRUE,
#                                        .normalization=1000000,
#                                        .ldascore=2,
#                                        .bootnums=30, 
#                                        .ci=0.95, 
#                                        type="species",
#                                        ...
#                                        ){
#     standardGeneric("mp_diff_analysis")
#})
#
#
#setMethod("mp_diff_analysis", signature(.data="MPSE"),
#          function(
#              .data,
#              .abundance,
#              .group,
#              .sec.group=NULL,
#              .tiplevel=OTU,
#              .first.test.method=kruskal.test,
#              .first.test.alpha=0.05,
#              .p.adjust=fdr,
#              .strict = TRUE,
#              .second.test.method=wilcox.test,
#              .second.test.alpha=0.05,
#              .cl.min=5
#              .cl.test=TRUE,
#              .subcl.min = 3,
#              .subcl.test=TRUE,
#              .normalization=1000000,
#              .ldascore=2,
#              .bootnums=30,
#              .ci=0.95,
#              type="species",
#              ...
#              ){
#
#     .abundance <- rlang::enquo(.abundance)
#     .group <- rlang::enquo(.group)
#     .sec.group <- rlang::enquo(.sec.group)
#     .tiplevel <- rlang::enquo(.tiplevel)
#     .first.test.method <- rlang::enquo(.first.test.method)
#     .first.test.alpha <- rlang::enquo(.first.test.alpha)
#     .p.adjust <- rlang::enquo(.p.adjust)
#     .second.test.method <- rlang::enquo(.second.test.method)
#     .second.test.alpha <- rlang::enquo(.second.test.alpha)
#     
#     if (rlang::quo_is_missing(.abundance)){
#         
#     }
#
#
#
#})

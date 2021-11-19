#' @title Dropping Species with Few abundance and Few Occurrences
#' 
#' @description
#' Drop species or features from the feature data frame or phyloseq that occur
#' fewer than or equal to a threshold number of occurrences and fewer 
#' abundance than to a threshold abundance.
#'
#' @param obj object, phyloseq or a dataframe of species (n_sample, n_feature).
#' @param minocc numeric, the threshold number of occurrences to be 
#' dropped, if < 1.0,it will be the threshold ratios of occurrences, 
#' default is 0.
#' @param minabu numeric, the threshold abundance, if fewer than the 
#' threshold will be dropped, default is 0.
#' @param ..., additional parameters.
#' @return dataframe of new features.
#' @rdname drop_taxa
#' @export
#' @author Shuangbin Xu
#' @examples
#' \dontrun{
#' otudafile <- system.file("extdata", "otu_tax_table.txt",
#'                          package="MicrobiotaProcess")
#' otuda <- read.table(otudafile, sep="\t", 
#'                     header=TRUE, row.names=1, 
#'                     check.names=FALSE, skip=1, 
#'                     comment.char="")
#' otuda <- otuda[sapply(otuda, is.numeric)]
#' otuda <- data.frame(t(otuda), check.names=FALSE)
#' dim(otuda)
#' otudat <- drop_taxa(otuda, minocc=0.1, minabu=1)
#' dim(otudat)
#' data(test_otu_data)
#' test_otu_data %<>% as.phyloseq()
#' keepps <- drop_taxa(test_otu_data, minocc=0.1, minabu=0)
#' }
setGeneric("drop_taxa",function(obj, ...){standardGeneric("drop_taxa")})

#' @aliases drop_taxa,data.frame
#' @rdname drop_taxa
#' @export
setMethod("drop_taxa", "data.frame",
    function(obj, minocc=0, minabu=0,...){
    if (minocc < 1.0){
    	minocc <- round(dim(obj)[1]*minocc, 0)
    }
    obj <- obj[,apply(obj>minabu,2,sum)>=minocc,drop=FALSE]
    return (obj)
})

#' @aliases drop_taxa,phyloseq
#' @rdname drop_taxa
#' @export
setMethod("drop_taxa", "phyloseq", function(obj, ...){
    otuda <- checkotu(obj)
    otuda <- drop_taxa(otuda, ...)
    keepotu <- colnames(otuda)
    if (!is.null(obj@phy_tree)){
        tmptree <- obj@phy_tree
        removetip <- setdiff(tmptree$tip.label, keepotu)
        tmptree <- ape::drop.tip(tmptree, removetip)
        obj@phy_tree <- tmptree
    }
    if (!is.null(obj@refseq)){
        tmpseq <- obj@refseq
        keepseq <- intersect(keepotu, names(tmpseq))
        tmpseq <- tmpseq[keepseq]
        obj@refseq <- tmpseq
    }
    if (!is.null(obj@tax_table)){
        tmptax <- obj@tax_table
        keeptax <- intersect(keepotu, rownames(tmptax))
        tmptax <- tmptax[match(keeptax,rownames(tmptax)),,drop=FALSE]
        obj@tax_table <- tmptax
    }
    obj@otu_table <- phyloseq::otu_table(otuda,taxa_are_rows=FALSE)
    return(obj)
})

#' Filter OTU (Features) By Abundance Level
#' @rdname mp_filter_taxa-methods
#' @param .data MPSE or tbl_mpse or grouped_df_mpse object.
#' @param .abundance the column names of abundance, default is NULL, 
#' meaning the 'Abundance' column.
#' @param min.abun numeric minimum abundance required for each one sample
#' default is 0 (.abundance=Abundance or NULL), meaning the abundance of
#' OTU (Features) for each one sample should be >= 0.
#' @param min.prop numeric minimum proportion of samples that contains the OTU (Features)
#' when min.prop larger than 1, meaning the minimum number of samples that contains 
#' the OTU (Features).
#' @param include.lowest logical whether include the lower boundary of \code{min.abun}
#' default is FALSE ( > \code{min.abun}), if it is TRUE, meaning (>= \code{min.abun}).
#' @param ... additional parameters, meaningless now.
#' @author Shuangbin Xu
#' @export
#' @examples
#' data(mouse.time.mpse)
#' mouse.time.mpse %>% mp_filter_taxa(.abundance=Abundance, min.abun=1, min.prop=1)
#' # For tbl_mpse object.
#' mouse.time.mpse %>% as_tibble %>% mp_filter_taxa(.abundance=Abundance, min.abun=1, min.prop=1)
#' # This also can be done using group_by, filter of dplyr.
#' mouse.time.mpse %>% 
#'  dplyr::group_by(OTU) %>% 
#'  dplyr::filter(sum(Abundance>=1)>=1)
setGeneric("mp_filter_taxa", 
  function(
    .data, 
    .abundance = NULL, 
    min.abun = 0,
    min.prop = 0.05, 
    include.lowest = FALSE,
    ...)
  standardGeneric("mp_filter_taxa")
)

#' @rdname mp_filter_taxa-methods
#' @aliases mp_filter_taxa,MPSE
#' @exportMethod mp_filter_taxa
setMethod("mp_filter_taxa", 
  signature(.data="MPSE"), 
  function(
    .data, 
    .abundance = NULL, 
    min.abun = 0, 
    min.prop = 0.05, 
    include.lowest = FALSE,
    ...){

    .abundance <- rlang::enquo(.abundance)

    if (rlang::quo_is_null(.abundance)){
        .abundance <- as.symbol("Abundance")
    }

    assayda <- .data %>% 
               mp_extract_assays(.abundance=!!.abundance)

    if (min.prop < 1){
        min.prop <- ncol(assayda) * min.prop
    }
    
    if (include.lowest){
        assayda <- assayda[rowSums(assayda >= min.abun) >= min.prop, ,drop=FALSE]
    }else{
        assayda <- assayda[rowSums(assayda > min.abun) >= min.prop, ,drop=FALSE]
    }
    .data <- .data[rownames(.data) %in% rownames(assayda),,drop=FALSE]
    
    return(.data)
})

.internal_filter_taxa <- function(.data, .abundance = NULL, min.abun = 0, min.prop = 0.05, include.lowest = FALSE, ...){

    .abundance <- rlang::enquo(.abundance)
    flag.group <- inherits(.data, "grouped_df_mpse")

    if (rlang::quo_is_null(.abundance)){
        .abundance <- as.symbol("Abundance")
    }

    if (min.prop < 1){
        min.prop <- .data %>% 
                    pull("Sample") %>% 
                    unique() %>% 
                    length() %>%
                    magrittr::multiply_by(min.prop)
    }
    
    if (flag.group){
        tmpgroups <- attr(.data, "groups")
        groupvars <- names(tmpgroups)[names(tmpgroups) != ".rows"]
        groupvars <- lapply(groupvars, function(i) as.symbol(i))
        .data %<>% ungroup
    }

    if (include.lowest){
        .data %<>% 
              dplyr::group_by(!!rlang::sym("OTU")) %>%
              dplyr::filter(sum(!!.abundance >= min.abun)>= min.prop) %>% 
              ungroup()
    }else{
        .data %<>% 
              dplyr::group_by(!!rlang::sym("OTU")) %>%
              dplyr::filter(sum(!!.abundance > min.abun) >= min.prop) %>%
              ungroup()
    }

    if (flag.group){
        .data <- do.call(group_by, c(list(.data), groupvars)) 
    }
    return(.data)
}

#' @rdname mp_filter_taxa-methods
#' @aliases mp_filter_taxa,tbl_mpse
#' @exportMethod mp_filter_taxa
setMethod("mp_filter_taxa", signature(.data="tbl_mpse"), .internal_filter_taxa)

#' @rdname mp_filter_taxa-methods
#' @aliases mp_filter_taxa,grouped_df_mpse
#' @exportMethod mp_filter_taxa
setMethod("mp_filter_taxa", signature(.data="grouped_df_mpse"), .internal_filter_taxa)

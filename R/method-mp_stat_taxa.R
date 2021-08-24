#' Count the number and total number taxa for each sample at different taxonomy levels
#' @rdname mp_stat_taxa-methods
#' @param .data MPSE or tbl_mpse object
#' @param .abundance the column name of abundance to be calculated
#' @param action a character "get" returns a table only contained the number and 
#' total number for each sample at different taxonomy levels, "only" returns a 
#' non-redundant tibble contained a nest column (StatTaxaInfo) and other sample information, 
#' "add" returns a update object (.data) contained a nest column (StatTaxaInfo). 
#' @param ... additional parameter
#' @return update object or tbl_df according action argument
#' @author Shuangbin Xu
#' @export
#' @examples
#' data(mouse.time.mpse)
#' mouse.time.mpse %>%
#'  mp_stat_taxa(.abundance=Abundance, action="only")
setGeneric("mp_stat_taxa", function(.data, .abundance, action="add", ...) standardGeneric("mp_stat_taxa"))

.internal_stat_taxa <- function(.data, .abundance, action="add", ...){
    action %<>% match.arg(c("add", "only", "get"))
    .abundance <- rlang::enquo(.abundance)
    if (rlang::quo_is_missing(.abundance)){
        .abundance <- as.symbol("Abundance")
    }

    clnm1 <- paste0("NumBy", rlang::as_name(.abundance), "EachTaxonomy")
    clnm2 <- paste0("TotalNumBy", rlang::as_name(.abundance), "EachTaxonomy")

    dat <- .data %>% 
           mp_cal_abundance(.abundance=!!.abundance, force=TRUE, action="only") %>%
           tidyr::unnest(cols=paste0(rlang::as_name(.abundance) %>% 
                                     gsub("^Rel", "", .) %>%
                                     gsub("BySample$", "", .), 
                                 "BySample")
             ) %>%
           dplyr::group_by(.data$Sample, .data$nodeClass) %>%
           dplyr::summarize(!!clnm1:=sum(!!.abundance>0), !!clnm2:=sum(!!.abundance))

    if (action=="get"){
        return(dat)
    }

    dat %<>% 
        tidyr::nest(!.data$Sample) %>% 
        suppressWarnings() %>%
        rename(StatTaxaInfo="data")

    if (action == "only"){
        return(dat)
    }else if (action=="add"){
        if (inherits(.data, "MPSE")){
            .data@colData <- .data %>% 
                             mp_extract_sample() %>%
                             dplyr::left_join(dat, by="Sample", suffix=c("", ".y")) %>%
                             column_to_rownames(var="Sample") %>% 
                             S4Vectors::DataFrame(check.names=FALSE)
        }else{
            samplevar <- .data %>% attr("samplevar")
            assaysvar <- .data %>% attr("assaysvar")
            othernm <- colnames(.data)[!colnames(.data) %in% c("OTU", "Sample", assaysvar, samplevar)]
            .data %<>% left_join(dat, by="Sample", suffix=c("", ".y")) %>%
                       select(c("OTU", "Sample", assaysvar, samplevar, colnames(dat), othernm))
        }
        return(.data)
    }
}

#' @rdname mp_stat_taxa-methods
#' @aliases mp_stat_taxa,MPSE
#' @exportMethod mp_stat_taxa
setMethod("mp_stat_taxa", signature(.data="MPSE"), .internal_stat_taxa)

#' @rdname mp_stat_taxa-methods
#' @aliases mp_stat_taxa,tbl_mpse
#' @exportMethod mp_stat_taxa
setMethod("mp_stat_taxa", signature(.data="tbl_mpse"), .internal_stat_taxa)

#' @rdname mp_stat_taxa-methods
#' @aliases mp_stat_taxa,grouped_df_mpse
#' @exportMethod mp_stat_taxa
setMethod("mp_stat_taxa", signature(.data="grouped_df_mpse"), .internal_stat_taxa)

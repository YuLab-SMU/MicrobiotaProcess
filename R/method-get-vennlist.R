#' @title generate a vennlist for VennDiagram 
#' @param obj phyloseq, phyloseq class or data.frame
#' a dataframe contained one character column and the others are numeric.
#' or all columns should be numeric if sampleinfo isn't NULL.
#' @param sampleinfo dataframe; a sample information, default is NULL.
#' @param  factorNames character, a column name of sampleinfo, 
#' when sampleinfo isn't NULL, factorNames shouldn't be NULL, default is NULL,
#' when the input is phyloseq, the factorNames should be provided. 
#' @param ..., additional parameters
#' @return return a list for VennDiagram.
#' @author Shuangbin Xu
#' @export 
#' @examples
#' \dontrun{
#' data(test_otu_data)
#' vennlist <- get_vennlist(test_otu_data, 
#'                  factorNames="group")
#' vennlist
#' library(VennDiagram)
#' venn.diagram(vennlist, height=5, 
#'              width=5, filename = "./test_venn.pdf", 
#'              alpha = 0.85, fontfamily = "serif", 
#'              fontface = "bold",cex = 1.2, 
#'              cat.cex = 1.2, cat.default.pos = "outer",
#'              cat.dist = c(0.22,0.22,0.12,0.12), 
#'              margin = 0.1, lwd = 3, 
#'              lty ='dotted', 
#'              imagetype = "pdf")
#' }
setGeneric("get_vennlist", function(obj, ...)standardGeneric("get_vennlist"))

#' @aliases get_vennlist,phyloseq
#' @rdname get_vennlist
#' @export 
setMethod("get_vennlist", "phyloseq", function(obj, factorNames, ...){
    otuda <- checkotu(obj)
    sampleda <- checksample(obj)
    #tmpfactors <- colnames(sampleda)[factorNamesIndex]
    if(is.null(factorNames)){stop("The object is phyloseq, factorNames should not be NULL.")}
    vennlist <- get_vennlist(obj=otuda,
                             sampleinfo=sampleda,
                             factorNames=factorNames, ...)
    return(vennlist)
})

#' @aliases get_vennlist,data.framet
#' @rdname get_vennlist
#' @export
setMethod("get_vennlist", "data.frame", function(obj,
    sampleinfo=NULL,
    factorNames=NULL,...){
    if (!is.null(sampleinfo) && !is.null(factorNames)){
    	sampleinfo <- sampleinfo[,match(factorNames, colnames(sampleinfo)), 
    							 drop=FALSE]
    }
    if (!is.null(sampleinfo) && is.null(factorNames)){
    	stop("when sampleinfo isn't NULL, factorNames shouldn't be NULL")
    }
    obj <- get_count(data=obj, 
                     featurelist=sampleinfo)
    vennlist <- apply(obj, 1, function(x){names(x[x>0])})
    return(vennlist)
})


#' Calculating the OTU for each sample or group, the result can be visualized by 'ggVennDiagram'
#' @rdname mp_cal_venn-methods
#' @param .data MPSE or tbl_mpse object
#' @param .group the name of group to be calculated.
#' if it is no provided, the sample will be used.
#' @param .abundance the name of otu abundance to be calculated.
#' if it is null, the rarefied abundance will be used.
#' @param action character, "add" joins the new information to the tibble of tbl_mpse or 
#' rowData of MPSE. "only" and "get" return a non-redundant tibble with the just new information. 
#' @param force logical whether calculate the relative abundance forcibly when the abundance
#' is not be rarefied, default is FALSE.
#' @param ... additional parameters.
#' @return update object or tibble according the 'action'
#' @export
#' @author Shuangbin Xu
#' @examples
#' data(mouse.time.mpse)
#' mouse.time.mpse %>%
#' mp_rrarefy() %>%
#' mp_cal_venn(.abundance=RareAbundance, .group=time, action="add") -> mpse
#' mpse
#' library(ggplot2)
#' mpse %>% 
#'   mp_extract_sample() %>% 
#'   select(time, vennOftime) %>%
#'   distinct() %>%
#'   pull(var=vennOftime, name=time) %>%
#'   ggVennDiagram::ggVennDiagram()
setGeneric("mp_cal_venn", function(.data, .group, .abundance=NULL, action="add", force=FALSE, ...)standardGeneric("mp_cal_venn"))

#' @rdname mp_cal_venn-methods
#' @aliases mp_cal_venn,MPSE
#' @exportMethod mp_cal_venn
setMethod("mp_cal_venn", signature(.data="MPSE"), function(.data, .group, .abundance=NULL, action="add", force=FALSE, ...){
    
    .abundance <- rlang::enquo(.abundance)
    .group <- rlang::enquo(.group)

    action %<>% match.arg(c("add", "get", "only"))

    if (rlang::quo_is_null(.abundance)){
        .abundance <- as.symbol("RareAbundance")
    }

    if (rlang::quo_is_missing(.group)){
        .group <- as.symbol("Sample")
    }

    if (!valid_rare(.data, .abundance=.abundance) && !force){
        glue::glue("The rarefied abundance of species might not be provided. Rarefaction of all
                    observations is performed automatically. If you still want to calculate the
                    alpha index with the '.abundance', you can set 'force=TRUE'. ")
        .data <- mp_rrarefy(.data=.data, ...)
        .abundance <- as.symbol("RareAbundance")
    }

    xx <- SummarizedExperiment::assays(.data)@listData

    da <- xx[[rlang::as_name(.abundance)]] %>% 
          tibble::as_tibble(rownames="OTU") %>%
          tidyr::pivot_longer(!as.symbol("OTU"), 
                              names_to="Sample", 
                              values_to=rlang::as_name(.abundance))

    sampleda <- .data %>%
                mp_extract_sample()

    if (ncol(sampleda)>1){
        da %<>% left_join(sampleda, by="Sample")
    }

    dat <- da %>% 
           .internal_cal_venn(.abundance=.abundance, .group=.group)

    if (action == "add"){
        .data@colData <- 
             sampleda %>%
             dplyr::left_join(dat, by=rlang::as_name(.group)) %>%
             tibble::column_to_rownames(var="Sample") %>%
             S4Vectors::DataFrame(check.names=FALSE)

        return(.data)
    
    }else if (action == "only"){
        return (dat)
    }else if (action == "get"){
        return (dat)
    }
})


.internal_mp_cal_venn <- function(.data, .group, .abundance=NULL, action="add", force=FALSE, ...){
    .abundance <- rlang::enquo(.abundance)
    .group <- rlang::enquo(.group)

    action %<>% match.arg(c("add", "get", "only"))

    if (rlang::quo_is_null(.abundance)){
        .abundance <- as.symbol("RareAbundance")
    }

    if (rlang::quo_is_missing(.group)){
        .group <- as.symbol("Sample")
    }

    if (!valid_rare(.data, .abundance=.abundance) && !force){
        glue::glue("The rarefied abundance of species might not be provided. Rarefaction of all
                    observations is performed automatically. If you still want to calculate the
                    alpha index with the '.abundance', you can set 'force=TRUE'. ")
        .data <- mp_rrarefy(.data=.data, ...)
        .abundance <- as.symbol("RareAbundance")
    }

    dat <- .data %>% 
          dplyr::ungroup() %>%
          dplyr::select(!!as.symbol("OTU"), !!.group, !!.abundance) %>%
          .internal_cal_venn(.abundance=.abundance, .group=.group) 
    
    if (action=="add"){
        .data %<>% 
            dplyr::left_join(dat, by=rlang::as_name(.group))
        return(.data)
    }else if(action == "get"){
        return (dat)
    }else if(action == "get"){
        return (dat)
    }
}


.internal_cal_venn <- function(.data, .abundance, .group){
    vennnm <- paste0("vennOf", rlang::as_name(.group))

    dat <- .data %>% 
        dplyr::group_by(!!as.symbol("OTU"), !!.group) %>%
        dplyr::summarize(AbundanceBy=sum(!!.abundance)) %>%
        suppressMessages() %>% 
        dplyr::filter(!!as.symbol("AbundanceBy")>0) %>%
        dplyr::select(- !!as.symbol("AbundanceBy")) %>%
        group_by(!!.group) %>% 
        dplyr::summarize(across(!!as.symbol("OTU"), list, .names=vennnm))
    return(dat)
}

#' @rdname mp_cal_venn-methods
#' @aliases mp_cal_venn,tbl_mpse
#' @exportMethod mp_cal_venn
setMethod("mp_cal_venn", signature(.data = "tbl_mpse"), .internal_mp_cal_venn)

#' @rdname mp_cal_venn-methods
#' @aliases mp_cal_venn,grouped_df_mpse
#' @exportMethod mp_cal_venn
setMethod("mp_cal_venn", signature(.data = "grouped_df_mpse"), .internal_mp_cal_venn)

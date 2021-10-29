#' @title generate the dataset for upset of UpSetR
#' @param obj object, phyloseq or data.frame, if it is data.frame, 
#' the shape of it should be row sample * columns features.
#' @param sampleda data.frame, if the obj is data.frame, the sampleda
#' should be provided.
#' @param factorNames character, the column names of factor in sampleda
#' @param threshold integer, default is 0.
#' @param ..., additional parameters.
#' @return a data.frame for the input of `upset` of `UpSetR`.
#' @author Shuangbin Xu
#' @rdname get_upset
#' @export
#' @examples
#' \dontrun{
#' data(test_otu_data)
#' upsetda <- get_upset(test_otu_data, factorNames="group")
#' otudafile <- system.file("extdata", "otu_tax_table.txt",
#'                          package="MicrobiotaProcess")
#' samplefile <- system.file("extdata","sample_info.txt", 
#'                          package="MicrobiotaProcess")
#' otuda <- read.table(otudafile, sep="\t", header=TRUE, 
#'                     row.names=1, check.names=FALSE,
#'                     skip=1, comment.char="")
#' sampleda <- read.table(samplefile,sep="\t", 
#'                        header=TRUE, row.names=1)
#' head(sampleda)
#' otuda <- otuda[sapply(otuda, is.numeric)]
#' otuda <- data.frame(t(otuda), check.names=FALSE)
#' head(otuda[1:5, 1:5])
#' upsetda2 <- get_upset(obj=otuda, sampleda=sampleda, 
#'                      factorNames="group")
#' #Then you can use `upset` of `UpSetR` to visualize the results.
#' library(UpSetR)
#' upset(upsetda, sets=c("B","D","M","N"), sets.bar.color = "#56B4E9",
#'       order.by = "freq", empty.intersections = "on")
#' }
setGeneric("get_upset", function(obj, ...)standardGeneric("get_upset"))

#' @aliases get_upset,data.frame
#' @rdname get_upset
#' @importFrom stats na.omit
#' @export
setMethod("get_upset", "data.frame", function(obj, sampleda, factorNames, threshold=0){
    flaglen <- length(na.omit(match(rownames(obj),rownames(sampleda))))
    sampleda <- sampleda[,match(factorNames, colnames(sampleda)),drop=FALSE]
    if (flaglen==0){
        stop("The sample names of obj and sampleda should be consistent!
              Please check the rownames of obj and rownames of sampleda!") 
    }
    if (flaglen!=0 & flaglen < nrow(obj)){
        message("There are some sample names are not consistent!")
    }
    dameta <- merge(obj, sampleda, by=0)
    rownames(dameta) <- as.vector(dameta$Row.names)
    dameta$Row.names <- NULL
    dameta <- get_count(dameta)
    daupset <- apply(dameta, 1, 
                     function(x){unlist(lapply(x, function(x){if(x>threshold){1}else{0}}))})
    daupset <- data.frame(daupset, check.names=FALSE)
    return(daupset)
})

#' @aliases get_upset,phyloseq
#' @rdname get_upset
#' @export
setMethod("get_upset", "phyloseq", function(obj,...){
    otuda <- checkotu(obj)
    sampledata <- checksample(obj)
    daupset <- get_upset(obj=otuda, sampleda=sampledata,...)
    return(daupset)
})


#' Calculating the samples or groups for each OTU, the result can be visualized by 'ggupset'
#' @rdname mp_cal_upset-methods
#' @param .data MPSE or tbl_mpse object
#' @param .group the name of group to be calculated.
#' if it is no provided, the sample will be used.
#' @param .abundance the name of otu abundance to be calculated.
#' if it is null, the rarefied abundance will be used.
#' @param action character, "add" joins the new information to the tibble of tbl_mpse or 
#' rowData of MPSE. "only" and "get" return a non-redundant tibble with the just new information. 
#' which is a treedata object.
#' @param force logical whether calculate the relative abundance forcibly when the abundance
#' is not be rarefied, default is FALSE.
#' @param ... additional parameters.
#' @return update object or tibble according the 'action'
#' @seealso [mp_plot_upset()] 
#' @export
#' @author Shuangbin Xu
#' @examples
#' data(mouse.time.mpse)
#' mpse <- mouse.time.mpse %>%
#'         mp_rrarefy() %>%
#'         mp_cal_upset(.abundance=RareAbundance, .group=time, action="add")
#' mpse
#' library(ggplot2)
#' library(ggupset)
#' p <- mpse %>% mp_plot_upset(.group=time, .upset=ggupsetOftime)
#' p
#' # or set action="only"
#' \dontrun{
#' tbl <- mouse.time.mpse %>% 
#'        mp_rrarefy() %>% 
#'        mp_cal_upset(.abundance=RareAbundance, .group=time, action="only") 
#' tbl
#' p2 <- tbl %>%
#'       ggplot(aes(x=ggupsetOftime)) +
#'       geom_bar() +
#'       ggupset::scale_x_upset() +
#'       ggupset::theme_combmatrix(combmatrix.label.extra_spacing=30)
#'}
setGeneric("mp_cal_upset", function(.data, .group, .abundance=NULL, action="add", force=FALSE, ...)standardGeneric("mp_cal_upset"))

#' @rdname mp_cal_upset-methods
#' @aliases mp_cal_upset,MPSE
#' @exportMethod mp_cal_upset
setMethod("mp_cal_upset", signature(.data="MPSE"), function(.data, .group, .abundance=NULL, action="add", force=FALSE, ...){
    
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
        trash <- try(silent = TRUE,
                     expr = {
                         .data <- mp_rrarefy(.data = .data, ...)
                     }
                 )
        if (inherits(trash, "try-error")){
            stop_wrap("The 'Abundance' column cannot be rarefied, please check whether it is integer (count).
                       Or you can set 'force=TRUE' to calculate the result of 'upset' without rarefaction.
                      ")
        }

        message_wrap("The rarefied abundance of species might not be provided. Rarefaction of all
                      observations is performed automatically using 'Abundance' column.
                      If you still want to calculate the result of 'upset' with the specified '.abundance',
                      you can set 'force=TRUE'. ")
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
        da %<>% left_join(sampleda, by="Sample", suffix=c("", ".y"))
    }

    dat <- da %>% 
           .internal_cal_upset(.abundance=.abundance, .group=.group)

    if (action == "add"){
        SummarizedExperiment::rowData(.data) <- 
             SummarizedExperiment::rowData(.data) %>%
             avoid_conflict_names() %>%
             tibble::as_tibble(rownames="OTU") %>%
             dplyr::left_join(dat, by="OTU", suffix=c("", ".y")) %>%
             tibble::column_to_rownames(var="OTU")

        return(.data)
    
    }else if (action == "only"){
        return (dat)
    }else if (action == "get"){
        return (dat)
    }
})


.internal_mp_cal_upset <- function(.data, .group, .abundance=NULL, action="add", force=FALSE, ...){
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
          .internal_cal_upset(.abundance=.abundance, .group=.group) 
    
    if (action=="add"){
        .data %<>% 
            dplyr::left_join(dat, by="OTU", suffix=c("", ".y"))
        return(.data)
    }else if(action == "only"){
        return (dat)
    }else if(action == "get"){
        return (dat)
    }
}


.internal_cal_upset <- function(.data, .abundance, .group){
    upsetnm <- paste0("ggupsetOf", rlang::as_name(.group))

    dat <- .data %>% 
        dplyr::mutate(across(!!.group, as.vector)) %>%
        dplyr::group_by(!!as.symbol("OTU"), !!.group) %>%
        dplyr::summarize(AbundanceBy=sum(!!.abundance)) %>%
        suppressMessages() %>% 
        dplyr::filter(!!as.symbol("AbundanceBy")>0) %>%
        dplyr::select(- !!as.symbol("AbundanceBy")) %>%
        group_by(!!as.symbol("OTU")) %>% 
        dplyr::summarize(across(!!.group, list, .names=upsetnm))
    return(dat)
}

#' @rdname mp_cal_upset-methods
#' @aliases mp_cal_upset,tbl_mpse
#' @exportMethod mp_cal_upset
setMethod("mp_cal_upset", signature(.data = "tbl_mpse"), .internal_mp_cal_upset)

#' @rdname mp_cal_upset-methods
#' @aliases mp_cal_upset,grouped_df_mpse
#' @exportMethod mp_cal_upset
setMethod("mp_cal_upset", signature(.data = "grouped_df_mpse"), .internal_mp_cal_upset)

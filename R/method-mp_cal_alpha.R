#' @title alpha index
#' 
#' @description
#' calculate the alpha index (Obseve,Chao1,Shannon,Simpson) of sample
#' with \code{\link[vegan]{diversity}}
#' @param obj object, data.frame of (nrow sample * ncol taxonomy(feature)) 
#' or phyloseq.
#' @param mindepth numeric, Subsample size for rarefying community.
#' @param sampleda data.frame,sample information, row sample * column factors.
#' @param force logical whether calculate the alpha index even the count of otu is
#' not rarefied, default is FALSE. If it is TRUE, meaning the rarefaction is not be
#' performed automatically.
#' @param ... additional arguments.
#' @return data.frame contained alpha Index.
#' @author Shuangbin Xu
#' @rdname get_alphaindex
#' @export
#' @examples
#' \dontrun{
#' otudafile <- system.file("extdata", "otu_tax_table.txt", 
#'                         package="MicrobiotaProcess")
#' otuda <- read.table(otudafile, sep="\t", 
#'              header=TRUE, row.names=1, 
#'              check.names=FALSE, skip=1, comment.char="")
#' otuda <- otuda[sapply(otuda, is.numeric)] %>% t() %>% 
#'           data.frame(check.names=FALSE)
#' set.seed(1024)
#' alphatab <- get_alphaindex(otuda)
#' head(as.data.frame(alphatab))
#' data(test_otu_data)
#' class(test_otu_data)
#' set.seed(1024)
#' alphatab2 <- get_alphaindex(test_otu_data)
#' head(as.data.frame(alphatab2))
#' }
setGeneric("get_alphaindex", function(obj, ...){standardGeneric("get_alphaindex")})

#' @aliases get_alphaindex,matrix
#' @rdname get_alphaindex
#' @importFrom vegan estimateR diversity specnumber
#' @export
setMethod("get_alphaindex", "matrix", function(obj, mindepth, sampleda, force=FALSE, ...){
    #if (!identical(all.equal(obj, round(obj)),TRUE)){
    #    stop("the data should be integer (counts)!")
    #}
    if (missing(mindepth) || is.null(mindepth)){
        mindepth <- min(rowSums(obj))
    }
    if (obj %>% rowSums %>% var != 0 && !force){
        obj <- vegan::rrarefy(obj, mindepth)
    }
    Shannon <- diversity(obj)
    Simpson <- diversity(obj, index="simpson")
    J <- Shannon/log(specnumber(obj))
    if (identical(all.equal(obj, round(obj)),TRUE)){
        spn <- estimateR(obj)
        Observe <- spn[1,]
        Chao1 <- spn[2, ] 
        ACE <- spn[4, ]
        alpha <- data.frame(Observe=Observe, Chao1=Chao1, ACE=ACE, Shannon, Simpson, J)
    }else{
        Observe <- apply(obj, 1, function(x)sum(x>0))
        Chao1 <- NULL
        ACE <- NULL
        alpha <- data.frame(Observe=Observe, Shannon, Simpson, J)
    }
    if (missing(sampleda)){
        sampleda <- NULL
    }
    alpha <- new("alphasample",
                 alpha=alpha,
		 sampleda=sampleda)
    return(alpha)
})     

#' @aliases get_alphaindex,data.frame
#' @rdname get_alphaindex
#' @export
setMethod("get_alphaindex", "data.frame", function(obj, ...){
    obj <- obj[,colSums(obj)>0,drop=FALSE]
    obj <- as.matrix(obj)
    alpha <- get_alphaindex(obj, ...)
    return(alpha)
})

#' @aliases get_alphaindex,integer
#' @rdname get_alphaindex
#' @export
setMethod("get_alphaindex", "integer", function(obj, ...){
    obj <- obj[obj>0]
    obj <- as.matrix(obj)
    alpha <- get_alphaindex(obj, ...)
    alpha <- alpha@alpha
    return(alpha)
})

#' @aliases get_alphaindex,numeric
#' @rdname get_alphaindex
#' @export
setMethod("get_alphaindex", "numeric",function(obj, ...){
    obj <- obj[obj>0]
    obj <- as.matrix(obj)
    alpha <- get_alphaindex(obj, ...)
    alpha <- alpha@alpha
    return(alpha)
})

#' @aliases get_alphaindex,phyloseq
#' @rdname get_alphaindex
#' @export
setMethod("get_alphaindex", "phyloseq", function(obj, ...){
    otuda <- checkotu(obj)
    sampleda <- checksample(obj)
    alpha <- get_alphaindex(obj=otuda, sampleda=sampleda,...)
    return(alpha)
})

#' @title calculate the alpha index with MPSE or tbl_mpse
#' @docType methods
#' @rdname mp_cal_alpha-methods
#' @param .data MPSE or tbl_mpse object
#' @param .abundance The column name of OTU abundance column to be calculate
#' @param action character it has three options, "add" joins the new information 
#' to the input tbl (default), "only" return a non-redundant tibble with the just 
#' new information, ang 'get' return a 'alphasample' object.
#' @param force logical whether calculate the alpha index even the '.abundance' is 
#' not rarefied, default is FALSE.
#' @param ... additional arguments
#' @return update object or other (refer to action)
#' @export
#' @author Shuangbin Xu
#' @examples
#' data(mouse.time.mpse)
#' mpse <- mouse.time.mpse %>% 
#'         mp_rrarefy() %>%
#'         mp_cal_alpha(.abundance=RareAbundance)
#' tbl <- mpse %>% 
#'        mp_extract_sample
#' tbl
#' tbl %<>% 
#'   tidyr::pivot_longer(cols=!c("Sample", "time"), names_to="measure", values_to="alpha")
#' tbl
#' library(ggplot2)
#' library(ggsignif)
#' library(gghalves)
#' p <- ggplot(data=tbl, aes(x=time, y=alpha, fill=time)) + 
#'      geom_half_violin(color=NA, side="l", trim=FALSE) + 
#'      geom_boxplot(aes(color=time), fill=NA, position=position_nudge(x=.22), width=0.2) + 
#'      geom_half_point(side="r", shape=21) + 
#'      geom_signif(comparisons=list(c("Early", "Late")), test="wilcox.test", textsize=2) + 
#'      facet_wrap(facet=vars(measure), scales="free_y", nrow=1) +
#'      scale_fill_manual(values=c("#00A087FF", "#3C5488FF")) + 
#'      scale_color_manual(values=c("#00A087FF", "#3C5488FF"))
#' p
setGeneric("mp_cal_alpha", function(.data, .abundance=NULL, 
                                    action=c("add", "only", "get"), 
                                    force=FALSE, ...){
    standardGeneric("mp_cal_alpha")
})


#' @rdname mp_cal_alpha-methods
#' @aliases mp_cal_alpha,MPSE
#' @exportMethod mp_cal_alpha 
setMethod("mp_cal_alpha", signature(.data="MPSE"),
          function(.data, .abundance=NULL, action="add", force=FALSE, ...){
    action %<>% match.arg(c("add", "only", "get"))
    
    .abundance <- rlang::enquo(.abundance)
    if (rlang::quo_is_null(.abundance)){
        .abundance <- as.symbol("RareAbundance")
    }

    if (!valid_rare(.data, .abundance=.abundance) && !force){
        glue::glue("The rarefied abundance of species might not be provided. Rarefaction of all
                    observations is performed automatically. If you still want to calculate the
                    alpha index with the '.abundance', you can set 'force=TRUE'. ")
        .data <- mp_rrarefy(.data=.data, ...)
        .abundance <- as.symbol("RareAbundance")
    }
    
    alphada <- SummarizedExperiment::assays(.data)@listData[[rlang::as_name(.abundance)]] %>% 
               t() %>% 
               get_alphaindex(force=TRUE)
    if (action=="get"){
        alphada@sampleda <- .data@colData %>% data.frame(check.names=FALSE)
        return(alphada)
    }
    da <- alphada@alpha %>% as_tibble(rownames="Sample")
    da <- .data %>%
          mp_extract_sample() %>%
          left_join(da, by="Sample", suffix=c("", ".y")) 
    if (action=="add"){
        .data@colData <- da %>% 
                         column_to_rownames(var="Sample") %>% 
                         S4Vectors::DataFrame(check.names=FALSE)
        return(.data)
    }else if (action=="only"){
        return(da)
    }else{
        stop("Note: The action must be either \"add\" for adding this information to your data frame or \"get\" to just get the information")
    }
})

.internal_mp_cal_alpha <- function(.data, .abundance=NULL, action="add", force=FALSE, ...){
    action %<>% match.arg(c("add", "only", "get"))
    
    .abundance <- rlang::enquo(.abundance)

    if (rlang::quo_is_null(.abundance)){
        .abundance <- as.symbol("RareAbundance")
    }

    if (!valid_rare(.data, .abundance=.abundance) && !force){
        glue::glue("The rarefied abundance of species might not be provided. Rarefaction of all 
                    observations is performed automatically. If you still want to calculate the 
                    alpha index with the '.abundance', you can set 'force=TRUE'. ")
        .data <- mp_rrarefy(.data=.data, ...)
        .abundance <- as.symbol("RareAbundance")
    }

    alphada <- .data %>% 
               select(c("Sample", "OTU", rlang::as_name(.abundance))) %>%
               tidyr::pivot_wider(id_cols="Sample", 
                                  names_from="OTU", 
                                  values_from=rlang::as_name(.abundance)) %>%
               column_to_rownames(var="Sample") %>%
               get_alphaindex(force=TRUE)

    if (action=="get"){
        samplevar <- attr(.data, "samplevar")
        alphada@sampleda <- .data %>% 
                            select(samplevar) %>% 
                            column_to_rownames(var="Sample")
        return(alphada)
    }
    da <- alphada@alpha %>% as_tibble(rownames="Sample")
    if (action=="add"){
        samplevar <- .data %>% attr("samplevar")
        assaysvar <- .data %>% attr("assaysvar")
        othernm <- colnames(.data)[!colnames(.data) %in% c("OTU", "Sample", assaysvar, samplevar)]
        .data %<>% left_join(da, by="Sample", suffix=c("", ".y")) %>%
                   select(c("OTU", "Sample", assaysvar, samplevar, colnames(da), othernm))
        return(.data)
    }else if (action=="only"){
        return(da)
    }
}

#' @rdname mp_cal_alpha-methods
#' @aliases mp_cal_alpha,tbl_mpse
#' @exportMethod mp_cal_alpha
setMethod("mp_cal_alpha", signature(.data="tbl_mpse"), .internal_mp_cal_alpha)


#' @rdname mp_cal_alpha-methods
#' @aliases mp_cal_alpha,grouped_df_mpse
#' @exportMethod mp_cal_alpha
setMethod("mp_cal_alpha", signature(.data="grouped_df_mpse"), .internal_mp_cal_alpha)

valid_rare <- function(.data, ...){
    UseMethod("valid_rare")
}

valid_rare.MPSE <- function(.data, .abundance=NULL){
    assaysvar <- names(.data@assays)
    #.abundance <- ifelse(is.null(.abundance), "RareAbundance", .abundance)
    .abundance <- rlang::as_name(.abundance)
    if (.abundance %in% assaysvar){
        flag <- SummarizedExperiment::assays(.data)@listData[[.abundance]] %>% colSums %>% var == 0
        return(flag)
        #return(TRUE)
    }else{
        return(FALSE)
    }
}

.internal_valid_rare <- function(.data, .abundance=NULL){
    assaysvar <- attr(.data, "assaysvar")
    .abundance <- rlang::as_name(.abundance)
    #.abundance <- ifelse(is.null(.abundance), "RareAbundance", .abundance)
    if (.abundance %in% assaysvar){
        flag <- .data %>%
            group_by(!!as.symbol("Sample")) %>%
            dplyr::summarize(Total=sum(!!as.symbol(.abundance))) %>%
            select("Total") %>% unlist() %>% var == 0
        return(flag)
        #return(TRUE)
    }else{
        return(FALSE)
    }
}

valid_rare.tbl_mpse <- function(.data, .abundance=NULL){
    .internal_valid_rare(.data=.data, .abundance=.abundance)
}

valid_rare.grouped_df_mpse <- function(.data, .abundance=NULL){
    .internal_valid_rare(.data=.data, .abundance=.abundance)
}


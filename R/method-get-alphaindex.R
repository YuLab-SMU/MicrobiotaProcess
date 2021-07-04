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
setGeneric("get_alphaindex", function(obj, ...){standardGeneric("get_alphaindex")})

#' @aliases get_alphaindex,matrix
#' @rdname get_alphaindex
#' @importFrom vegan estimateR diversity specnumber
#' @export
setMethod("get_alphaindex", "matrix", function(obj, mindepth, sampleda, force=FALSE, ...){
    if (!identical(all.equal(obj, round(obj)),TRUE)){
        stop("the data should be integer (counts)!")
    }
    if (missing(mindepth) || is.null(mindepth)){
        mindepth <- min(rowSums(obj))
    }
    if (obj %>% rowSums %>% var != 0 && !force){
        obj <- vegan::rrarefy(obj, mindepth)
    }
    Chao <- estimateR(obj)
    Shannon <- diversity(obj)
    Simpson <- diversity(obj, index="simpson")
    J <- Shannon/log(specnumber(obj))
    alpha <- data.frame(Observe=Chao[1,],
                        Chao1=Chao[2,],
                        ACE=Chao[4,],
                        Shannon,
                        Simpson,
                        J)
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
        writeLines(c("The rarefied abundance of species might not be provided",
                    "rarefaction of all observations is performed automatically."))
        .data <- mp_rrarefy(.data=.data, ...)
    }
    
    alphada <- SummarizedExperiment::assays(.data)@listData[[rlang::as_name(.abundance)]] %>% 
               t() %>% 
               get_alphaindex(force=TRUE)
    if (action=="get"){
        alphada@sampleda <- .data@colData %>% as.data.frame()
        return(alphada)
    }
    da <- alphada@alpha %>% as_tibble(rownames="Sample")
    da <- .data@colData %>%
          as_tibble(rownames="Sample") %>%
          left_join(da, by="Sample") %>% 
          column_to_rownames(var="Sample")
    if (action=="add"){
        .data@colData <- S4Vectors::DataFrame(da)
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
        message("The rarefied abundance of species might not be provided. Rarefaction of all observations is performed automatically.")
        .data <- mp_rrarefy(.data=.data, ...)
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
        .data %<>% left_join(da, by="Sample") %>%
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


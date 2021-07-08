#' @title get the data of specified taxonomy
#' @param obj phyloseq, phyloseq class or data.frame
#' the shape of data.frame (nrow sample * column feature
#' taxa_are_rows set FALSE, nrow feature * ncol sample, 
#' taxa_are_rows set TRUE).
#' @param taxa_are_rows logical, if the column of data.frame
#' are features, it should be set FALSE.
#' @param taxda data.frame, the classifies of feature contained 
#' in obj(data.frame).
#' @param taxlevel character, the column names of taxda that you want to get.
#' when the input is phyloseq class, you can use 1 to 7.
#' @param sampleda data.frame, the sample information.
#' @param type character, the type of datasets, default is "species", 
#' if the dataset is not about species, such as dataset of kegg function, 
#' you should set it to "others". 
#' @param ..., additional parameters.
#' @return phyloseq class contained tax data.frame and sample information.
#' @author Shuangbin Xu
#' @export
#' @rdname get_taxadf
#' @examples
#' library(ggplot2)
#' data(test_otu_data)
#' phytax <- get_taxadf(test_otu_data, taxlevel=2)
#' phytax
#' head(phyloseq::otu_table(phytax))
#' phybar <- ggbartax(phytax) + 
#'          xlab(NULL) + ylab("relative abundance (%)")
setGeneric("get_taxadf", function(obj, ...)standardGeneric("get_taxadf"))

#' @aliases get_taxadf,phyloseq
#' @rdname get_taxadf
#' @export
setMethod("get_taxadf", "phyloseq", function(obj, taxlevel=2, type="species",...){
    if (is.null(obj@tax_table)){
    	stop("The tax table is empty!")
    }else{
        taxdf <- obj@tax_table
        
        if (!"fillNAtax" %in% names(attributes(taxdf))){
            taxdf <- fillNAtax(taxdf, type)
        }
    }
    otuda <- checkotu(obj)
    sampleda <- get_sample(obj)
    taxanames <- colnames(obj@tax_table)
    if (inherits(taxlevel, 'numeric')){taxlevel <- taxanames[taxlevel]}
    if (inherits(taxlevel, 'character')){
        if (!taxlevel %in% taxanames){
            stop("the taxlevel should be among the values of rank_names(phyloseq)")
        }else{
            taxlevel <- taxanames[match(taxlevel, taxanames)]
        }
    }
    taxdf <- get_taxadf(obj=otuda, 
                        taxda=taxdf, 
                        taxlevel=taxlevel,
                        sampleda=sampleda,
                        taxa_are_rows=FALSE,
                        type=type, 
                        ...)
    return(taxdf)
})

#' @aliases get_taxadf,data.frame
#' @rdname get_taxadf
#' @export 
setMethod("get_taxadf", "data.frame", 
          function(obj, taxda, 
                   taxa_are_rows,
                   taxlevel,
                   sampleda=NULL, 
                   type="species", ...){
    if (!taxa_are_rows){
        obj <- data.frame(t(obj), check.names=FALSE)
    }
    if(!is.null(sampleda) && !inherits(sampleda, "sample_data")){
        sampleda <- phyloseq::sample_data(sampleda)
    }
    if (!"fillNAtax" %in% names(attributes(taxda))){
        taxda <- fillNAtax(taxda, type=type)
    }
    if (inherits(taxlevel, "numeric")){taxlevel <- colnames(taxda)[taxlevel]}
    tmptax <- taxda[, match(taxlevel, colnames(taxda)), drop=FALSE]
    tmptaxda <- taxda[, seq(from=1, to=match(taxlevel, colnames(taxda))), drop=FALSE]
    tmptaxda <- tmptaxda[!duplicated(tmptaxda),,drop=FALSE]
    rownames(tmptaxda) <- as.vector(tmptaxda[, match(taxlevel, colnames(tmptaxda))])
    taxdf <- phyloseq::otu_table(get_count(data=obj, 
                                 featurelist=tmptax), 
                       taxa_are_rows=TRUE)
    taxdf <- phyloseq::phyloseq(taxdf, sampleda, phyloseq::tax_table(as.matrix(tmptaxda)))
    attr(taxdf@tax_table, "fillNAtax") <- TRUE
    return(taxdf)
})


#' Calculate the (relative) abundance of each taxonomy class for each sample or group.
#' @rdname mp_cal_abundance-methods
#' @param .data MPSE or tbl_mpse object
#' @param .abundance the name of otu abundance to be calculated
#' @param .group the name of group to be calculated.
#' @param relative logical whether calculate the relative abundance.
#' @param action character, "add" joins the new information to the taxatree and 
#' otutree if they exists (default). In addition, All taxonomy class will be added 
#' the taxatree, and OTU (tip) information will be added to the otutree."only" return 
#' a non-redundant tibble with the just new information. "get" return 'taxatree' slot
#' which is a treedata object.
#' @param force logical whether calculate the relative abundance forcibly when the abundance 
#' is not be rarefied, default is FALSE.
#' @param ... additional parameters.
#' @return update object or tibble according the 'action'
#' @export
setGeneric("mp_cal_abundance", 
           function(.data, 
                    .abundance = NULL, 
                    .group = NULL, 
                    relative = TRUE, 
                    action = "add",
                    force = FALSE, 
                    ...){
    standardGeneric("mp_cal_abundance")
})

#' @rdname mp_cal_abundance-methods
#' @aliases mp_cal_abundance,MPSE 
#' @importFrom dplyr across
#' @exportMethod mp_cal_abundance
setMethod("mp_cal_abundance", signature(.data="MPSE"), 
          function(.data, 
                   .abundance = NULL, 
                   .group = NULL, 
                   relative = TRUE, 
                   action = "add", 
                   force = FALSE,
                   ...){

    action %<>% match.arg(c("add", "get", "only"))
          
    .abundance <- rlang::enquo(.abundance)
    .group <- rlang::enquo(.group)

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
    assaysvar <- .data %>% SummarizedExperiment::assayNames()
    xx <- SummarizedExperiment::assays(.data)@listData

    da <- xx[[rlang::as_name(.abundance)]] %>%
          tibble::as_tibble(rownames="OTU") %>%
          tidyr::pivot_longer(!as.symbol("OTU"), names_to="Sample", values_to=rlang::as_name(.abundance))

    sampleda <- .data@colData %>%
                avoid_conflict_names() %>%
                tibble::as_tibble(rownames="Sample")

    if (ncol(sampleda)>1){
        da %<>% left_join(sampleda, by="Sample")
    }

    otumeta <-
        SummarizedExperiment::rowData(.data) %>%
        avoid_conflict_names() %>%
        tibble::as_tibble(rownames="OTU")
    
    if (ncol(otumeta) > 1){
        da %<>% dplyr::left_join(otumeta, by="OTU")
    }
    
    if (!is.null(.data@taxatree)){
        taxada <- taxatree_to_tb(.data@taxatree) %>%
                  tibble::as_tibble(rownames="OTU")
        taxada <- taxada[ ,!colnames(taxada) %in% 
                           c(colnames(.data@taxatree@data), 
                             colnames(.data@taxatree@extraInfo)),
                           drop=FALSE]
        da %<>% dplyr::left_join(taxada, by="OTU")
        taxavar <- colnames(taxada)
    }else{
        taxavar <- "OTU"
    }
    
    if (!rlang::quo_is_null(.group)){
        da1 <- lapply(taxavar, function(x) 
                               .internal_cal_feature_abun(da=da, 
                                         .abundance=.abundance, 
                                         feature=x, 
                                         byID=.group,
                                         relative=relative))
    }else{
        da1 <- lapply(taxavar, function(x)
                      .internal_cal_feature_abun(da=da,
                                         .abundance=.abundance,
                                         feature=x,
                                         byID=as.symbol("Sample"),
                                         relative=relative))
    }

    if (action %in% c("add", "get")){
        if (rlang::quo_is_null(.group) && relative){
            newRelabun <- paste0(c("Rel", rlang::as_name(.abundance), "BySample"), collapse="")
            otuRelabun <- da1[[1]] %>% 
                          select(-!!.abundance) %>%
                          tidyr::pivot_wider(id_cols="OTU", names_from="Sample", values_from=as.symbol(newRelabun)) %>%
                          tibble::column_to_rownames(var="OTU")
            SummarizedExperiment::assays(.data)@listData <- c(xx, list(otuRelabun)) %>% 
                setNames(c(assaysvar, newRelabun))
        }
        
        da1 %<>% 
             dplyr::bind_rows() %>% 
             nest_internal() %>% 
             rename(label="OTU")

        if (!is.null(.data@taxatree)){
            extranm <- setdiff(colnames(da1), c(colnames(.data@taxatree@data), colnames(.data@taxatree@extraInfo)))
            da1 %<>% dplyr::select(extranm)
            .data@taxatree %<>% treeio::full_join(da1, by="label")
        }

        if (!is.null(.data@otutree)){
            da1 %<>% dplyr::filter(.data$label %in% .data@otutree@phylo$tip.label)
            extranm <- setdiff(colnames(da1), c(colnames(.data@otutree@data), colnames(.data@otutree@extraInfo)))
            da1 %<>% dplyr::select(extranm)
            .data@otutree %<>% treeio::full_join(da1, by="label")
        }
        
        if (action=="add"){
            return(.data)
        }else{
            if (is.null(.data@taxatree)){
                message("The taxatree of the MPSE object is empty!")
            }
            return(.data@taxatree)
        }

    }else if(action=="only"){
        da1 %<>% 
            setNames(taxavar) %>%
            dplyr::bind_rows(.id="TaxaClass") %>% 
            dplyr::rename(AllTaxa="OTU")
        
        if (ncol(sampleda)>1){
            da1 %<>% dplyr::left_join(sampleda, by="Sample")
        }

        return(da1)
    }
    
})

.internal_cal_feature_abun <- function(da, .abundance, feature, byID, relative){
    Totalnm <- paste0("TotalNumsBy", rlang::as_name(byID))
    if(rlang::as_name(byID)=="Sample"){
        newabun <- rlang::as_name(.abundance)
        sampleind <- NULL
    }else{
        newabun <- paste0(c(rlang::as_name(.abundance), "By", rlang::as_name(byID)), collapse="")
        sampleind <- as.symbol("Sample")
    }

    da %<>%
        dplyr::group_by(!!byID) %>%
        dplyr::mutate(across(!!.abundance, sum, .names=Totalnm)) %>%
        dplyr::group_by(across(c(!!as.symbol(feature), !!byID))) %>%
        dplyr::mutate(across(!!.abundance, sum, .names=newabun))

    if (relative){
        newRelabun <- paste0(c("Rel", rlang::as_name(.abundance), "By", rlang::as_name(byID)), collapse="")
        da %<>%
            dplyr::mutate(across(!!as.symbol(newabun), ~ .x/!!as.symbol(Totalnm) * 100, .names=newRelabun))
        if(is.null(sampleind)){
            da %<>% select(c(as.symbol(feature), !!byID, as.symbol(newabun), as.symbol(newRelabun)))
        }else{
            da %<>% select(c(as.symbol(feature), !!sampleind, !!byID, as.symbol(newabun), as.symbol(newRelabun)))
        }
        da %<>% 
            ungroup() %>%
            distinct() 

    }else{
        if (is.null(sampleind)){
            da %<>%
            select(c(as.symbol(feature), !!byID, as.symbol(newabun)))
        }else{
            da %<>%
            select(c(as.symbol(feature), !!sampleind, !!byID, as.symbol(newabun)))
        }
        da %<>%
            ungroup() %>%
            distinct()
    }
    colnames(da)[1] <- "OTU"
    return(da)
}

.internal_mp_cal_abundance <- function(.data, .abundance=NULL, .group=NULL, relative=TRUE, action="add", force=FALSE, ...){

    action %<>% match.arg(c("add", "get", "only"))

    .abundance <- rlang::enquo(.abundance)
    .group <- rlang::enquo(.group)
    
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

    assaysvar <- .data %>% attr("assaysvar")
    taxavar <- .data %>% attr("taxavar")

    if (!is.null(taxavar)){
        taxavar <- c("OTU", taxavar)
    }else{
        taxavar <- "OTU"    
    }
    
    if (!rlang::quo_is_null(.group)){
        da1 <- lapply(taxavar, function(x)
                               .internal_cal_feature_abun(da=.data,
                                         .abundance=.abundance,
                                         feature=x,
                                         byID=.group,
                                         relative=relative))
    }else{
        da1 <- lapply(taxavar, function(x)
                      .internal_cal_feature_abun(da=.data,
                                         .abundance=.abundance,
                                         feature=x,
                                         byID=as.symbol("Sample"),
                                         relative=relative))
    }
    
    if (action %in% c("add", "get")){
        if (rlang::quo_is_null(.group) && relative){
            newRelabun <- paste0(c("Rel", rlang::as_name(.abundance), "BySample"), collapse="")
            dx1 <- da1[[1]] %>% select(c("OTU", "Sample", as.symbol(newRelabun)))
            othernm <- colnames(.data)[!colnames(.data) %in% c("OTU", "Sample", assaysvar)]
            attr(.data, "assaysvar") <- c(assaysvar, newRelabun)
            .data %<>% 
                  left_join(dx1, by=c("OTU", "Sample")) %>%
                  select(c("OTU", "Sample", assaysvar, newRelabun, othernm))
        }

        da1 %<>%
             dplyr::bind_rows() %>%
             nest_internal() %>%
             rename(label="OTU")

        taxatree <- .data %>% attr("taxatree")

        if (!is.null(taxatree)){
            dat1 <- da1 %>% select(setdiff(colnames(da1), c(colnames(taxatree@data), colnames(taxatree@extraInfo))))
            attr(.data, "taxatree") <- taxatree %>% treeio::full_join(dat1, by="label")
        }
     
        otutree <- .data %>% attr("otutree")

        if (!is.null(otutree)){
            dat2 <- da1 %>% dplyr::filter(.data$label %in% .data@otutree@phylo$tip.label) %>%
                     select(setdiff(colnames(da1), c(colnames(otutree@data), colnames(otutree@extraInfo))))
            attr(.data, "otutree") <- otutree %>% treeio::full_join(dat2, by="label")
        }
        
        if (action=="get"){
            if (is.null(taxatree)){
                message("The taxatree of the MPSE object is empty")
            }

            return(attr(.data, "taxatree"))

        }else{
            return(.data)
        }

    }else if(action=="only"){
        da1 %<>%
            setNames(taxavar) %>%
            dplyr::bind_rows(.id="TaxaClass") %>%
            dplyr::rename(AllTaxa="OTU")
        
        samplevar <- .data %>% attr("samplevar")

        if (length(samplevar)>1){
            sampleda <- .data %>% ungroup() %>% select(samplevar)
            da1 %<>% dplyr::left_join(sampleda, by="Sample")
        }

        return(da1)
    }
}

#' @rdname mp_cal_abundance-methods
#' @aliases mp_cal_abundance,tbl_mpse
#' @exportMethod mp_cal_abundance
setMethod("mp_cal_abundance", signature(.data = "tbl_mpse"), .internal_mp_cal_abundance)

#' @rdname mp_cal_abundance-methods
#' @aliases mp_cal_abundance,grouped_df_mpse
#' @exportMethod mp_cal_abundance
setMethod("mp_cal_abundance", signature(.data = "grouped_df_mpse"), .internal_mp_cal_abundance)

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
#' \dontrun{
#' library(ggplot2)
#' data(test_otu_data)
#' phytax <- get_taxadf(test_otu_data, taxlevel=2)
#' phytax
#' head(phyloseq::otu_table(phytax))
#' phybar <- ggbartax(phytax) + 
#'          xlab(NULL) + ylab("relative abundance (%)")
#' }
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
#' @seealso [mp_plot_abundance()] and [mp_extract_abundance()]
#' @author Shuangbin Xu
#' @export
#' @examples
#' data(mouse.time.mpse)
#' mouse.time.mpse %<>%
#'   mp_rrarefy() 
#' mouse.time.mpse
#' mouse.time.mpse %<>%
#'   mp_cal_abundance(.abundance=RareAbundance, action="add") %>% 
#'   mp_cal_abundance(.abundance=RareAbundance, .group=time, action="add") 
#' mouse.time.mpse
#' p1 <- mouse.time.mpse %>% 
#'       mp_plot_abundance(.abundance=RelRareAbundanceBySample, 
#'                         .group=time, taxa.class="Phylum", topn=20)
#' p2 <- mouse.time.mpse %>% 
#'       mp_plot_abundance(.abundance = RareAbundance, 
#'                         .group = time, 
#'                         taxa.class = Phylum, 
#'                         topn = 20, 
#'                         relative = FALSE, 
#'                         force = TRUE)
#' p1 / p2
#' # Or you can also extract the result and visulize it with ggplot2 and ggplot2-extension
#' \dontrun{
#' tbl <- mouse.time.mpse %>%
#'        mp_extract_abundance(taxa.class="Class", topn=10)
#' tbl
#' library(ggplot2)
#' library(ggalluvial)
#' library(dplyr)
#' tbl %<>%
#'   tidyr::unnest(cols=RareAbundanceBySample) 
#' tbl
#' p <- ggplot(data=tbl,
#'             mapping=aes(x=Sample, 
#'                         y=RelRareAbundanceBySample, 
#'                         alluvium=label,
#'                         fill=label)
#'      ) + 
#'      geom_flow(stat="alluvium", lode.guidance = "frontback", color = "darkgray") +
#'      geom_stratum(stat="alluvium") +
#'      labs(x=NULL, y="Relative Abundance (%)") +
#'      scale_fill_brewer(name="Class", type = "qual", palette = "Paired") +
#'      facet_grid(cols=vars(time), scales="free_x", space="free") +
#'      theme(axis.text.x=element_text(angle=-45, hjust=0))
#' p
#' }
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
        trash <- try(silent = TRUE,
                     expr = {
                         .data <- mp_rrarefy(.data = .data, ...)
                     }
                 )
        if (inherits(trash, "try-error")){
            stop_wrap("The 'Abundance' column cannot be rarefied, please check whether it is integer (count).
                       Or you can set 'force=TRUE' to calculate the (relative) abundance without rarefaction.
                      ")
        }

        message_wrap("The rarefied abundance of species might not be provided. Rarefaction of all
                      observations is performed automatically using 'Abundance' column. If you still 
                      want to calculate the (relative) 'abundance' with the specified '.abundance',
                      you can set 'force=TRUE'. ")
        .abundance <- as.symbol("RareAbundance")           
    }
    assaysvar <- .data %>% SummarizedExperiment::assayNames()
    xx <- SummarizedExperiment::assays(.data)@listData

    da <- xx[[rlang::as_name(.abundance)]] %>%
          tibble::as_tibble(rownames="OTU") %>%
          tidyr::pivot_longer(!as.symbol("OTU"), names_to="Sample", values_to=rlang::as_name(.abundance)) %>%
          dtplyr::lazy_dt()

    sampleda <- .data %>% mp_extract_sample()
    if (ncol(sampleda)==1){
        sampleda %<>% dplyr::mutate(.DTPLYREXTRA=0)
    }

    da %<>% left_join(sampleda, by="Sample", suffix=c("", ".y"))
    if (".DTPLYREXTRA" %in% colnames(sampleda)){
        sampleda %<>% select(-".DTPLYREXTRA")
    }
    otumeta <-
        SummarizedExperiment::rowData(.data) %>%
        avoid_conflict_names() %>%
        tibble::as_tibble(rownames="OTU")
    
    if (ncol(otumeta) > 1){
        da %<>% dplyr::left_join(otumeta, by="OTU", suffix=c("", ".y"))
    }
    
    if (!is.null(.data@taxatree)){
        #taxada <- taxatree_to_tb(.data@taxatree) %>%
        #          tibble::as_tibble(rownames="OTU")
        #taxada <- taxada[ ,!colnames(taxada) %in% 
        #                   c(colnames(.data@taxatree@data), 
        #                     colnames(.data@taxatree@extraInfo)),
        #                   drop=FALSE]
        taxada <- .data %>% mp_extract_taxonomy()
        da %<>% dplyr::left_join(taxada, by="OTU", suffix=c("", ".y"))
        taxavar <- colnames(taxada)
    }else{
        taxavar <- "OTU"
    }
    if (!rlang::quo_is_null(.group)){
        da1 <- lapply(rlang::syms(taxavar), 
                      .internal_cal_feature_abun,
                                         da = da, 
                                         .abundance = .abundance, 
                                         byID = .group,
                                         relative = relative,
                                         sampleda = NULL
               )
    }else{
        sampledat <- sampleda[, !vapply(sampleda, function(x)is.list(x)||is.numeric(x), logical(1))]
        da1 <- lapply(rlang::syms(taxavar),
                      .internal_cal_feature_abun,
                                         da = da,
                                         .abundance = .abundance,
                                         byID = as.symbol("Sample"),
                                         relative = relative,
                                         sampleda = sampledat
               )
    }

    if (rlang::quo_is_null(.group) && relative){
        newRelabun <- paste0("Rel", rlang::as_name(.abundance), "BySample")
        otuRelabun <- da1[[1]] %>% 
				      tidyr::unnest(cols=paste0(rlang::as_name(.abundance),"BySample")) %>%
                      select(-!!.abundance) %>%
                      tidyr::pivot_wider(id_cols="OTU", names_from="Sample", values_from=as.symbol(newRelabun)) %>%
                      tibble::column_to_rownames(var="OTU")
        SummarizedExperiment::assays(.data)@listData <- c(xx, list(otuRelabun)) %>% 
            setNames(c(assaysvar, newRelabun))
    }
    
    	
    da1 %<>% 
         dplyr::bind_rows() %>%
         #nest_internal() 
         rename(label="OTU")

    if (!is.null(.data@taxatree)){
        extranm <- intersect(colnames(da1), c(colnames(.data@taxatree@data), colnames(.data@taxatree@extraInfo)))
        .data@taxatree %<>% tidytree::select(-c(extranm), keep.td=TRUE) %>%  
            treeio::full_join(da1, by="label")
    }else{
        .data %<>% left_join(da1, by=c("OTU"="label"))
    }
    otutree <- .data %>% mp_extract_tree(type="otutree")
    if (!is.null(otutree)){
        da1 %<>% dplyr::filter(!!as.symbol("label") %in% otutree@phylo$tip.label)
        extranm <- intersect(colnames(da1), c(colnames(otutree@data), colnames(otutree@extraInfo)))
        otutree(.data) <- otutree %>% tidytree::select(-c(extranm), keep.td=TRUE) %>% treeio::full_join(da1, by="label")
    }
    
    if (action=="add"){
        return(.data)
    }else if (action=="get"){
        if (is.null(.data@taxatree)){
            message("The taxatree of the MPSE object is empty!")
        }
        return(.data@taxatree)
    }else if (action=="only"){
       if (is.null(.data@taxatree)){
           da1 %<>% 
               tidyr::unnest () %>% 
               suppressWarnings()
           
           if (ncol(sampleda)>1){
               sampleda %<>% dplyr::select(c("Sample", setdiff(colnames(sampleda), colnames(da1))))
               da1 %<>% dplyr::left_join(sampleda, by="Sample", suffix=c("", ".y")) %>% dplyr::distinct()
           }
       }else{
           da1 <- .data@taxatree %>% 
                   as_tibble() %>%
                   dplyr::select(-c("parent", "node", "nodeDepth")) %>%
                   dplyr::filter(.data$nodeClass != "Root")
       }
       return(da1)
    }
    
})

.internal_cal_feature_abun <- function(feature, da, .abundance, byID, relative, sampleda){
    Totalnm <- paste0("TotalNumsBy", rlang::as_name(byID))
    if(rlang::as_name(byID)=="Sample"){
        newabun <- rlang::as_name(.abundance)
        bygroup <- paste0(newabun, "BySample")
    }else{
        newabun <- paste0(rlang::as_name(.abundance), "By", rlang::as_name(byID))
        bygroup <- newabun
    }
    da %<>%
        dplyr::group_by(!!byID) %>%
        dplyr::mutate(across(!!.abundance, sum, .names=Totalnm)) %>%
        dplyr::group_by(!!feature, .add = TRUE) %>%
        dplyr::mutate(across(!!.abundance, sum, .names=newabun))

    if (relative){
        newRelabun <- paste0(c("Rel", rlang::as_name(.abundance), "By", rlang::as_name(byID)), collapse="")
        da %<>%
            dplyr::mutate(across(!!as.symbol(newabun), ~ .x/!!as.symbol(Totalnm) * 100, .names=newRelabun)) %>%
            select(!!feature, !!byID, !!as.symbol(newabun), !!as.symbol(newRelabun))
    }else{
        da %<>% select(!!feature, !!byID, !!as.symbol(newabun)) 
    }
    da %<>% dplyr::ungroup() %>% dplyr::distinct() %>% dplyr::rename(OTU=1)

    if (!is.null(sampleda) && ncol(sampleda)>1){
        da <- da %>% dplyr::left_join(sampleda, by="Sample", suffix=c("", ".y"))
    }

    da %<>% as_tibble() %>% tidyr::nest(!!bygroup:=colnames(.)[ colnames(.) !="OTU"])
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
        trash <- try(silent = TRUE,
                     expr = {
                         .data <- mp_rrarefy(.data = .data, ...)
                     }
                 )
        if (inherits(trash, "try-error")){
            stop_wrap("The 'Abundance' column cannot be rarefied, please check whether it is integer (count).
                       Or you can set 'force=TRUE' to calculate the (relative) abundance without rarefaction.
                      ")
        }

        message_wrap("The rarefied abundance of species might not be provided. Rarefaction of all
                      observations is performed automatically using 'Abundance' column. If you still
                      want to calculate the (relative) 'abundance' with the specified '.abundance',
                      you can set 'force=TRUE'. ")
        .abundance <- as.symbol("RareAbundance")        
    }

    assaysvar <- .data %>% attr("assaysvar")
    taxavar <- .data %>% attr("taxavar")
    sampleda <- .data %>% mp_extract_sample()
    if (!is.null(taxavar)){
        taxavar <- c("OTU", taxavar)
    }else{
        taxavar <- "OTU"    
    }
    
    if (!rlang::quo_is_null(.group)){
        da1 <- lapply(rlang::syms(taxavar), 
                      .internal_cal_feature_abun,
                                         da = dtplyr::lazy_dt(.data),
                                         .abundance = .abundance,
                                         byID = .group,
                                         relative = relative,
                                         sampleda = NULL
                      )
    }else{
        sampledat <- sampleda[, !vapply(sampleda, function(x)is.list(x)||is.numeric(x), logical(1))]
        da1 <- lapply(rlang::syms(taxavar), 
                      .internal_cal_feature_abun,
                                         da = dtplyr::lazy_dt(.data),
                                         .abundance = .abundance,
                                         byID = as.symbol("Sample"),
                                         relative = relative,
                                         sampleda = sampledat
                                         )
    }
    
    if (rlang::quo_is_null(.group) && relative){
        newRelabun <- paste0("Rel", rlang::as_name(.abundance), "BySample")
        dx1 <- da1[[1]] %>% 
               tidyr::unnest(cols=paste0(rlang::as_name(.abundance),"BySample")) %>% 
               select(c("OTU", "Sample", as.symbol(newRelabun)))
        othernm <- colnames(.data)[!colnames(.data) %in% c("OTU", "Sample", assaysvar)]
        attr(.data, "assaysvar") <- c(assaysvar, newRelabun)
        .data %<>% 
              left_join(dx1, by=c("OTU", "Sample"), suffix=c("", ".y")) %>%
              select(c("OTU", "Sample", assaysvar, newRelabun, othernm))
    }
    
    da1 %<>%
         dplyr::bind_rows() %>%
         drop_class("tbl_mpse") %>%
         #nest_internal() %>%
         dplyr::rename(label="OTU")

    taxatree <- .data %>% attr("taxatree")

    if (!is.null(taxatree)){
        extranm <- intersect(colnames(da1), c(colnames(taxatree@data), colnames(taxatree@extraInfo)))
        attr(.data, "taxatree") <- taxatree %>% 
                                   tidytree::select(-c(extranm), keep.td=TRUE) %>% 
                                   treeio::full_join(da1, by="label")
    }else{
        .data %<>% left_join(da1, by=c("OTU"="label"))
    }
    
    otutree <- .data %>% attr("otutree")

    if (!is.null(otutree)){
        dat2 <- da1 %>% dplyr::filter(!!as.symbol("label") %in% otutree@phylo$tip.label)
        extranm <- intersect(colnames(da1), c(colnames(otutree@data), colnames(otutree@extraInfo)))
        attr(.data, "otutree") <- otutree %>% 
            tidytree::select(-c(extranm), keep.td=TRUE) %>%
            treeio::full_join(dat2, by="label")
    }
    
    if (action=="get"){
        if (is.null(taxatree)){
            message("The taxatree of the object is empty")
        }

        return(attr(.data, "taxatree"))

    }else if (action=="add"){
        return(.data)
    }else if(action=="only"){
        if (is.null(taxatree)){
            da1 %<>%
                  tidyr::unnest() %>%
                  suppressWarnings() 

            samplevar <- .data %>% attr("samplevar")

            if (length(samplevar)>1){
                sampleda <- .data %>% 
                            ungroup() %>% 
                            select(c("Sample", setdiff(samplevar, colnames(da1))))
                da1 %<>% dplyr::left_join(sampleda, by="Sample", suffix=c("", ".y"))
            }
        }else{
            da1 <- .data %>% 
                   attr("taxatree") %>% 
                   as_tibble %>%
                   dplyr::select(-c("parent", "node", "nodeDepth")) %>%
                   dplyr::filter(.data$nodeClass != "Root")
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

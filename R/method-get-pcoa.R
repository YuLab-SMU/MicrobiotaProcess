#' @title performs principal coordinate analysis (PCoA)
#' @param data data.frame, numeric data.frame nrow sample * ncol features.
#' @param obj phyloseq, the phyloseq class or dist class.
#' @param sampleda data.frame, nrow sample * ncol factor, default is NULL.
#' @param distmethod character, the method of distance, 
#' see also \code{\link[phyloseq]{distance}}
#' @param taxa_are_rows logical, if feature of data is column, 
#' it should be set FALSE.
#' @param tree phylo, the phylo class, default is NULL, 
#' when use unifrac method, it should be required.
#' @param method character, the standardization method for 
#' community ecologists, default is hellinger, if the data 
#' has be normlized, it shound be set NULL.
#' @param ..., additional parameter, see also
#' \code{\link[MicrobiotaProcess]{get_dist}}.
#' @return pcasample object, contained prcomp or 
#' pcoa and sampleda (data.frame).
#' @author Shuangbin Xu
#' @export
#' @examples
#' \dontrun{
#'     library(phyloseq)
#'     data(GlobalPatterns)
#'     subGlobal <- subset_samples(GlobalPatterns, 
#'                   SampleType %in% c("Feces", "Mock", "Ocean", "Skin"))
#'     pcoares <- get_pcoa(subGlobal, 
#'                        distmethod="euclidean",
#'                        method="hellinger")
#'     pcoaplot <- ggordpoint(pcoares, biplot=FALSE,
#'                             speciesannot=FALSE,
#'                             factorNames=c("SampleType"), 
#'                             ellipse=FALSE)
#' }
get_pcoa <- function(obj, ...){
    UseMethod("get_pcoa")
}

#' @method get_pcoa data.frame
#' @rdname get_pcoa
#' @export
get_pcoa.data.frame <- function(obj, 
    distmethod="euclidean", 
    taxa_are_rows=FALSE,
    sampleda=NULL, 
    tree=NULL,
    #type="sample",
    method="hellinger",
    ...){
    tmpdist <- get_dist.data.frame(obj, 
                                distmethod=distmethod,
                                taxa_are_rows=taxa_are_rows,
                                sampleda=sampleda,
                                tree=tree,
                                #type=type,
                                method=method,
                                ...)
    data <- attr(tmpdist, "originalD")
    pcoares <- get_pcoa.dist(tmpdist, 
                             distmethod=distmethod, 
                             data=data,
                             sampleda=sampleda)
    return(pcoares)
}

#' @method get_pcoa dist
#' @importFrom ape pcoa
#' @importFrom stats cov
#' @rdname get_pcoa
#' @export
get_pcoa.dist <- function(obj, distmethod, data=NULL, sampleda=NULL, method="hellinger", ...){
    if (missing(distmethod)){
        if (!is.null(attr(obj, "distmethod"))){
            distmethod <- attr(obj, "distmethod")
        }else{
            distmethod <- NULL
        }
    }
    pcoares <- pcoa(obj, ...)
    attr(pcoares, "distmethod") <- distmethod
    if (!is.null(data)){
        n <- nrow(data)
        points.stand <- scale(pcoares$vectors)
        S <- cov(data, points.stand)
        tmpEig <- pcoares$values$Eigenvalues
        tmpEig <- tmpEig[seq_len(dim(S)[2])]
        diagtmp <- diag((tmpEig/(n-1))^(-0.5))
        U <- S %*% diagtmp
        colnames(U) <- colnames(pcoares$vectors)
        attr(pcoares, "varcorr") <- U
    }
    res <- new("pcasample", 
               pca=pcoares,
               sampleda=sampleda)
    return(res)
}

#' @method get_pcoa phyloseq
#' @rdname get_pcoa
#' @export
get_pcoa.phyloseq <- function(obj,distmethod="euclidean",...){
    sampleda <- checksample(obj)
    tmpdist <- get_dist.phyloseq(obj, distmethod=distmethod,...)
    otuda <- attr(tmpdist, "originalD")
    pcoares <- get_pcoa.dist(tmpdist, 
                             distmethod=distmethod, 
                             data=otuda, 
                             sampleda=sampleda)
    return(pcoares)
}

#' @method get_coord pcoa
#' @rdname get_coord
#' @export
get_coord.pcoa <- function(obj, pc){
    coord <- obj$vector[,pc]
    vp <- obj$values$Relative_eig
    vp <- vp*100/sum(vp)
    tmpvp1 <- round(vp[pc[1]], 2)
    tmpvp2 <- round(vp[pc[2]], 2)
    xlab_text <- paste0("PCoA", pc[1], "(", tmpvp1, "%)")
    ylab_text <- paste0("PCoA", pc[2], "(", tmpvp2, "%)")
    title_text <- paste0("PCoA - PCoA",pc[1], " VS PCoA",pc[2], " (", attr(obj, "distmethod"), ")")
    ordplotclass <- new("ordplotClass",
                        coord=coord,
                        xlab=xlab_text,
                        ylab=ylab_text,
                        title=title_text
                        )
    return(ordplotclass)
}

#' @method get_varct pcoa
#' @rdname get_varct
#' @export
get_varct.pcoa <- function(obj,...){
    if (is.null(attr(obj,"varcorr"))){
        stop("The pcoa class have not `varcorr` attr")
    }
    else{
        varcorr <- attr(obj, "varcorr")
        varcorr2 <- varcorr^2
        componentcos <- apply(varcorr2, 2, sum)
        varcontrib <- data.frame(t(apply(varcorr2,
                                         1, 
                                         contribution,
                                         componentcos)),
                                         check.names=FALSE)
        res <- list(VarContribution=varcontrib,
                    VarCoordinates=varcorr)
        attr(res, "class") <- "VarContrib"
        return(res)
    }
}

#' Principal Coordinate Analysis with MPSE or tbl_mpse object
#' @rdname mp_cal_pcoa-methods
#' @param .data MPSE or tbl_mpse object
#' @param .abundance the name of abundance to be calculated.
#' @param distmethod character the method to calculate distance.
#' @param .dim integer The number of dimensions to be returned, default is 3.
#' @param action character "add" joins the pca result to the object, "only" return
#' a non-redundant tibble with the pca result. "get" return 'pcasample' object can
#' be visualized with 'ggordpoint'.
#' @param ... additional parameters see also 'mp_cal_dist'.
#' @return update object or tbl according to the action.
#' @export
#' @author Shuangbin Xu
#' @examples
#' data(mouse.time.mpse)
#' mouse.time.mpse %>% 
#'   mp_decostand(.abundance=Abundance) %>% 
#'   mp_cal_pcoa(.abundance=hellinger, distmethod="bray", .dim=2, action="only") -> tbl
#' tbl
#' x <- names(tbl)[grepl("PCo1 ", names(tbl))] %>% as.symbol()
#' y <- names(tbl)[grepl("PCo2 ", names(tbl))] %>% as.symbol()
#' library(ggplot2)
#' tbl %>% 
#'  ggplot(aes(x=!!x, y=!!y, color=time)) + 
#'  geom_point() +
#'  geom_vline(xintercept=0, color="grey20", linetype=2) + 
#'  geom_hline(yintercept=0, color="grey20", linetype=2) +
#'  theme_bw() +
#'  theme(panel.grid=element_blank())
setGeneric("mp_cal_pcoa", function(.data, .abundance, distmethod="bray", .dim=3, action="add", ...)standardGeneric("mp_cal_pcoa"))

#' @rdname mp_cal_pcoa-methods
#' @aliases mp_cal_pcoa,MPSE
#' @exportMethod mp_cal_pcoa
setMethod("mp_cal_pcoa", signature(.data="MPSE"), function(.data, .abundance, distmethod="bray", .dim=3, action="add", ...){
    
    action %<>% match.arg(c("add", "only", "get"))

    .abundance <- rlang::enquo(.abundance)
    
    if(! distmethod %in% colnames(.data@colData)){
        if (rlang::quo_is_missing(.abundance)){
            rlang::abort("The .abundance must be required, when the distmethod is not present in the object.")
        }else{
            .data %<>% mp_cal_dist(.abundance=!!.abundance, distmethod=distmethod, action="add", ...)
        }
    }
    
    distobj <- .data %>% mp_extract_dist(distmethod=distmethod)
    
    pcoa <- ape::pcoa(distobj)

    if (action=="get"){
        #sampleda <- .data %>%
        #            mp_extract_sample() %>%
        #            tibble::column_to_rownames(var="Sample")
        #res <- new("pcasample", pca=pcoa, sampleda=sampleda)
        return(pcoa)         
    }
    dat <- pcoa %>% tidydr(display=c("sites", "features"))

    da <- .data %>%
          mp_extract_sample() %>%
          dplyr::left_join(
                 dat[, seq_len(.dim+1)],
                 by=c("Sample"="sites")
          )

    if (action=="only"){
        da %<>%
             add_attr(dat %>% attr("features_tb"), name="features_tb") %>%
             add_internal_attr(object=pcoa, name="PCoA")
        return(da)
    }else if (action=="add"){
        .data@colData <- da %>%
                         tibble::column_to_rownames(var="Sample") %>%
                         S4Vectors::DataFrame(check.names=FALSE)
        .data %<>% add_internal_attr(object=pcoa, name="PCoA")
        return(.data)
    }
})


.internal_cal_pcoa <- function(.data, .abundance, distmethod="bray", .dim=3, action="add", ...){

    action %<>% match.arg(c("add", "only", "get"))

    .abundance <- rlang::enquo(.abundance)

    if(! distmethod %in% colnames(.data)){
        if (rlang::quo_is_missing(.abundance)){
            rlang::abort("The .abundance must be required, when the distmethod is not present in the object.")
        }else{
            .data %<>% mp_cal_dist(.abundance=!!.abundance, distmethod=distmethod, action="add", ...)
        }
    }

    distobj <- .data %>% mp_extract_dist(distmethod=distmethod)

    pcoa <- ape::pcoa(distobj)

    if (action=="get"){
        sampleda <- .data %>%
                    mp_extract_sample() %>%
                    tibble::column_to_rownames(var="Sample")
        res <- new("pcasample", pca=pcoa, sampleda=sampleda)
        return(res)
    }else if (action=="only"){
        da <- .data %>%
              mp_extract_sample() %>%
              dplyr::left_join(
                  pcoa$vectors[,seq_len(.dim)] %>%
                  as_tibble(rownames="Sample") %>%
                  setNames(c("Sample", paste0("PCoA", seq_len(.dim)))),
                 by="Sample"
              ) %>%
              add_internal_attr(object=pcoa, name="PCoA")
        return(da)
    }else if (action=="add"){
        .data %<>%
            dplyr::left_join(
                pcoa$vectors[,seq_len(.dim)] %>%
                as_tibble(rownames="Sample") %>%
                setNames(c("Sample", paste0("PCoA", seq_len(.dim)))),
                by="Sample"
            ) %>%
            add_internal_attr(object=pcoa, name="PCoA")

        return(.data)
    }
}

#' @rdname mp_cal_pcoa-methods
#' @aliases mp_cal_pcoa,tbl_mpse
#' @exportMethod mp_cal_pcoa
setMethod("mp_cal_pcoa", signature(.data="tbl_mpse"), .internal_cal_pcoa)

#' @rdname mp_cal_pcoa-methods
#' @aliases mp_cal_pcoa,grouped_df_mpse
#' @exportMethod mp_cal_pcoa
setMethod("mp_cal_pcoa", signature(.data="grouped_df_mpse"), .internal_cal_pcoa)

#' Nonmetric Multidimensional Scaling Analysis with MPSE or tbl_mpse object
#' @rdname mp_cal_nmds-methods
#' @param .data MPSE or tbl_mpse object
#' @param .abundance the name of abundance to be calculated.
#' @param distmethod character the method to calculate distance.
#' @param .dim integer The number of dimensions to be returned, default is 2.
#' @param action character "add" joins the NMDS result to the object, "only" return
#' a non-redundant tibble with the NMDS result. "get" return 'metaMDS' object can
#' be analyzed with related 'vegan' function.
#' @param seed a random seed to make this analysis reproducible, default is 123.
#' @param ... additional parameters see also 'mp_cal_dist'.
#' @return update object or tbl according to the action.
#' @author Shuangbin Xu
#' @export
#' @examples
#' data(mouse.time.mpse)
#' mouse.time.mpse %>%
#'   mp_decostand(.abundance=Abundance) %>%
#'   mp_cal_nmds(.abundance=hellinger, distmethod="bray", .dim=2, action="only") -> tbl
#' tbl
#' x <- names(tbl)[grepl("NMDS1", names(tbl))] %>% as.symbol()
#' y <- names(tbl)[grepl("NMDS2", names(tbl))] %>% as.symbol()
#' library(ggplot2)
#' tbl %>%
#'  ggplot(aes(x=!!x, y=!!y, color=time)) +
#'  geom_point() +
#'  geom_vline(xintercept=0, color="grey20", linetype=2) +
#'  geom_hline(yintercept=0, color="grey20", linetype=2) +
#'  theme_bw() +
#'  theme(panel.grid=element_blank())
setGeneric("mp_cal_nmds", function(.data, .abundance, distmethod="bray", .dim=2, action="add", seed=123, ...)standardGeneric("mp_cal_nmds"))

#' @rdname mp_cal_nmds-methods
#' @aliases mp_cal_nmds,MPSE
#' @exportMethod mp_cal_nmds
setMethod("mp_cal_nmds", signature(.data="MPSE"), function(.data, .abundance, distmethod="bray", .dim=2, action="add", seed=123, ...){

    action %<>% match.arg(c("add", "only", "get"))

    .abundance <- rlang::enquo(.abundance)

    if (distmethod %in% distMethods$vegdist ){
        x <- .data %>% mp_extract_abundance(.abundance=!!.abundance, byRow=FALSE)
        
        nmds <- withr::with_seed(seed, vegan::metaMDS(x, distance=distmethod, k=.dim, ...))

    }else{
        if(! distmethod %in% colnames(.data@colData)){
            if (rlang::quo_is_missing(.abundance)){
                rlang::abort("The .abundance must be required, when the distmethod is not present in the object.")
            }else{
                .data %<>% mp_cal_dist(.abundance=!!.abundance, distmethod=distmethod, action="add")
            }
        }

        distobj <- .data %>%
                   mp_extract_dist(distmethod=distmethod) 

        nmds <- withr::with_seed(seed, vegan::metaMDS(distobj, distance=distmethod, k=.dim, ...))
    }

    if (action=="get"){
        return(nmds)
    }

    dat <- nmds %>% 
           tidydr(display=c("sites", "features"), 
                  choices=seq_len(.dim))

    da <- .data %>%
          mp_extract_sample() %>%
          dplyr::left_join(
                 dat,
                 by=c("Sample"="sites")
          )

    if (action=="only"){
        da %<>%
             add_attr(dat %>% attr("features_tb"), name="features_tb") %>%
             add_internal_attr(object=nmds, name="NMDS")
        return(da)
    }else if (action=="add"){
        .data@colData <- da %>%
                         tibble::column_to_rownames(var="Sample") %>%
                         S4Vectors::DataFrame(check.names=FALSE)
        .data %<>% add_internal_attr(object=nmds, name="NMDS")
        return (.data)
    }
})

.internal_cal_nmds <- function(.data, .abundance, distmethod="bray", .dim=2, action="add", seed=123, ...){

    action %<>% match.arg(c("add", "only", "get"))

    .abundance <- rlang::enquo(.abundance)

    if (distmethod %in% distMethods$vegdist ){
        x <- .data %>% mp_extract_abundance(.abundance=!!.abundance, byRow=FALSE)

        nmds <- withr::with_seed(seed, vegan::metaMDS(x, distance=distmethod, k=.dim, ...))

    }else{
        if(! distmethod %in% colnames(.data)){
            if (rlang::quo_is_missing(.abundance)){
                rlang::abort("The .abundance must be required, when the distmethod is not present in the object.")
            }else{
                .data %<>% mp_cal_dist(.abundance=!!.abundance, distmethod=distmethod, action="add")
            }
        }

        distobj <- .data %>%
                   mp_extract_dist(distmethod=distmethod)

        nmds <- withr::with_seed(seed, vegan::metaMDS(distobj, distance=distmethod, k=.dim, ...))
    }

    dat <- nmds %>%
           tidydr(display=c("sites", "features"), choices=seq_len(.dim))

    if (action=="get"){
        return(nmds)
    }else if (action=="only"){
        da <- .data %>%
              mp_extract_sample() %>%
              dplyr::left_join(
                  dat,
                  by=c("Sample"="sites")
              ) %>%
              add_attr(dat %>% attr("features_tb"), name="features_tb") %>%
              add_internal_attr(object=nmds, name="NMDS")
        return(da)
    }else if (action=="add"){
        .data %<>%
            dplyr::left_join(
                dat,
                by=c("Sample"="sites")
            ) %>%
            add_attr(dat %>% attr("features_tb"), name="features_tb") %>%
            add_internal_attr(object=nmds, name="NMDS")
        return(.data)
    }
}

#' @rdname mp_cal_nmds-methods
#' @aliases mp_cal_nmds,tbl_mpse
#' @exportMethod mp_cal_nmds
setMethod("mp_cal_nmds", signature(.data="tbl_mpse"), .internal_cal_nmds)

#' @rdname mp_cal_nmds-methods
#' @aliases mp_cal_nmds,grouped_df_mpse
#' @exportMethod mp_cal_nmds
setMethod("mp_cal_nmds", signature(.data="grouped_df_mpse"), .internal_cal_nmds) 

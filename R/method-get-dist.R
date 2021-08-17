#' @title calculate distance
#' @param obj phyloseq, phyloseq class or data.frame
#' nrow sample * ncol feature. 
#' @param method character, default is hellinger, 
#' see alse \code{\link[vegan]{decostand}} 
#' @param distmethod character, default is "euclidean", 
#' see also \code{\link[phyloseq]{distanceMethodList}}
#' @param taxa_are_rows logical, default is FALSE.
#' @param sampleda data.frame, nrow sample * ncol factors.
#' @param tree object, the phylo class, see also \code{\link[ape]{as.phylo}}.
#' @param ..., additional parameters.
#' @return distance class contianed distmethod and originalD attr
#' @export
#' @examples
#' \dontrun{
#' data(test_otu_data)
#' distclass <- get_dist(test_otu_data)
#' hcsample <- get_clust(distclass)
#' }
get_dist <- function(obj,...){
    UseMethod("get_dist")
}

#' @method get_dist data.frame
#' @rdname get_dist
#' @importFrom vegan decostand
#' @export
get_dist.data.frame <- function(obj, 
    distmethod="euclidean",
    taxa_are_rows=FALSE,	
    sampleda=NULL,
    tree=NULL,
    method="hellinger",
    ...){
    objphyloseq <- phyloseq::phyloseq(
                       otu_table=phyloseq::otu_table(obj, 
                       taxa_are_rows=taxa_are_rows),
                       sam_data=phyloseq::sample_data(sampleda),
                       phy_tree=tree)
    return(get_dist.phyloseq(objphyloseq, 
                             distmethod=distmethod, 
                             method=method,
                             ...))
    
}

#' @method get_dist phyloseq
#' @seealso \code{\link[phyloseq]{distance}}
#' @rdname get_dist
#' @export
get_dist.phyloseq <- function(obj, distmethod="euclidean", method="hellinger",...){
    tmpmethod <- gsub("^(u.*)*unifrac$", "unifrac", distmethod, ignore.case = TRUE)
    tmpmethod <- gsub("^w.*unifrac$", "wunifrac", distmethod, ignore.case = TRUE) 
    tree <- obj@phy_tree
    if (tmpmethod=="unifrac" || tmpmethod=="wunifrac"){
    	if(is.null(tree)){
    		stop("The tree should be required when the distmethod is `WeightUniFrac` or `UnWeightUniFrac`")
    	}
    }
    if (!is.null(method)){
    	if (obj@otu_table@taxa_are_rows){
    		tmpotu <- t(obj@otu_table)
    	}else{
    		tmpotu <- data.frame(obj@otu_table)
    	}
    	obj@otu_table <- phyloseq::otu_table(decostand(tmpotu, method=method), 
                                   taxa_are_rows=FALSE)
    }
    disres <- phyloseq::distance(obj, method=distmethod, type="sample", ...)
    attr(disres, "distmethod") <- distmethod
    attr(disres, "originalD") <- data.frame(obj@otu_table, check.names=FALSE)
    return(disres)
}


#' Calculate the distances between the samples with specified abundance.
#'
#' @rdname mp_cal_dist-methods
#' @param .data MPSE or tbl_mpse object
#' @param .abundance the name of otu abundance to be calculated
#' @param .env the column names of continuous environment factors,
#' default is NULL.
#' @param distmethod character the method to calculate distance.
#' option is "manhattan", "euclidean", "canberra", "bray", "kulczynski", 
#' "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup",
#' "binomial", "chao", "cao" (implemented in vegdist of vegan), and
#' "w", "-1", "c", "wb", "r", "I", "e", "t", "me", "j", "sor", "m", "-2", "co"
#' "cc", "g", "-3", "l", "19", "hk", "rlb", "sim", "gl", "z" (implemented in 
#' betadiver of vegan), "maximum", "binary", "minkowski" (implemented in dist 
#' of stats), "unifrac", "weighted unifrac" (implemented in phyloseq),
#' @param action character, "add" joins the distance data to the object, "only" return
#' a non-redundant tibble with the distance information. "get" return 'dist' object.
#' @param scale logical whether scale the metric of environment (.env is provided) before
#' the distance was calculated, default is FALSE. The environment matrix can be processed
#' when it was joined to the MPSE or tbl_mpse object.
#' @param ... additional parameters.
#'
#' some dot arguments if \code{distmethod} is \code{unifrac} or \code{weighted unifrac}:
#' \itemize{
#'   \item \code{weighted} logical, whether to use weighted-UniFrac calculation, which considers the relative 
#'     abundance of taxa, default is FALSE, meaning unweightrd-UniFrac, which only considers 
#'     presence/absence of taxa.
#'   \item \code{normalized} logical, whether normaized the branch length of tree to the range between 0 and 1 
#'     when the \code{weighted}=\code{TRUE}.
#'   \item \code{parallel} logical, whether to execute the calculation in parallel, default is FALSE.
#' }
#' @return update object or tibble according the 'action'
#' @author Shuangbin Xu
#' @export
#' @examples
#' data(mouse.time.mpse)
#' mouse.time.mpse %<>%
#'     mp_decostand(.abundance=Abundance) %>% 
#'     mp_cal_dist(.abundance=hellinger, distmethod="bray")
#' # Visualization
#' library(ggplot2)
#' tbl <- mouse.time.mpse %>%
#'        mp_extract_dist(distmethod="bray", .group=time)
#' tbl
#' tbl %>% 
#'   ggplot(aes(x=GroupsComparison, y=bray)) + 
#'   geom_boxplot(aes(fill=GroupsComparison)) + 
#'   geom_jitter(width=0.1) + 
#'   xlab(NULL) + 
#'   theme(legend.position="none")
setGeneric("mp_cal_dist", function(.data, .abundance, .env=NULL, distmethod="bray", action="add", scale=FALSE, ...)standardGeneric("mp_cal_dist"))

#' @rdname mp_cal_dist-methods
#' @aliases mp_cal_dist,MPSE
#' @importFrom rlang :=
#' @exportMethod mp_cal_dist
setMethod("mp_cal_dist", signature(.data="MPSE"), function(.data, .abundance, .env=NULL, distmethod="bray", action="add", scale=FALSE, ...){
    
    action %<>% match.arg(c("add", "get", "only"))

    .abundance <- rlang::enquo(.abundance)
    .env <- rlang::enquo(.env)

    dotargs <- list(...)

    otutree <- .data@otutree
    if (is.character(distmethod)){
        distmethod <- gsub("^(u.*)*unifrac$", "unifrac", distmethod, ignore.case = TRUE)
        distmethod <- gsub("^w.*unifrac$", "wunifrac", distmethod, ignore.case = TRUE)        
    }

    if (is.function(distmethod)){
        distfun <- distmethod
        f <- match.call()
        distmethod <- gsub(".*::", "\\2", rlang::call_args(f)$distmethod)
    }else if (distmethod %in% distMethods$vegdist){
        distfun <- vegan::vegdist 
    }else if (distmethod %in% distMethods$betadiver){
        distfun <- vegan::betadiver
    }else if (distmethod %in% distMethods$dist){
        distfun <- stats::dist
    }else if (distmethod %in% distMethods$UniFrac){
        distfun <- cal_Unifrac_dist
        if (is.null(otutree)){
            rlang::abort("The otutree slot is required for the UniFrac calculation.")
        }
        dotargs <- .internal_append(dotargs, list(tree=otutree@phylo))
        if (distmethod == "unifrac"){
            dotargs <- .internal_append(dotargs, list(weighted = FALSE))
        }
        if (distmethod == "wunifrac"){
            dotargs <- .internal_append(dotargs, list(weighted = TRUE))
        }
    }else{
        distfun <- vegan::designdist
    }

    if (rlang::quo_is_missing(.abundance) && rlang::quo_is_null(.env)){
        rlang::abort(c("The one of .abundance and .env should be provided", 
                     "The .abundance should be specified one column name of abundance of feature or",
                     "The .env should be specified names of continuous sample environment factor."))
    }else if(!rlang::quo_is_null(.env)){
        da <- .data %>% 
              mp_extract_sample() %>%
              column_to_rownames(var="Sample") %>%
              dplyr::select(!!.env)
        if (ncol(da)==1 && da %>% pull(!!.env) %>% rlang::is_list()){
            da %<>%
                as_tibble(rownames="Sample") %>%
                tidyr::unnest() %>%
                suppressWarnings() %>% 
                tidyr::pivot_wider(id_cols="Sample", 
                                   names_from=vapply(., is.character, logical(1)) %>% 
                                              select_true_nm(rm="Sample"), 
                                   values_from=vapply(., is.numeric, logical(1)) %>% 
                                               select_true_nm()
                                  ) %>%
                column_to_rownames(var="Sample")
        }
        if (scale){
            da %<>% scale()
        }
        distsampley <- paste0("Env_", distmethod, "Sampley")
    }else{

        da <- .data %>% 
              mp_extract_assays(.abundance=!!.abundance, byRow=FALSE)
        distsampley <- paste0(distmethod, "Sampley")
    }

    if (distmethod %in% distMethods$UniFrac){
        if (!rlang::quo_is_null(.env)){
            rlang::abort("The distance of sample based on environment factor is not calculated via UniFrac method.")
        }
        params <- c(list(da), dotargs)
    }else{
        params <- c(list(da, method=distmethod), dotargs)
    }

    da <- do.call(distfun, params)

    if (action=="get"){
        return(da)
    }
    
    if (!rlang::quo_is_null(.env)){
        distmethod <- paste0("Env_", distmethod)
    }

    dat <- da %>% 
        as.matrix %>% 
        corrr::as_cordf(diagonal=0) %>% 
        corrr::stretch(na.rm=FALSE, remove.dups=TRUE) %>%
        dplyr::rename(!!distmethod:="r", !!distsampley:="y") %>% 
        tidyr::nest(!!distmethod:=c(!!as.symbol(distsampley), !!as.symbol(distmethod)))

    dat <- .data %>% 
           mp_extract_sample() %>%
           left_join(dat, by=c("Sample"="x")) 

    if (action=="only"){
        return(dat)   

    }else if (action=="add"){
        .data@colData <- dat %>%
            column_to_rownames(var="Sample") %>%
            S4Vectors::DataFrame(check.names=FALSE)
        return(.data)
    }
          
})


.internal_cal_dist <- function(.data, .abundance, .env=NULL, distmethod="bray", action="add", scale=FALSE, ...){
    action %<>% match.arg(c("add", "get", "only"))

    .abundance <- rlang::enquo(.abundance)
    .env <- rlang::enquo(.env)

    if (is.character(distmethod)){
        distmethod <- gsub("^(u.*)*unifrac$", "unifrac", distmethod, ignore.case = TRUE)
        distmethod <- gsub("^w.*unifrac$", "wunifrac", distmethod, ignore.case = TRUE)
    }

    dotargs <- list(...)
    otutree <- .data %>% attr("otutree")

    if (is.function(distmethod)){
        distfun <- distmethod
        f <- match.call()
        distmethod <- gsub(".*::", "\\2", rlang::call_args(f)$distmethod)
    }else if (distmethod %in% distMethods$vegdist){
        distfun <- vegan::vegdist
    }else if (distmethod %in% distMethods$betadiver){
        distfun <- vegan::betadiver
    }else if (distmethod %in% distMethods$dist){
        distfun <- stats::dist
    }else if (distmethod %in% distMethods$UniFrac){
        distfun <- cal_Unifrac_dist
        if (is.null(otutree)){
            rlang::abort("The otutree slot is required for the UniFrac calculation.")
        }
        dotargs <- .internal_append(dotargs, list(tree=otutree@phylo))
        if (distmethod == "unifrac"){
            dotargs <- .internal_append(dotargs, list(weighted = FALSE))
        }
        if (distmethod == "wunifrac"){
            dotargs <- .internal_append(dotargs, list(weighted = TRUE))        
        }
    }else{
        distfun <- vegan::designdist
    }

    if (rlang::quo_is_missing(.abundance) && rlang::quo_is_null(.env)){
        rlang::abort(c("The one of .abundance and .env should be provided",
                     "The .abundance should be specified one column name of abundance of feature or",
                     "The .env should be specified names (>2) of continuous sample environment factor."))
    }else if(!rlang::quo_is_null(.env)){
        da <- .data %>%
              mp_extract_sample() %>%
              column_to_rownames(var="Sample") %>%
              dplyr::select(!!.env)
        if (ncol(da)==1 && da %>% pull(!!.env) %>% rlang::is_list()){
            da %<>%
                as_tibble(rownames="Sample") %>%
                tidyr::unnest() %>%
                suppressWarnings() %>%
                tidyr::pivot_wider(id_cols="Sample",
                                   names_from=vapply(., is.character, logical(1)) %>%
                                              select_true_nm(rm="Sample"),
                                   values_from=vapply(., is.numeric, logical(1)) %>%
                                               select_true_nm()
                                  ) %>%
                column_to_rownames(var="Sample")
        }
        if (scale){
            da %<>% scale()
        }
        distsampley <- paste0("Env_", distmethod, "Sampley")
    }else{

        da <- .data %>%
              mp_extract_assays(.abundance=!!.abundance, byRow=FALSE)
        distsampley <- paste0(distmethod, "Sampley")
    }    

    if (distmethod %in% distMethods$UniFrac){
        if (!rlang::quo_is_null(.env)){
            rlang::abort("The distance of sample based on environment factor is not calculated via UniFrac method.")
        }
        params <- c(list(da), dotargs)
    }else{
        params <- c(list(da, method=distmethod), dotargs)
    }

    da <- do.call(distfun, params)
    
    if (action=="get"){
        return(da)
    }

    if (!rlang::quo_is_null(.env)){
        distmethod <- paste0("Env_", distmethod)
    }    

    dat <- da %>%
        as.matrix %>%
        corrr::as_cordf(diagonal=0) %>%
        corrr::stretch(na.rm=FALSE, remove.dups=TRUE) %>%
        dplyr::rename(!!distmethod:="r", !!distsampley:="y") %>%
        tidyr::nest(!!distmethod:=c(!!as.symbol(distsampley), !!as.symbol(distmethod)))

    if (action=="only"){
        dat <- .data %>% 
               mp_extract_sample() %>%
               dplyr::left_join(dat, by=c("Sample"="x")) 
        return(dat)
    }else if (action=="add"){
        .data %<>% 
            dplyr::left_join(dat, by=c("Sample"="x"))
        return(.data)
    }
}

#' @rdname mp_cal_dist-methods
#' @aliases mp_cal_dist,tbl_mpse
#' @exportMethod mp_cal_dist
setMethod("mp_cal_dist", signature(.data="tbl_mpse"), .internal_cal_dist)

#' @rdname mp_cal_dist-methods
#' @aliases mp_cal_dist,grouped_df_mpse
#' @exportMethod mp_cal_dist
setMethod("mp_cal_dist", signature(.data="grouped_df_mpse"), .internal_cal_dist)

################################################################################
# Fast UniFrac for R.
# Adapted from The ISME Journal (2010) 4, 17-27; doi:10.1038/ismej.2009.97
#
# This is from original implementation in phyloseq implemented by
# Paul J. McMurdie (https://github.com/joey711/phyloseq)
################################################################################

#' @importFrom foreach %dopar%
cal_Unifrac_dist <- function(x, tree, weighted = FALSE, normalized = TRUE, parallel=FALSE){
    if (is.null(tree$edge.length)){
        rlang::abort("The branch length of tree is none, the UniFrac distance is not be calculated.")
    }

    x <- x %>% t()
    tree <- .internal_reroot(tree)
    
    if (!all(tree$tip.label == rownames(x))){
        x <- x[tree$tip.label,]
    }

    ntip <- ape::Ntip(tree)
    node_num <- ape::Nnode(tree, internal.only=FALSE)
    # The original method might generate error when a node has three or larger than tree children nodes.
    #node.desc <- matrix(tree$edge[order(tree$edge[,1]), 2], byrow=TRUE, ncol=2) 
    edge <- tree$edge

    edge_array <- matrix(0, nrow=node_num, ncol=ncol(x),
                         dimnames=list(NULL, sample_names=colnames(x)))

    edge_array[1:ntip, ] <- x

    ord.node <- order(ape::node.depth(tree))[seq(from=ntip+1, to=node_num, by=1)]

    for(i in ord.node){
        # The original method might generate error when a node has three or larger than tree children nodes.
        #edge_array[i,] <- colSums(edge_array[node.desc[i-ntip,], , drop=FALSE], na.rm = TRUE)
        edge_array[i, ] <- colSums(edge_array[edge[edge[,1]==i, 2], ,drop=FALSE], na.rm=TRUE)
    }
    
    spn <- utils::combn(colnames(x), 2, simplify=FALSE)
    edge_array <- edge_array[tree$edge[,2], ]
    samplesums <- colSums(x)

    if(!weighted){
        # For unweighted UniFrac, convert the edge_array to an occurrence (presence/absence binary) array
        edge_occ <- (edge_array > 0) - 0
    }
    if( weighted & normalized ){
        # This is only relevant to weighted-UniFrac.
        # For denominator in the normalized distance, we need the age of each tip.
        # 'z' is the tree in postorder order used in calls to .C
        # Descending order of left-hand side of edge (the ancestor to the node)
        z = ape::reorder.phylo(tree, order="postorder")
        # Call phyloseq-internal function that in-turn calls ape's internal
        # horizontal position function, in C, using the re-ordered phylo object, `z`
        tipAges = ape::node.depth.edgelength(tree)
        # Keep only the tips, and add the tip labels in case `z` order differs from `tree`
        tipAges <- tipAges[seq_len(ntip)]
        names(tipAges) <- z$tip.label
        # Explicitly re-order tipAges to match x
        tipAges <- tipAges[rownames(x)]
    }

    if( !parallel ){ foreach::registerDoSEQ() }

    distlist <- foreach::foreach( i = spn, .packages="MicrobiotaProcess") %dopar% {
        A  <- i[1]
        B  <- i[2]
        AT <- samplesums[A]
        BT <- samplesums[B]
        if( weighted ){
            # weighted UniFrac
            wUF_branchweight <- abs(edge_array[, A]/AT - edge_array[, B]/BT)
            # calculate the w-UF numerator
            numerator <- sum({tree$edge.length * wUF_branchweight}, na.rm = TRUE)
            # if not-normalized weighted UniFrac, just return "numerator";
            # the u-value in the w-UniFrac description
            if(!normalized){
                return(numerator)
            } else {
                # denominator (assumes tree-indices and otu_table indices are same order)
                denominator <- sum({tipAges * (x[, A]/AT + x[, B]/BT)}, na.rm = TRUE)
                # return the normalized weighted UniFrac values
                return(numerator / denominator)
            }
        } else {
            # Unweighted UniFrac
            # Subset matrix to just columns A and B
            edge_occ_AB <- edge_occ[, c(A, B)]
            # Keep only the unique branches. Sum the lengths
            edge_uni_AB_sum <- sum((tree$edge.length * edge_occ_AB)[rowSums(edge_occ_AB, na.rm=TRUE) < 2, ], na.rm=TRUE)
            # Normalize this sum to the total branches among these two samples, A and B
            uwUFpairdist <- edge_uni_AB_sum / sum(tree$edge.length[rowSums(edge_occ_AB, na.rm=TRUE) > 0])
            return(uwUFpairdist)
        }
    }
    # Initialize UniFracMat with NAs
    UniFracMat <- matrix(NA_real_, ncol(x), ncol(x))
    rownames(UniFracMat) <- colnames(UniFracMat) <- colnames(x)
    # Matrix-assign lower-triangle of UniFracMat. Then coerce to dist and return.
    matIndices <- do.call(rbind, spn)[, 2:1]
    # Take care of edge case where there are two samples -> 1 pair of indices -> rbind doesn't return a matrix
    if(!is.matrix(matIndices)) matIndices <- matrix(matIndices, ncol=2)
    UniFracMat[matIndices] <- unlist(distlist)
    return(stats::as.dist(UniFracMat))
}

.internal_reroot <- function(tree){
    if (ape::is.rooted(tree)){
        return(tree)
    }
    ranlabel <- withr::with_seed(123, sample(tree$tip.label, 1))

    tree <- ape::root(tree, 
                      outgroup=ranlabel, 
                      edgelabel = TRUE, 
                      resolve.root = TRUE)

    if (!ape::is.rooted(tree)){
        rlang::abort("The rooted tree is required for UniFrac calculation. 
                     You can refer to treeio::root or ape::root.")
    }else{
        return(tree)
    }
}

.internal_append <- function(x, value){
    newarg <- setdiff(names(value), names(x))

    if (length(newarg)==0){
        return(x)
    }

    x <- append(x, value[newarg])
    return(x)
}

distMethods <- list(
    vegdist    = c("manhattan", "euclidean", "canberra", "bray",
                   "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn",
                   "mountford", "raup" , "binomial", "chao", "cao"),
    betadiver  = c("w", "-1", "c", "wb", "r", "I", "e", "t", "me", "j",
                   "sor", "m", "-2", "co", "cc", "g", "-3", "l", "19", "hk", "rlb",
                   "sim", "gl", "z"),
    dist       = c("maximum", "binary", "minkowski"),
    UniFrac    = c("unifrac", "wunifrac"),
    JSD        = "jsd",
    designdist = "ANY"
)

select_true_nm <- function(x, rm=NULL){
    dat <- x[x] %>% names()
    if (!is.null(rm)){
       dat <- dat[!dat %in% rm]
    }
    return(dat)
}

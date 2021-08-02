#' Mantel and Partial Mantel Tests for MPSE or tbl_mpse Object
#' @rdname mp_mantel-methods
#' @param .data MPSE or tbl_mpse object
#' @param .abundance the name of otu abundance to be calculated
#' @param .y.env the column names of continuous environment factors to
#' perform Mantel statistic, it is required.
#' @param .z.env the column names of continuous environment factors to 
#' perform Partial Mantel statistic based on this, default is NULL.
#' @param distmethod character the method to calculate distance based on .abundance.
#' @param distmethod.y character the method to calculate distance based on .y.env. 
#' @param distmethod.z character the method of calculated distance based on .z.env
#' @param method character Correlation method, options is "pearson", "spearman" or
#' "kendall"
#' @param permutations the number of permutations required, default is 999.
#' @param action character, "add" joins the mantel result to the internal attributes
#' of the object, "only" and "get" return 'mantel' or 'mantel.partial' (if .z.env is
#' provided) object.
#' @param seed a random seed to make the analysis reproducible, default is 123.
#' @param scale.y logical whether scale the environment matrix (.y.env) before 
#' the distance is calculated, default is FALSE
#' @param scale.z logical whether scale the environment matrix (.z.env) before
#' the distance is calculated, default is FALSE
#' @param ... additional parameters, see also \code{\link[vegan]{mantel}}.
#' @return update object or tibble according the 'action'
#' @seealso \code{\link[vegan]{mantel}}
#' @export
#' @examples
#' library(vegan)
#' data(varespec, varechem)
#' mpse <- MPSE(assays=list(Abundance=t(varespec)), colData=varechem)
#' mpse %>% mp_mantel(.abundance=Abundance, 
#'                    .y.env=colnames(varechem),
#'                    distmethod.y="euclidean",
#'                    scale.y = TRUE
#'                    )
setGeneric("mp_mantel", 
           function(.data, 
                    .abundance, 
                    .y.env, 
                    .z.env=NULL, 
                    distmethod="bray", 
                    distmethod.y="euclidean", 
                    distmethod.z="euclidean", 
                    method="pearson", 
                    permutations=999,
                    action="get",
                    seed=123,
                    scale.y = FALSE,
                    scale.z = FALSE,
                    ...){
               standardGeneric("mp_mantel")
})

.internal_mp_mantel <- function(
                   .data,
                   .abundance,
                   .y.env,
                   .z.env=NULL,
                   distmethod="bray",
                   distmethod.y="euclidean",
                   distmethod.z="euclidean",
                   method="pearson", 
                   permutations=999, 
                   action="get",
                   seed = 123,
                   scale.y = FALSE,
                   scale.z = FALSE,
                   ...
                   ){

    action %<>% match.arg(c("add", "only", "get"))
    .abundance <- rlang::enquo(.abundance)
    .y.env <- rlang::enquo(.y.env)
    .z.env <- rlang::enquo(.z.env)

    if (rlang::quo_is_missing(.y.env)){
        rlang::abort(c("The .y.env is required to calculated the distance via distmethod.y method",
                       "The name of continuous environment factor should be provided."))
    }

    if(! distmethod %in% colnames(.data %>% mp_extract_sample())){
        if (rlang::quo_is_missing(.abundance)){
            rlang::abort("The .abundance must be required, when the distmethod is not present in the object.")
        }else{
            .data %<>% mp_cal_dist(.abundance=!!.abundance, distmethod=distmethod, action="add")
        }
    }

    xdist <- .data %>% mp_extract_dist(distmethod=distmethod)

    ydist <- .data %>% mp_cal_dist(.env=!!.y.env, distmethod=distmethod.y, action="get", scale=scale.y)
    
    if (!rlang::quo_is_null(.z.env)){
        zdist <- .data %>% mp_cal_dist(.env=!!.z.env, distmethod=distmethod.z, action="get", scale=scale.z)
        res <- withr::with_seed(seed, vegan::mantel.partial(xdis=xdist, ydis=ydist, zdis=zdist, method=method, permutations=permutations, ...))
    }else{
        res <- withr::with_seed(seed, vegan::mantel(xdis=xdist, ydis=ydist, method=method, permutations=permutations, ...))
    }

    if (action %in% c("get")){
        return(res)
    }else if(action=="add"){
        print(res)
        resnm <- ifelse(!rlang::quo_is_null(.z.env), "mantel.partial", "mantel")
        message(paste0("The result of ",resnm," has been saved to the internal attribute !"))
        message(paste0("It can be extracted using this-object %>% mp_extract_internal_attr(name='",resnm,"')"))
        .data %<>%
            add_internal_attr(object=res, name=toupper(resnm))
        return(.data)    
    }else{
        print(res)
        da <- mp_fortify(res)
        return(da)
    }
}

#' @rdname mp_mantel-methods
#' @aliases mp_mantel,MPSE
#' @exportMethod mp_mantel
setMethod("mp_mantel", signature(.data="MPSE"), .internal_mp_mantel)

#' @rdname mp_mantel-methods
#' @aliases mp_mantel,tbl_mpse
#' @exportMethod mp_mantel
setMethod("mp_mantel", signature(.data="tbl_mpse"), .internal_mp_mantel)

#' @rdname mp_mantel-methods
#' @aliases mp_mantel,grouped_df_mpse
#' @exportMethod mp_mantel
setMethod("mp_mantel", signature(.data="grouped_df_mpse"), .internal_mp_mantel)

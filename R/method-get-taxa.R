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
#' @importFrom phyloseq otu_table tax_table taxa_are_rows rank_names
#' @export
setMethod("get_taxadf", "phyloseq", function(obj, taxlevel=2, type="species",...){
    if (is.null(obj@tax_table)){
    	stop("The tax table is empty!")
    }else{
    	taxdf <- tax_table(obj)
    }
    otuda <- checkotu(obj)
    sampleda <- get_sample(obj)
    if (inherits(taxlevel, 'numeric')){taxlevel <- rank_names(obj)[taxlevel]}
    if (inherits(taxlevel, 'character')){
    	if (!taxlevel %in% rank_names(obj)){
    		stop("the taxlevel should be among the values of rank_names(phyloseq)")
    	}else{
    		taxlevel <- rank_names(obj)[match(taxlevel,rank_names(obj))]
    	}
    }
    #taxlevel <- rank_names(obj)[taxlevel]
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
#' @importFrom phyloseq phyloseq otu_table tax_table
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
        sampleda <- sample_data(sampleda)
    }
    taxda <- fillNAtax(taxda, type=type)
    if (inherits(taxlevel, "numeric")){taxlevel <- colnames(taxda)[taxlevel]}
    tmptax <- taxda[, match(taxlevel, colnames(taxda)), drop=FALSE]
    tmptaxda <- taxda[, seq(from=1, to=match(taxlevel, colnames(taxda))), drop=FALSE]
    tmptaxda <- tmptaxda[!duplicated(tmptaxda),,drop=FALSE]
    rownames(tmptaxda) <- as.vector(tmptaxda[, match(taxlevel, colnames(tmptaxda))])
    taxdf <- otu_table(get_count(data=obj, 
                                 featurelist=tmptax), 
                       taxa_are_rows=TRUE)
    taxdf <- phyloseq(taxdf, sampleda, tax_table(as.matrix(tmptaxda)))
    return(taxdf)
})

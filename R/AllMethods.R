#' @title taxonomy barplot
#' @param obj phyloseq, phyloseq class or data.frame, 
#' (nrow sample * ncol feature (factor)) or the data.frame for geom_bar.
#' @param mapping set of aesthetic mapping of ggplot2, default is NULL,
#' if the data is the data.frame for geom_bar, the mapping should be set.
#' @param position character, default is `stack`. 
#' @param stat character, default is `identity`.
#' @param width numeric, the width of bar, default is 0.7.
#' @param topn integer, the top number of abundance taxonomy(feature).
#' @param count logical, whether show the relative abundance.  
#' @param sampleda data.frame, (nrow sample * ncol factor), the sample 
#' information, if the data doesn't contain the information.
#' @param factorLevels list, the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this.
#' @param facetNames character, default is NULL.
#' @param plotgroup logical, whether calculate the mean or median etc 
#' for each group, default is FALSE.
#' @param groupfun character, how to calculate for feature in each group,
#' the default is `mean`, this will plot the mean of feature in each group.
#' @param ... additional parameters, see \code{\link[ggplot2]{ggplot}}
#' @return barplot of tax
#' @author Shuangbin Xu
#' @export
#' @examples
#' library(ggplot2)
#' data(test_otu_data)
#' otubar <- ggbartax(test_otu_data) + 
#'          xlab(NULL) + ylab("relative abundance(%)")
ggbartax <- function(obj,...){
    UseMethod("ggbartax")
}

#' @method ggbartax phyloseq
#' @importFrom phyloseq otu_table taxa_are_rows
#' @rdname ggbartax
#' @export
ggbartax.phyloseq <- function(obj, ...){
    if (is.null(obj@otu_table)){
    	stop("The otu table is empty!")
    }else{
    	otudata <- get_otudata(obj)
    }
    if (!is.null(obj@sam_data)){
    	sampleda <- data.frame(sample_data(obj), check.names=FALSE)
    	p <- ggbartax.default(obj=otudata, sampleda=sampleda, ...)
    }else{
    	p <- ggbartax.default(obj=otudata,...)
    }
    return(p)	
}

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

#' @title Rarefaction alpha index
#' @param obj phyloseq, phyloseq class or data.frame
#' shape of data.frame (nrow sample * ncol feature ( + factor)). 
#' @param shadow logical, whether merge samples with group (factorNames) and
#' display the ribbon of group, default is TRUE.
#' @param linesize integer, default is 0.5. 
#' @param chunks integer, the number of subsample in a sample,
#'  default is 400.
#' @param sampleda data.frame, (nrow sample * ncol factor)
#' @param factorNames character, default is missing.
#' @param facetnrow integer, the nrow of facet, default is 1.
#' @param factorLevels list, the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this.
#' @param indexNames character, default is "Observe",
#' only for "Observe", "Chao1", "ACE", "Shannon", "Simpson", "J".
#' @param se logical, default is FALSE.
#' @param method character, default is lm. 
#' @param formula formula, default is `y ~ log(x)`
#' @param ... additional parameters, 
#' see also \code{\link{ggplot2}{ggplot}}.
#' @return figure of rarefaction curves
#' @author Shuangbin Xu
#' @export
#' @examples
#' data(test_otu_data)
#' library(ggplot2)
#' prare <- ggrarecurve(test_otu_data,
#'                indexNames=c("Observe","Chao1","ACE"),
#'                shadow=FALSE,
#'                factorNames="group"
#'          ) +
#'          theme(legend.spacing.y=unit(0.02,"cm"),
#'                legend.text=element_text(size=6))
ggrarecurve <- function(obj, ...){
    UseMethod("ggrarecurve")
}

#' @method ggrarecurve phyloseq
#' @rdname ggrarecurve
#' @export
ggrarecurve.phyloseq <- function(obj, chunks=400, factorLevels=NULL, ...){
    #otuda <- checkotu(obj)
    #sampleda <- data.frame(get_sample(obj),check.names=FALSE)
    obj <- get_rarecurve(obj, chunks=chunks, factorLevels=factorLevels)
    p <- ggrarecurve.rarecurve(obj=obj, ...)
    return(p)	
}


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
#' data(test_otu_data)
#' vennlist <- get_vennlist(test_otu_data, 
#'                  factorNames="group")
#' vennlist
#' #library(VennDiagram)
#' #venn.diagram(vennlist, height=5, 
#' #             width=5, filename = "./test_venn.pdf", 
#' #             alpha = 0.85, fontfamily = "serif", 
#' #             fontface = "bold",cex = 1.2, 
#' #             cat.cex = 1.2, cat.default.pos = "outer",
#' #             cat.dist = c(0.22,0.22,0.12,0.12), 
#' #             margin = 0.1, lwd = 3, 
#' #             lty ='dotted', 
#' #             imagetype = "pdf")
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

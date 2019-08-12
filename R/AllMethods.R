#' @title taxonomy barplot
#' @param obj data.frame or phyloseq class. data.frame should be 
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
#' @param settheme logical, default is TRUE, or you can set FALSE, then
#' set the theme by youself.
#' @param facetNames character, default is NULL.
#' @param setColors logical, default is TRUE, or you can set FALSE, then 
#' set colors by `scale_fill_manual` of `ggplot2`.
#' @param ... additional parameters, see \code{\link[ggplot2]{ggplot}}
#' @return barplot of tax
#' @author ShuangbinXu
#' @export
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
		otudata <- getotudata(obj)
	}
	if (!is.null(obj@sam_data)){
		sampleda <- data.frame(sample_data(obj), check.names=FALSE)
		p <- ggbartax.default(otudata, sampleda=sampleda, ...)
	}else{
		p <- ggbartax.default(otudata,...)
	}
	return(p)	
}

#' @title get the data of specified taxonomy
#' @param obj phyloseq class or data.frame (default), the shape of data.frame should be
#' row sample * column feature.
#' @param taxda data.frame, the classifies of feature contained in obj.
#' @param taxlevel character, the column names of taxda that you want to get.
#' when the input is phyloseq class, the common names of tax rank is "Kingdom", 
#' "Phylum", "Class", "Order", "Family", "Genus", "Species", but you can use 1 to 7.
#' @param sampleda data.frame, the sample information.
#' @param ... additional parameters, see also \code{\link[MicrobiotaProcess]{CountOrRatios}}
#' @return phyloseq class contained tax data.frame and sample information.
#' @author ShuangbinXu
#' @export
gettaxdf <- function(obj,...){
	UseMethod("gettaxdf")
}

#' @method gettaxdf phyloseq
#' @importFrom phyloseq otu_table tax_table taxa_are_rows rank_names
#' @rdname gettaxdf
#' @export
gettaxdf.phyloseq <- function(obj, taxlevel=2, ...){
	if (is.null(obj@tax_table)){
		stop("The tax table is empty!")
	}else{
		taxdf <- tax_table(obj)
	}
	otuda <- checkotu(obj)
	sampleda <- getsample(obj)
	if (inherits(taxlevel, 'numeric')){taxlevel <- rank_names(obj)[taxlevel]}
	if (inherits(taxlevel, 'character')){
		if (!taxlevel %in% rank_names(obj)){
			stop("the taxlevel should be among the values of rank_names(phyloseq)")
		}else{
			taxlevel <- rank_names(obj)[match(taxlevel,rank_names(obj))]
		}
	}
	#taxlevel <- rank_names(obj)[taxlevel]
	taxdf <- gettaxdf.default(otuda, 
							  taxda=taxdf, 
							  taxlevel=taxlevel,
							  sampleda=sampleda,...)
	return(taxdf)
}

#' @method gettaxdf default
#' @importFrom phyloseq otu_table tax_table
#' @rdname gettaxdf
#' @export 
gettaxdf.default <- function(otuda, taxda, 
							 taxlevel,
							 sampleda=NULL,
							 ...){
	if (!isTRUE(taxa_are_rows)){
		otuda <- data.frame(t(otuda), check.names=FALSE)
	}
	if(!is.null(sampleda) && !inherits(sampleda, "sample_data")){
		sampleda <- sample_data(sampleda)
	}
	taxda <- fillNAtax(taxda)
	tmptax <- taxda[,match(taxlevel, colnames(taxda)), drop=FALSE]
	taxdf <- otu_table(CountOrRatios(otuda, 
									 tmptax, 
									 rownamekeep=FALSE,...), 
					   taxa_are_rows=TRUE)
	taxdf <- new("phyloseq",
				 otu_table=taxdf,
				 sam_data=sampleda)
	return(taxdf)

}

#' @title Rarefaction alpha index
#' @param obj data.frame or phyloseq class, shape of data.frame (nrow sample * ncol feature (factor)) or
#' the data.frame for stat_smooth.
#' @param mapping, set of aesthetic mapping of ggplot2, default is NULL,
#' if the data is the data.frame for stat_smooth, the mapping should be set. 
#' @param linesize integer, default is 0.5. 
#' @param chunks integer, the number of subsample in a sample,
#'  default is 400.
#' @param sampleda, data.frame, (nrow sample * ncol factor)
#' @param factorNames character, default is missing.
#' @param facetnrow, the nrow of facet, default is 1.
#' @param factorLevels list, the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this.
#' @param indexNames vector character, default is "Observe", only for "Observe",
#' "Chao1", "ACE", "Shannon", "Simpson", "J".
#' @param se logical, default is FALSE.
#' @param method character, default is lm. 
#' @param formula formula, default is `y ~ log(x)`
#' @param ... additional parameters, see \code{\link{ggplot2}{ggplot}}.
#' @author ShuangbinXu
#' @export
ggrarecurve <- function(obj, ...){
	UseMethod("ggrarecurve")
}

#' @method ggrarecurve phyloseq
#' @rdname ggrarecurve
#' @export
ggrarecurve.phyloseq <- function(obj, ...){
	otuda <- checkotu(obj)
	sampleda <- data.frame(getsample(obj),check.names=FALSE)
	p <- ggrarecurve.default(data=otuda, sampleda=sampleda, ...)
	return(p)	
}


#' @title generate a vennlist for VennDiagram 
#' @param obj data.frame or phyloseq class, a dataframe contained one character column and the others are numeric.
#' all columns should be numeric if sampleinfo isn't NULL.
#' @param sampleinfo dataframe; a sample information, default is NULL.
#' @param  factorNames character, a column name of sampleinfo, 
#' when sampleinfo isn't NULL, factorNames shouldn't be NULL, default is NULL,
#' when the input is phyloseq, the factorNames should be provided. 
#' @param ... additional parameters, see \code{\link[MicrobitaProcess]{CountOrRatios}}.
#' @return return a list for VennDiagram.
#' @author ShuangbinXu
#' @export 
getvennlist <- function(obj,...){
	UseMethod("getvennlist")
}

#' @method getvennlist phyloseq
#' @rdname getvennlist
#' @export 
getvennlist.phyloseq <- function(obj, ...){
	otuda <- checkotu(obj)
	sampleda <- checksample(obj)
	#tmpfactors <- colnames(sampleda)[factorNamesIndex]
	vennlist <- getvennlist.default(da=otuda,
									sampleinfo=sampleda,
									#factorNames=factorNames,
									...)
	return(vennlist)
}

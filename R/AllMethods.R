#' @title ggbartax
#'
#' @param obj data.frame or phyloseq.
#' @param ... additional parameters
#' @return barplot of tax
#' @export
ggbartax <- function(obj,...){
	UseMethod("ggbartax")
}

#' @title ggbartax
#' @param obj the phyloseq object.
#' @param sampleda data.frame, nrow sample * ncol factor. 
#' @param ... additional parameters.
#' @importFrom phyloseq otu_table taxa_are_rows
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

#' @title gettaxdf
#' @param obj the phyloseq object
#' @param ... additional parameters.
#' @export
gettaxdf <- function(obj,...){
	UseMethod("gettaxdf")
}

#' @title gettaxdf
#' @param obj the phyloseq object.
#' @param taxlevel character, default is "Phylum", options is "Kingdom",
#' "Phylum", "Class", "Order", "Family", "Genus", "Species".
#' @param ...  additional parameters.
#' @importFrom phyloseq otu_table tax_table taxa_are_rows
#' @export
gettaxdf.phyloseq <- function(obj, taxlevel="Phylum", ...){
	if (is.null(obj@tax_table)){
		stop("The tax table is empty!")
	}else{
		taxdf <- tax_table(obj)
		taxdf[is.na(taxdf)] <- "Unknown"
	}
	otuda <- checkotu(obj)
	sampleda <- getsample(obj)
	taxdf <- gettaxdf.default(otuda, 
							  taxda=taxdf, 
							  taxlevel=taxlevel,
							  sampleda=sampleda,...)
	#tmptax <- taxdf[,match(taxlevel,colnames(taxdf)), drop=FALSE]
	#taxdf <- otu_table(CountOrRatios(otuda, tmptax, rownamekeep=FALSE, ...), 
	#				   taxa_are_rows=TRUE)
	#taxdf <- new("phyloseq", 
	#			 otu_table=taxdf,
	#			 sam_data=sampleda)
	return(taxdf)
}

#' @param otuda data.frame, otu table.
#' @param taxda data.frame, the dataframe of taxonomy.
#' @param taxlevel character, the column names of taxda that you want to get.
#' @param sampleda data.frame, the sample information.
#' @importFrom phyloseq otu_table tax_table
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
	tmptax <- taxdf[,match(taxlevel, colnames(taxdf)), drop=FALSE]
	taxdf <- otu_table(CountOrRatios(otuda, 
									 tmptax, 
									 rownamekeep=FALSE,...), 
					   taxa_are_rows=TRUE)
	taxdf <- new("phyloseq",
				 otu_table=taxdf,
				 sam_data=sampleda)
	return(taxdf)

}


#' @title ggrarecurve
#' @param obj the phyloseq object.
#' @param ... additional parameters.
#' @export
ggrarecurve <- function(obj, ...){
	UseMethod("ggrarecurve")
}

#' @param obj the phyloseq object.
#' @param ... additional parameters.
#' @export
ggrarecurve.phyloseq <- function(obj, ...){
	otuda <- checkotu(obj)
	sampleda <- data.frame(getsample(obj),check.names=F)
	p <- ggrarecurve.default(data=otuda, sampleda=sampleda, ...)
	return(p)	
}

#' @title getvennlist
#' @param obj the phyloseq object.
#' @param ... additional parameters.
#' @export 
getvennlist <- function(obj,...){
	UseMethod("getvennlist")
}

#' @param obj the phyloseq object.
#' @param factorNamesIndex integer, the index of column names of sample_data.
#' @param ... additional parameters.
#' @export 
getvennlist.phyloseq <- function(obj, factorNames, ...){
	otuda <- checkotu(obj)
	sampleda <- checksample(obj)
	#tmpfactors <- colnames(sampleda)[factorNamesIndex]
	vennlist <- getvennlist.default(da=otuda,
									sampleinfo=sampleda,
									factorNames=factorNames,
									...)
	return(vennlist)
}



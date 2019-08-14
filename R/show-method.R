#' @title show method
#' @name show
#' @docType methods
#' @rdname show-methods
#' @aliases show, diffAnalysisClass
#' @param object object, `diffAnalysisClass` class
#' @importFrom methods show
#' @exportMethod show
#' @author Shuangbin Xu
#' @return print info
#' @examples
#' data(kostic2012crc)
#' kostic2012crc
#' head(phyloseq::sample_data(kostic2012crc),3)
#' kostic2012crc <- phyloseq::rarefy_even_depth(kostic2012crc,rngseed=1024)
#' table(sample_data(kostic2012crc)$DIAGNOSIS)
#' diffres <- diffAnalysis(kostic2012crc, class="DIAGNOSIS",
#'                         mlfun="lda", filtermod="fdr",
#'                         firstcomfun = "kruskal.test",
#'                         firstalpha=0.05, strictmod=TRUE, 
#'                         secondcomfun = "wilcox.test",
#'                         submin=3, subclwilc=TRUE,
#'                         secondalpha=0.01, lda=3)
#' show(diffres)
setMethod("show", 
		  "diffAnalysisClass",
		  function(object){
			originalD <- object@originalD
		  	cat(paste0("The original data: ", ncol(originalD),
						 " features and ", nrow(originalD)," samples"),
				  fill=TRUE)
			sampleda <- object@sampleda
			cat(paste0("The sample data: ", ncol(sampleda), " variables and ", nrow(sampleda), " samples"),
				fill=TRUE)
			taxda <- object@taxda
			if(!is.null(taxda)){cat(paste0("The taxda contained ", nrow(taxda), " by ",ncol(taxda), " rank"),
									fill=TRUE)}
			else{cat("The taxda is NULL",fill=TRUE)}
			kwres <- object@kwres
			numfirstf <- nrow(kwres[kwres$pvalue<=0.05 & !is.na(kwres$pvalue),])
			cat(paste0("after first test (default is kruskal.test) number of feature (pvalue <=0.05):", 
					   numfirstf),
				fill=TRUE)
			secondvars <- getsecondTRUEvar(object)
			cat(paste0("after second test (default is wilcox.test) number of significantly discriminative feature:", 
					   nrow(secondvars)),
					   fill=TRUE)
			mlres <- tidydiffAnalysis(object) 
			uncertain <- length(grep("__un_", mlres$f))
			cat(paste0("after LDA or rf, Number of discriminative features: ", 
					   nrow(mlres), " (certain taxonomy classification:", 
					   nrow(mlres) -uncertain , 
					   "; uncertain taxonomy classication: ",uncertain,")"), 
				fill=TRUE)
		  }
)




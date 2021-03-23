#' @title method extensions to show for diffAnalysisClass objects.
#' @rdname show-methods
#' @param object object, `diffAnalysisClass` class
#' @importFrom methods show
#' @author Shuangbin Xu
#' @return print info
#' @export
#' @examples
#' # don't run in examples
#' #data(kostic2012crc)
#' #kostic2012crc
#' #head(phyloseq::sample_data(kostic2012crc),3)
#' #kostic2012crc <- phyloseq::rarefy_even_depth(kostic2012crc,rngseed=1024)
#' #table(phyloseq::sample_data(kostic2012crc)$DIAGNOSIS)
#' #set.seed(1024)
#' #diffres <- diff_analysis(kostic2012crc, classgroup="DIAGNOSIS",
#' #                        mlfun="lda", filtermod="fdr",
#' #                        firstcomfun = "kruskal.test",
#' #                        firstalpha=0.05, strictmod=TRUE, 
#' #                        secondcomfun = "wilcox.test",
#' #                        subclmin=3, subclwilc=TRUE,
#' #                        secondalpha=0.01, lda=3)
#' #show(diffres)
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
      firstfun <- extract_args(object, "firstcomfun")
      filtermod <- extract_args(object, "filtermod")
      alphafold <- extract_args(object, "firstalpha")
      cat(paste0("after first test (",firstfun,") number of feature (", filtermod,"<=",alphafold,"):", 
      numfirstf),fill=TRUE)
      secondvars <- get_second_true_var(object)
      secondfun <- extract_args(object, "secondcomfun")
      cat(paste0("after second test (",secondfun,") number of significantly discriminative feature:", 
      		   nrow(secondvars)),
      		   fill=TRUE)
      mlres <- as.data.frame(object) 
      uncertain <- length(grep("__un_", mlres$f))
      mlmethod <- extract_args(object, "mlfun")
      cat(paste0("after ",mlmethod,", Number of discriminative features: ", 
      		   nrow(mlres), " (certain taxonomy classification:", 
      		   nrow(mlres) -uncertain , 
      		   "; uncertain taxonomy classication: ",uncertain,")"), 
      	fill=TRUE)
    }
)


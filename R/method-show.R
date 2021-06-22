#' @title method extensions to show for diffAnalysisClass or alphasample objects.
#' @rdname show-methods
#' @param object object, diffAnalysisClass or alphasample class
#' @importFrom methods show
#' @author Shuangbin Xu
#' @return print info
#' @export
#' @examples
#' \dontrun{
#' data(kostic2012crc)
#' kostic2012crc
#' head(phyloseq::sample_data(kostic2012crc),3)
#' kostic2012crc <- phyloseq::rarefy_even_depth(kostic2012crc,rngseed=1024)
#' table(phyloseq::sample_data(kostic2012crc)$DIAGNOSIS)
#' set.seed(1024)
#' diffres <- diff_analysis(kostic2012crc, classgroup="DIAGNOSIS",
#'                         mlfun="lda", filtermod="fdr",
#'                         firstcomfun = "kruskal.test",
#'                         firstalpha=0.05, strictmod=TRUE, 
#'                         secondcomfun = "wilcox.test",
#'                         subclmin=3, subclwilc=TRUE,
#'                         secondalpha=0.01, lda=3)
#' show(diffres)
#' }
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

      
      firstvars <- extract_first_vars(obj=object) 
      firstfun <- extract_args(object, "firstcomfun")
      filtermod <- extract_args(object, "filtermod")
      alphafold <- extract_args(object, "firstalpha")
      cat(paste0("after first test (",firstfun,") number of feature (", filtermod,"<=",alphafold,"):", 
      length(firstvars)),fill=TRUE)
      secondvars <- get_second_true_var(object)
      #secondfun <- extract_args(object, "secondcomfun")
      secondfun <- check_second_fun(obj=object)
      cat(paste0("after second test (",secondfun,") number of significantly discriminative feature:", 
      		   nrow(secondvars)),
          fill=TRUE)
      mlres <- object@result 
      uncertain <- length(grep("__un_", mlres$f))
      mlmethod <- extract_args(object, "mlfun")
      cat(paste0("after ",mlmethod,", Number of discriminative features: ", 
      		   nrow(mlres), " (certain taxonomy classification:", 
      		   nrow(mlres) -uncertain , 
      		   "; uncertain taxonomy classication: ",uncertain,")"), 
      	fill=TRUE)
    }
)


extract_first_vars <- function(obj){
    filtermod <- extract_args(obj, "filtermod")
    firstalpha <- extract_args(obj, "firstalpha")
    kwres <- obj@kwres
    if (filtermod !="pvalue"){
        varsfirst <- kwres[kwres$fdr<=firstalpha& !is.na(kwres$fdr),,drop=FALSE]
    }else{
        varsfirst <- kwres[kwres$pvalue<=firstalpha&!is.na(kwres$pvalue),,drop=FALSE]
    }
    return (as.vector(varsfirst$f))
}

check_second_fun <- function(obj){
    subclass <- extract_args(obj, "subclass")
    classgroup <- extract_args(obj, "classgroup")
    strictmod <- extract_args(obj, "strictmod")
    clmin <- extract_args(obj, "clmin")
    clwilc <- extract_args(obj, "clwilc")
    fcfun <- extract_args(obj, "fcfun")
    secondcomfun <- extract_args(obj, "secondcomfun")
    subclmin <- extract_args(obj, "subclmin")
    subclwilc <- extract_args(obj, "subclwilc")
    if (!is.null(subclass) && strictmod){
        submin <- min(table(obj@sampleda[[subclass]]))
        if (submin >= subclmin && subclwilc){
            secondfun <- paste(secondcomfun, "and", fcfun)
        }else{
            secondfun <- fcfun
        }
    }else{
        groupmin <- min(table(obj@sampleda[[classgroup]]))
        if (groupmin >= clmin && clwilc){
            secondfun <- paste(secondcomfun, "and", fcfun)
        }else{
            secondfun <- fcfun
        }
    }
    return (secondfun)
}

#' @rdname show-methods
#' @exportMethod show
setMethod("show", "alphasample", function(object){
    print.alphasample(object)
})

print.alphasample <- function(object){
    msg <- "'alphasample' S4 object"
    sampleda <- object@sampleda
    if (!is.null(sampleda)){
        smsg <- sampleda
    }else{
        smsg <- "The 'alphasample' does not have sampleda slot."
    }
    cat (msg, fill=TRUE)
    cat ("The alpha diversity index :", fill=TRUE)
    print (head(object@alpha))
    if (is.character(smsg)){
        cat (smsg, fill=TRUE)
    }else{
        cat ("The sample information is present in 'alphasample':", fill=TRUE)
        print (head(sampleda))
    }
}

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

#' @aliases get_vennlist,data.framet
#' @rdname get_vennlist
#' @export
setMethod("get_vennlist", "data.frame", function(obj,
    sampleinfo=NULL,
    factorNames=NULL,...){
    if (!is.null(sampleinfo) && !is.null(factorNames)){
    	sampleinfo <- sampleinfo[,match(factorNames, colnames(sampleinfo)), 
    							 drop=FALSE]
    }
    if (!is.null(sampleinfo) && is.null(factorNames)){
    	stop("when sampleinfo isn't NULL, factorNames shouldn't be NULL")
    }
    obj <- get_count(data=obj, 
                     featurelist=sampleinfo)
    vennlist <- apply(obj, 1, function(x){names(x[x>0])})
    return(vennlist)
})

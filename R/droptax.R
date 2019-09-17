#' @title Dropping Species with Few abundance and Few Occurrences
#' 
#' @description
#' Drop species or features from the feature data frame or phyloseq that occur
#' fewer than or equal to a threshold number of occurrences and fewer 
#' abundance than to a threshold abundance.
#'
#' @param obj object, phyloseq or a dataframe of species (n_sample, n_feature).
#' @param minocc numeric, the threshold number of occurrences to be 
#' dropped, if < 1.0,it will be the threshold ratios of occurrences, 
#' default is 0.
#' @param minabu numeric, the threshold abundance, if fewer than the 
#' threshold will be dropped, default is 0.
#' @param ..., additional parameters.
#' @return dataframe of new features.
#' @export
#' @author Shuangbin Xu
#' @examples
#' otudafile <- system.file("extdata", "otu_tax_table.txt",
#'                          package="MicrobiotaProcess")
#' otuda <- read.table(otudafile, sep="\t", 
#'                     header=TRUE, row.names=1, 
#'                     check.names=FALSE, skip=1, 
#'                     comment.char="")
#' otuda <- otuda[sapply(otuda, is.numeric)]
#' dim(otuda)
#' otudat <- droptax(otuda, minocc=0.1, minabu=1)
#' dim(otudat)
#' data(test_otu_data)
#' keepps <- droptax(test_otu_data, minocc=0.1, minabu=0)
setGeneric("droptax",function(obj, ...){standardGeneric("droptax")})

#' @aliases droptax,data.frame
#' @rdname droptax
#' @export
setMethod("droptax", "data.frame",
    function(obj, minocc=0, minabu=0,...){
    if (minocc < 1.0){
    	minocc <- round(dim(obj)[1]*minocc, 0)
    }
    obj <- obj[,apply(obj>minabu,2,sum)>=minocc]
    return (obj)
})

#' @aliases droptax,phyloseq
#' @rdname droptax
#' @importFrom phyloseq otu_table
#' @importFrom ape drop.tip
#' @export
setMethod("droptax", "phyloseq", function(obj, ...){
    otuda <- checkotu(obj)
    otuda <- droptax(otuda, ...)
    keepotu <- colnames(otuda)
    if (!is.null(obj@phy_tree)){
        tmptree <- obj@phy_tree
        removetip <- setdiff(tmptree$tip.label, keepotu)
        tmptree <- drop.tip(tmptree, removetip)
        obj@phy_tree <- tmptree
    }
    if (!is.null(obj@refseq)){
        tmpseq <- obj@refseq
        keepseq <- intersect(keepotu, names(tmpseq))
        tmpseq <- tmpseq[keepseq]
        obj@refseq <- tmpseq
    }
    if (!is.null(obj@tax_table)){
        tmptax <- obj@tax_table
        keeptax <- intersect(keepotu, rownames(tmptax))
        tmptax <- tmptax[match(keeptax,rownames(tmptax)),]
        obj@tax_table <- tmptax
    }
    obj@otu_table <- otu_table(otuda,taxa_are_rows=FALSE)
    return(obj)
})


#' @title Import function to load the feature table and taxanomy table of dada2
#'
#' @description
#' the function can import the ouput of dada2, and generated the phyloseq obj contained the
#' argument class.
#' @param seqtab matrix, feature table, the output of \code{\link[dada2]{removeBimeraDenovo}}.
#' @param taxatab matrix, a taxonomic table, the output of \code{\link[dada2]{assignTaxonomy}},
#' or the ouput of \code{\link[dada2]{addSpecies}}.
#' @param reftree phylo or character, the phylo class of tree, or the tree file.
#' @param sampleda data.frame or character, the data.frame of sample information,
#' or the file of sample information, nrow samples X ncol factors.
#' @param btree logical, whether building the tree, default is FALSE. 
#' @param ..., additional parameters, see also \code{\link{build_tree}}.
#' @importFrom phyloseq phyloseq otu_table sample_data tax_table phy_tree import_qiime_sample_data
#' @importFrom Biostrings DNAStringSet
#' @importFrom ape read.tree
#' @return phyloseq class contained the argument class.
#' @author Shuangbin Xu
#' @export
#' @examples
#' seqtabfile <- system.file("extdata", "seqtab.nochim.rds",
#'                           package="MicrobiotaProcess")
#' taxafile <- system.file("extdata", "taxa_tab.rds",
#'                         package="MicrobiotaProcess")
#' seqtab <- readRDS(seqtabfile)
#' taxa <- readRDS(taxafile)
#' sampleda <- system.file("extdata", "mouse.time.dada2.txt", 
#'                         package="MicrobiotaProcess")
#' ps <- import_dada2(seqtab=seqtab, taxatab=taxa,
#'                    sampleda=sampleda)
#' ps
import_dada2 <- function(seqtab, taxatab=NULL, reftree=NULL, 
                         sampleda=NULL, btree=FALSE, ...){
    refseq <- colnames(seqtab)
    refseqnm <- paste0("OTU_",seq_len(length(refseq)))
    colnames(seqtab) <- refseqnm
    names(refseq) <- refseqnm
    if (!is.null(taxatab)){
        if (!identical(colnames(seqtab), rownames(taxatab))){
            taxatab <- taxatab[match(refseq, rownames(taxatab)),,drop=FALSE]
        }
        rownames(taxatab) <- refseqnm
        taxatab <- tax_table(as.matrix(taxatab))
    }
    refseq <- DNAStringSet(refseq)
    if (is.null(reftree)){
        if (btree){
            reftree <- build_tree(refseq, ...)
        }
    }
    if (!is.null(reftree) & inherits(reftree,"character")){
        reftree <- read.tree(reftree)
    }
    if (!is.null(sampleda)){
        if (inherits(sampleda, "data.frame")){
            sampleda <- sample_data(sampleda)
        }
        if (inherits(sampleda, "character")){
            sampleda <- import_qiime_sample_data(sampleda)
        }
    }
    otuda <- otu_table(seqtab, taxa_are_rows=FALSE)
    arglist <- list(otuda, taxatab, sampleda, refseq, reftree)
    arglist <- arglist[!unlist(lapply(arglist,function(x)is.null(x)))]
    ps <- do.call("phyloseq", arglist)
    return(ps)
}

#' @title building tree
#'
#' @description
#' The function can be used to building tree.
#' @param seqs DNAStringSet or DNAbin, the object of R. 
#' @param ..., additional parameters, see also \code{\link[DECIPHER]{AlignSeqs}}.
#' @return the phylo class of tree.
#' @author Shuangbin Xu
#' @export
#' @examples
#' seqtabfile <- system.file("extdata", "seqtab.nochim.rds", 
#'                           package="MicrobiotaProcess")
#' seqtab <- readRDS(seqtabfile)
#' refseq <- colnames(seqtab)
#' names(refseq) <- paste0("OTU_",seq_len(length(refseq)))
#' # refseq <- Biostrings::DNAStringSet(refseq)
#' # tree <- build_tree(refseq)
#' # or
#' # tree <- build_tree(refseq) 
setGeneric("build_tree", function(seqs, ...){standardGeneric("build_tree")})

#' @aliases build_tree,DNAStringSet
#' @rdname build_tree
#' @importFrom Biostrings DNAStringSet
#' @importFrom DECIPHER AlignSeqs
#' @importFrom phangorn phyDat dist.ml NJ pml optim.pml pml.control
#' @importFrom stats update
#' @export
setMethod("build_tree", "DNAStringSet", function(seqs, ...){
    alignment <- AlignSeqs(seqs, anchor=NA, ...)
    phalign <- phyDat(as(alignment, "matrix"), type="DNA")
    dm <- dist.ml(phalign)
    treeNJ <- NJ(dm)
    fit <- pml(treeNJ, data=phalign)
    treeGTR <- update(fit, k=4, inv=0.2)
    tree <- optim.pml(treeGTR, 
                      model="GTR", 
                      optInv=TRUE, 
                      optGamma=TRUE,
                      rearrangement = "stochastic", 
                      control = pml.control(trace = 0))
    tree <- tree$tree
    return(tree)
})

#' @aliases build_tree,DNAbin
#' @rdname build_tree
#' @export
setMethod("build_tree", "DNAbin", function(seqs, ...){
    seqs <- DNAbin2DNASet(seqs)
    tree <- build_tree(seqs, ...)
    return(tree)
})

#' @aliases build_tree,character
#' @rdname build_tree
#' @importFrom Biostrings readDNAStringSet DNAStringSet
#' @export
setMethod("build_tree", "character", function(seqs,...){
    if (length(seqs)==1){
        seqs <- readDNAStringSet(seqs)
    }else{
        seqs <- DNAStringSet(seqs)
    }
    tree <- build_tree(seqs, ...)
    return (tree)
})

#' @importFrom Biostrings DNAStringSet
#' @keywords internal
DNAbin2DNASet <- function(seqs){
    seqs <- toupper(unlist(lapply(as.character(seqs), paste0, collapse="")))
    seqs <- DNAStringSet(seqs)
    return (seqs)
}


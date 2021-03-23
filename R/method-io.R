#' @title Import function to load the output of qiime2.
#'
#' @description
#' The function was designed to import the output of qiime2 and convert them to phyloseq
#' class.
#' @param otuqza character, the file contained otu table, the ouput of qiime2.
#' @param taxaqza character, the file contained taxonomy, the ouput of qiime2,
#' default is NULL.
#' @param mapfilename character, the file contained sample information,
#' the tsv format, default is NULL.
#' @param refseqqza character, the file contained refrentent sequences, default is
#' NULL.
#' @param treeqza character, the file contained the tree file, default is NULL.
#' @param build_tree logical, whether building the tree, when the rownames of 
#'  feature table contains the sequence, default is FALSE.
#' @param parallel logical, whether parsing the column of taxonomy multi-parallel, default is FALSE.
#' @param ..., additional parameters, see alos\code{\link{build_tree}}.
#' @return phyloseq-class contained the argument class.
#' @importFrom phyloseq import_qiime_sample_data tax_table otu_table phyloseq
#' @importFrom Biostrings DNAStringSet
#' @export
#' @author Shuangbin Xu
#' @examples
#' otuqzafile <- system.file("extdata", "table.qza",
#'                           package="MicrobiotaProcess")
#' taxaqzafile <- system.file("extdata", "taxa.qza",
#'                            package="MicrobiotaProcess")
#' mapfile <- system.file("extdata", "metadata_qza.txt",
#'                        package="MicrobiotaProcess")
#' ps <- import_qiime2(otuqza=otuqzafile, taxaqza=taxaqzafile,
#'                     mapfilename=mapfile)
#' ps
import_qiime2 <- function(otuqza, taxaqza=NULL, mapfilename=NULL, 
                          refseqqza=NULL, treeqza=NULL,
                          build_tree=FALSE, parallel=FALSE, ...){
    otux <- read_qza(otuqza, parallel=parallel)
    otutab <- otu_table(otux$otutab,taxa_are_rows=TRUE)
    refseq <- reftree <- taxtab <- NULL
    if (!is.null(taxaqza)){
        taxax <- read_qza(taxaqza, parallel=parallel)
        taxax <- taxax[match(rownames(otutab), rownames(taxax)),,drop=FALSE]
        if (!is.null(otux$taxtab)){
            taxax <- otux$taxtab
        }
        flag <- guess_rownames(rownames(otutab))
        if (flag=="DNA"){
            refseq <- rownames(otutab)
            refseqnm <- paste0("OTU_", seq_len(length(refseq)))
            rownames(otutab) <- refseqnm
            rownames(taxax) <- refseqnm
            names(refseq) <- refseqnm
            if (build_tree){
                reftree <- build_tree(DNAStringSet(refseq), ...)
            }
        }
        taxtab <- tax_table(taxax)
    }
    if (!is.null(mapfilename)){
        sampleda <- import_qiime_sample_data(mapfilename)
    }
    if (!is.null(refseqqza)){
        refseq <- read_qza(refseqqza)
    }
    if (!is.null(refseq)){
        refseq <- DNAStringSet(refseq)
    }
    if (!is.null(treeqza)){
        reftree <- read_qza(treeqza)
    }
    if (is.null(mapfilename)){
        sampleda <- NULL
    }
    arglist <- list(otutab, sampleda, taxtab, refseq, reftree)
    arglist <- arglist[!unlist(lapply(arglist, function(x)is.null(x)))]
    ps <- do.call("phyloseq", arglist)
    return(ps)
}

#' @title read the qza file, output of qiime2.
#' @description
#' the function was designed to read the ouput of qiime2.
#' @param qzafile character, the format of file should be one of
#' `BIOMV210DirFmt`, `TSVTaxonomyDirectoryFormat`, `NewickDirectoryFormat`
#' and `DNASequencesDirectoryFormat`.
#' @param parallel logical, whether parsing the taxonomy by multi-parallel, 
#' efault is FALSE.
#' @return list contained one or multiple object of feature table, 
#' taxonomy table, tree and represent sequences.
#' @importFrom ape read.tree
#' @importFrom Biostrings readDNAStringSet
#' @importFrom utils unzip
#' @export
#' @examples
#' otuqzafile <- system.file("extdata", "table.qza",
#'                           package="MicrobiotaProcess")
#' otuqza <- read_qza(otuqzafile)
#' str(otuqza)
read_qza <- function(qzafile, parallel=FALSE){
    tmpdir <- tempdir()
    unzipfiles <- unzip(qzafile, exdir=tmpdir)
    metadafile <- unzipfiles[grep("metadata.yaml", unzipfiles)[1]]
    metaflag <- yaml::read_yaml(metadafile)
    formatflag <- metaflag$format
    datafile <- unzipfiles[3]
    formats <- c("BIOMV210DirFmt", "TSVTaxonomyDirectoryFormat", 
                 "NewickDirectoryFormat","DNASequencesDirectoryFormat")
    formatflag <- match.arg(formatflag, formats)
    switch(formatflag,
           BIOMV210DirFmt={
               datafile <- unzipfiles[grep(".biom", unzipfiles)]
               x <- read.featuretab(datafile)},
           TSVTaxonomyDirectoryFormat={
               datafile <- unzipfiles[grep("data/taxonomy.tsv", unzipfiles)]
               x <- read.taxa(datafile, parallel=parallel)},
           DNASequencesDirectoryFormat={
               datafile <- unzipfiles[grep("data/.*\\.fasta", unzipfiles)]
               x <- readDNAStringSet(datafile)},
           NewickDirectoryFormat = {
               datafile <- unzipfiles[grep("data/tree.nwk", unzipfiles)]
               x <- read.tree(datafile)})
    return(x)
}

#' @importFrom phyloseq parse_taxonomy_greengenes
#' @keywords internal
read.featuretab <- function(file){
    biomobj <- suppressWarnings(biomformat::read_biom(file))
    x <- data.frame(as(biomformat::biom_data(biomobj),"matrix"), check.names=FALSE)
    taxflag <- all(unlist(lapply(biomobj$rows, function(i){length(i$metadata)}))==0)
    if (taxflag){
        taxtab <- NULL
    }else{
        taxtab <- lapply(biomobj$rows, function(i){parse_taxonomy_greengenes(i$metadata$taxonomy)})
        names(taxtab) <- rownames(x)
        taxtab <- build_tax_table(taxtab)
    }
    return(list(otutab = x, taxtab=taxtab))
}

#' @importFrom plyr llply
#' @importFrom phyloseq parse_taxonomy_qiime build_tax_table 
#' @importFrom utils read.table
#' @keywords internal 
read.taxa <- function(file, parallel=FALSE){
    x <- read.table(file, sep="\t", row.names=1, header=TRUE)
    taxstring <- as.vector(x[[1]])
    taxtab <- suppressWarnings(llply(taxstring, parse_taxonomy_qiime, .parallel = parallel))
    names(taxtab) <- rownames(x)
    taxtab <- build_tax_table(taxtab)
    return(taxtab)
}

#' @keywords internal
guess_rownames <- function(names){
    x <- unlist(strsplit(toupper(names[1]), split = ""))
    freq <- mean(x %in% c("A", "G","T","C", "N", "X", "-"))
    if (length(x) > 30 & freq > 0.9){
        return("DNA")
    }else{
        return("Others")    
    }
}

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
#' @importFrom stats update
#' @export
setMethod("build_tree", "DNAStringSet", function(seqs, ...){
    alignment <- DECIPHER::AlignSeqs(seqs, anchor=NA, ...)
    phalign <- phangorn::phyDat(as(alignment, "matrix"), type="DNA")
    dm <- phangorn::dist.ml(phalign)
    treeNJ <- phangorn::NJ(dm)
    fit <- phangorn::pml(treeNJ, data=phalign)
    treeGTR <- update(fit, k=4, inv=0.2)
    tree <- phangorn::optim.pml(treeGTR, 
                      model="GTR", 
                      optInv=TRUE, 
                      optGamma=TRUE,
                      rearrangement = "stochastic", 
                      control = phangorn::pml.control(trace = 0))
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


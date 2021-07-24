#' @title Import function to load the output of qiime2.
#' 
#' @description
#' The function was designed to import the output of qiime2 and convert them to phyloseq
#' class.
#' @name ImportQiime2
#' @param otuqza character, the file contained otu table, the ouput of qiime2.
#' @param taxaqza character, the file contained taxonomy, the ouput of qiime2,
#' default is NULL.
#' @param mapfilename character, the file contained sample information,
#' the tsv format, default is NULL.
#' @param refseqqza character, the file contained reference sequences or the XStringSet object,
#' default is NULL.
#' @param treeqza character, the file contained the tree file or treedata object, which is the result 
#' by parsing function of treeio, default is NULL.
#' @param parallel logical, whether parsing the column of taxonomy multi-parallel, default is FALSE.
#' @param ..., additional parameters.
#' @return MPSE-class or phyloseq-class contained the argument class.
#' @importFrom phyloseq phyloseq
#' @export
#' @author Shuangbin Xu
#' @examples
#' otuqzafile <- system.file("extdata", "table.qza",
#'                           package="MicrobiotaProcess")
#' taxaqzafile <- system.file("extdata", "taxa.qza",
#'                            package="MicrobiotaProcess")
#' mapfile <- system.file("extdata", "metadata_qza.txt",
#'                        package="MicrobiotaProcess")
#' mpse <- mp_import_qiime2(otuqza=otuqzafile, taxaqza=taxaqzafile,
#'                          mapfilename=mapfile)
#' mpse
import_qiime2 <- function(otuqza, taxaqza=NULL, mapfilename=NULL, 
                          refseqqza=NULL, treeqza=NULL,
                          parallel=FALSE, ...){
    res <- .internal_import_qiime2(otuqza=otuqza, 
                                   taxaqza=taxaqza, 
                                   mapfilename=mapfilename, 
                                   refseqqza=refseqqza, 
                                   treeqza=treeqza, 
                                   parallel=parallel, 
                                   ...)
    ps <- .build_ps(res)
}

#' @rdname ImportQiime2
#' @export
mp_import_qiime2 <- function(otuqza, taxaqza=NULL, mapfilename=NULL,
                          refseqqza=NULL, treeqza=NULL, parallel=FALSE, ...){

    res <- .internal_import_qiime2(otuqza=otuqza,
                                   taxaqza=taxaqza,
                                   mapfilename=mapfilename,
                                   refseqqza=refseqqza,
                                   treeqza=treeqza,
                                   parallel=parallel, 
                                   ...)
    mpse <- .build_mpse(res=res)
    return(mpse)
}

#' @title Import function to load the feature table and taxanomy table of dada2
#'
#' @description
#' the function can import the ouput of dada2, and generated the phyloseq obj contained the
#' argument class.
#' @name ImportDada2
#' @param seqtab matrix, feature table, the output of \code{\link[dada2]{removeBimeraDenovo}}.
#' @param taxatab matrix, a taxonomic table, the output of \code{\link[dada2]{assignTaxonomy}},
#' or the ouput of \code{\link[dada2]{addSpecies}}.
#' @param reftree phylo, treedata or character, the treedata or phylo class of tree, or the tree file.
#' @param sampleda data.frame or character, the data.frame of sample information,
#' or the file of sample information, nrow samples X ncol factors.
#' @param ..., additional parameters.
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
#' mpse <- mp_import_dada2(seqtab=seqtab, taxatab=taxa,
#'                    sampleda=sampleda)
#' mpse
import_dada2 <- function(seqtab, taxatab=NULL, reftree=NULL, 
                         sampleda=NULL, 
                         ...){

    res <- .internal_import_dada2(seqtab=seqtab,
                                  taxatab=taxatab,
                                  reftree=reftree,
                                  sampleda=sampleda)

    ps <- .build_ps(res)

    return(ps)
}

#' @rdname ImportDada2
#' @export
mp_import_dada2 <- function(seqtab,
                             taxatab=NULL,
                             reftree=NULL,
                             sampleda=NULL,
                             ...){

    res <- .internal_import_dada2(seqtab=seqtab,
                                  taxatab=taxatab,
                                  reftree=reftree,
                                  sampleda=sampleda)

    mpse <- .build_mpse(res)

    return(mpse)
}

#' @title Import function to load the output of qiime.
#'
#' @description
#' The function was designed to import the output of qiime and convert them to MPSE
#' class.
#' @param otufilename character, the file contained otu table, the ouput of qiime.
#' @param mapfilename character, the file contained sample information,
#' the tsv format, default is NULL.
#' @param otutree treedata, phylo or character, the file contained reference sequences, or
#' treedata object, which is the result by parsing function of treeio, default is NULL.
#' @param refseq XStringSet or character, the file contained the tree file or XStringSet, default is NULL.
#' @param ..., additional parameters.
#' @return MPSE-class.
#' @export
#' @author Shuangbin Xu
mp_import_qiime <- function(otufilename, 
                            mapfilename = NULL,
                            otutree = NULL, 
                            refseq = NULL,
                            ...){

    res <- .internal_import_qiime(
                                  otufilename = otufilename,
                                  mapfilename = mapfilename,
                                  otutree = otutree,
                                  refseq = refseq,
                                  ...
                                  )

    mpse <- .build_mpse(res)

    return(mpse)
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
               x <- Biostrings::readDNAStringSet(datafile)},
           NewickDirectoryFormat = {
               datafile <- unzipfiles[grep("data/tree.nwk", unzipfiles)]
               x <- read.tree(datafile)})
    return(x)
}

#' @keywords internal
read.featuretab <- function(file){
    biomobj <- suppressWarnings(biomformat::read_biom(file))
    x <- data.frame(as(biomformat::biom_data(biomobj),"matrix"), check.names=FALSE)
    taxflag <- all(unlist(lapply(biomobj$rows, function(i){length(i$metadata)}))==0)
    if (taxflag){
        taxatab <- NULL
    }else{
        taxatab <- lapply(biomobj$rows, function(i){phyloseq::parse_taxonomy_greengenes(i$metadata$taxonomy)})
        names(taxatab) <- rownames(x)
        taxatab <- phyloseq::build_tax_table(taxatab)
        #taxatab <- remove_na_taxonomy_rank(taxatab)
        taxatab <- fillNAtax(taxatab)
    }
    return(list(otutab = x, taxatab=taxatab))
}

#' @keywords internal
read.taxa <- function(file, parallel=FALSE){
    x <- utils::read.table(file, sep="\t", row.names=1, header=TRUE, comment.char="")
    taxstring <- as.vector(x[[1]])
    taxatab <- suppressWarnings(plyr::llply(taxstring, phyloseq::parse_taxonomy_qiime, .parallel = parallel))
    names(taxatab) <- rownames(x)
    taxatab <- phyloseq::build_tax_table(taxatab)
    #taxatab <- remove_na_taxonomy_rank(taxatab)
    taxatab <- fillNAtax(taxatab)
    return(taxatab)
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
#' @export
setMethod("build_tree", "character", function(seqs,...){
    if (length(seqs)==1){
        seqs <- Biostrings::readDNAStringSet(seqs)
    }else{
        seqs <- Biostrings::DNAStringSet(seqs)
    }
    tree <- build_tree(seqs, ...)
    return (tree)
})

#' @keywords internal
DNAbin2DNASet <- function(seqs){
    seqs <- toupper(unlist(lapply(as.character(seqs), paste0, collapse="")))
    seqs <- Biostrings::DNAStringSet(seqs)
    return (seqs)
}


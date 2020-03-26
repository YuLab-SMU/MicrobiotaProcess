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
    otux <- read.qza(otuqza, parallel=parallel)
    otutab <- otu_table(otux$otutab,taxa_are_rows=TRUE)
    if (!is.null(taxaqza)){
        taxax <- read.qza(taxaqza, parallel=parallel)
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
            }else{
                reftree <- NULL
            }
        }else{
            refseq <- NULL
            reftree <- NULL
        }
        taxtab <- tax_table(taxax)
    }
    if (!is.null(mapfilename)){
        sampleda <- import_qiime_sample_data(mapfilename)
    }
    if (!is.null(refseqqza)){
        refseq <- read.qza(refseqqza)
    }
    if (!is.null(refseq)){
        refseq <- DNAStringSet(refseq)
    }
    if (!is.null(treeqza)){
        reftree <- read.qza(treeqza)
    }
    if (is.null(mapfilename)){
        sampleda <- NULL
    }
    if (is.null(taxaqza) && is.null(taxtab)){
        taxtab <- NULL
    }
    if (is.null(refseqqza) && is.null(refseq)){
        refseq <- NULL
    }
    if (is.null(treeqza) && is.null(reftree)){
        reftree <- NULL
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
#' @importFrom yaml read_yaml
#' @importFrom ape read.tree
#' @importFrom Biostrings readDNAStringSet
#' @importFrom utils unzip
#' @export
#' @examples
#' otuqzafile <- system.file("extdata", "table.qza",
#'                           package="MicrobiotaProcess")
#' otuqza <- read.qza(otuqzafile)
#' str(otuqza)
read.qza <- function(qzafile, parallel=FALSE){
    tmpdir <- tempdir()
    unzipfiles <- unzip(qzafile, exdir=tmpdir)
    metadafile <- unzipfiles[grep("metadata.yaml", unzipfiles)[1]]
    metaflag <- read_yaml(metadafile)
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

#' @importFrom biomformat read_biom biom_data
#' @importFrom phyloseq parse_taxonomy_greengenes
#' @keywords internal
read.featuretab <- function(file){
    biomobj <- suppressWarnings(read_biom(file))
    x <- data.frame(as(biom_data(biomobj),"matrix"), check.names=FALSE)
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


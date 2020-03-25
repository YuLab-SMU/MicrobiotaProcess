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
#' @param build_tree logical, whether building the tree, default is FALSE.
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
    otux <- read.qza(otuqza, build_tree=build_tree, parallel=parallel, ...)
    otutab <- otu_table(otux$otutab,taxa_are_rows=TRUE)
    if (!is.null(taxaqza)){
        taxax <- read.qza(taxaqza, build_tree=FALSE, parallel=parallel)
        if (!is.null(taxax$refseq) & !is.null(otux$refseq)){
            matchnames <- names(taxax$refseq[match(otux$refseq, taxax$refseq)])
            rownames(taxax$taxtab[match(matchnames, rownames(taxax$taxtab)),]) <- names(otux$refseq)
        }
        taxtab <- tax_table(taxax$taxtab)
    }
    if (!is.null(otux$taxtab)){
        taxtab <- tax_table(otux$taxtab)
    }
    if (!is.null(mapfilename)){
        sampleda <- import_qiime_sample_data(mapfilename)
    }
    if (!is.null(refseqqza)){
        refseq <- read.qza(refseqqza)
    }
    if (!is.null(otux$refseq)){
        refseq <- DNAStringSet(otux$refseq)
    }
    if (!is.null(otux$reftree)){
        reftree <- otux$reftree
    }
    if (!is.null(treeqza)){
        reftree <- read.qza(treeqza)
    }
    if (is.null(mapfilename)){
        sampleda <- NULL
    }
    if (is.null(taxaqza) & is.null(otux$taxtab)){
        taxtab <- NULL
    }
    if (is.null(refseqqza) & is.null(otux$refseq)){
        refseq <- NULL
    }
    if (is.null(treeqza) & is.null(otux$reftree)){
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
#' @param build_tree logical, whether building the tree, default is FALSE.
#' @param parallel logical, whether parsing the taxonomy by multi-parallel, 
#' efault is FALSE.
#' @param ..., additional parameters, see also \code{\link{build_tree}}.
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
read.qza <- function(qzafile, build_tree=FALSE, parallel=FALSE, ...){
    tmpdir <- tempdir()
    unzipfiles <- unzip(qzafile, exdir=tmpdir)
    tmp <- grep("yaml", unzipfiles)[1]
    metaflag <- read_yaml(unzipfiles[tmp])
    formatflag <- metaflag$format
    formats <- c("BIOMV210DirFmt", "TSVTaxonomyDirectoryFormat", 
                 "NewickDirectoryFormat","DNASequencesDirectoryFormat")
    formatflag <- match.arg(formatflag, formats)
    switch(formatflag,
           BIOMV210DirFmt={x <- grep("biom", unzipfiles)},
           TSVTaxonomyDirectoryFormat={x <- grep("tsv", unzipfiles)},
           DNASequencesDirectoryFormat={x <- grep("fasta", unzipfiles)},
           NewickDirectoryFormat = {x <- grep("nwk", unzipfiles)})
    datafile <- unzipfiles[x]
    switch(formatflag,
           BIOMV210DirFmt={x <- read.featuretab(datafile, build_tree=build_tree, ...)},
           TSVTaxonomyDirectoryFormat={x <- read.taxa(datafile, parallel=parallel)},
           DNASequencesDirectoryFormat={x <- readDNAStringSet(datafile)},
           NewickDirectoryFormat = {x <- read.tree(datafile)})
    return(x)
}

#' @importFrom biomformat read_biom biom_data
#' @importFrom phyloseq parse_taxonomy_greengenes
#' @keywords internal
read.featuretab <- function(file, build_tree=FALSE, ...){
    biomobj <- suppressWarnings(read_biom(file))
    x <- data.frame(as(biom_data(biomobj),"matrix"), check.names=FALSE)
    taxflag <- all(unlist(lapply(biomobj$rows, function(i){length(i$metadata)}))==0)
    flag <- guess_rownames(rownames(x))
    if (flag=="DNA"){
        refseq <- rownames(x)
        refseqnm <- paste0("OTU_", seq_len(length(refseq)))
        rownames(x) <- refseqnm
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
    if (taxflag){
        taxtab <- NULL
    }else{
        taxtab <- lapply(biomobj$rows, function(i){parse_taxonomy_greengenes(i$metadata$taxonomy)})
        names(taxtab) <- rownames(x)
        taxtab <- build_tax_table(taxtab)
    }
    return(list(otutab = x, taxtab=taxtab, refseq = refseq, reftree= reftree))
}

#' @importFrom plyr llply
#' @importFrom phyloseq parse_taxonomy_qiime build_tax_table 
#' @importFrom utils read.table
#' @keywords internal 
read.taxa <- function(file, parallel=FALSE){
    x <- read.table(file, sep="\t", row.names=1, header=TRUE)
    flag <- guess_rownames(rownames(x))
    if (flag=="DNA"){
        refseq <- rownames(x)
        refseqnm <- paste0("OTU_", seq_len(length(refseq)))
        rownames(x) <- refseqnm
        names(refseq) <- refseqnm
    }else{
        refseq <- NULL
    }
    taxstring <- as.vector(x[[1]])
    taxtab <- suppressWarnings(llply(taxstring, parse_taxonomy_qiime, .parallel = parallel))
    names(taxtab) <- rownames(x)
    taxtab <- build_tax_table(taxtab)
    return(list(taxtab=taxtab, refseq=refseq))
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


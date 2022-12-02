#' @keywords internal 
.internal_import_qiime2 <- function(otuqza, taxaqza=NULL, mapfilename=NULL, 
                                    refseqqza=NULL, treeqza=NULL, 
                                    parallel=FALSE,
                                    ...){
    otux <- read_qza(otuqza, parallel=parallel)
    otutab <- otux$otutab
    sampleda <- otux$sampleda
    refseq <- reftree <- taxatab <- otu.metada <- NULL
    if (!is.null(taxaqza)){
        otu.taxa.meta <- read_qza(taxaqza, parallel=parallel)
        taxax <- otu.taxa.meta$taxatab
        otu.metada <- otu.taxa.meta$otu.metada
        taxax <- taxax[match(rownames(otutab), rownames(taxax)),,drop=FALSE]
        otu.metada <- otu.metada[match(rownames(otutab), rownames(otu.metada)),,drop=FALSE]
        if (!is.null(otux$taxatab)){
            taxax <- otux$taxatab
            otu.metada <- NULL
        }
        flag <- guess_rownames(rownames(otutab))
        if (flag!="Other"){
            refseq <- rownames(otutab)
            refseqnm <- paste0("OTU_", seq_len(length(refseq)))
            rownames(otutab) <- refseqnm
            rownames(taxax) <- refseqnm
            names(refseq) <- refseqnm
        }
        taxatab <- as.matrix(taxax)
    }
    if (!is.null(mapfilename)){
        if (inherits(mapfilename, "character") && file.exists(mapfilename)){
            sampleda <- read_qiime_sample(mapfilename)
        }else if (!inherits(mapfilename, "data.frame")){
            rlang::abort("The mapfilename should be a file or data.frame contained sample information!")
        }else{
            sampleda <- mapfilename
        }
    }
    if (!is.null(refseqqza)){
        if (inherits(refseqqza, "XStringSet")){
            refseq <- refseqqza
        }else if (inherits(refseqqza, "character") && file.exists(refseqqza)){
            refseq <- read_qza(refseqqza)
        }
    }
    if (!is.null(refseq)){
        refseq <- build_refseq(refseq)
    }
    if (!is.null(treeqza)){
        if (inherits(treeqza, "treedata")){
            reftree <- treeqza
        }else if (inherits(treeqza, "character") && file.exists(treeqza)){
            rlang::warn(c("The treeqza is a tree file, it is parsing by read.tree function.",
                     "It is better to parse it with the function of treeio",
                     "then the treedata or phylo result all can be supported."))
            reftree <- treeio::read.treeqza(treeqza) %>% 
                       treeio::as.treedata()
        }else if (inherits(treeqza, "phylo")){
            reftree <- treeio::as.treedata(treeqza)
        }
    }

    return (list(otutab=otutab, 
                sampleda=sampleda, 
                taxatab=taxatab, 
                otu.metada = otu.metada,
                refseq=refseq, 
                otutree=reftree))
}

.internal_import_qiime <- function(otufilename, mapfilename = NULL, otutree = NULL, refseq = NULL){
    x <- read_qiime_otu(otufilename)

    if (!is.null(mapfilename)){
        if (file.exists(mapfilename)){
            sampleda <- read_qiime_sample(mapfilename)
        }else if (!inherits(mapfilename, "data.frame")){
            sampleda <- NULL
        }else{
            sampleda <- mapfilename
        }
    }else{
        sampleda <- NULL
    }

    if (!is.null(otutree)){
        if (file.exists(otutree)){
            rlang::warn(c("The reftree is a tree file, it is parsing by read.tree function.",
                     "It is better to parse it with the function of treeio",
                     "then the treedata or phylo result all can be supported."))
            otutree <- ape::read.tree(otutree) %>% treeio::as.treedata()
        }else if (inherits(otutree, "phylo")){
            otutree <- otutree %>% treeio::as.treedata()
        }else if (!inherits(otutree, "treedata")){
            otutree <- NULL
        }
    }

    if (!is.null(refseq)){
        if (file.exists(refseq)){
            refseq <- seqmagick::fa_read(refseq)
        }else if (inherits(refseq, "character")){
            refseq <- build_tree(refseq)
        }else if (!inherits(refseq, "XStringSet")){
            refseq <- NULL
        }
    }

    return (list(otutab = x$otuda,
                 sampleda = sampleda,
                 taxatab = x$taxada,
                 refseq = refseq,
                 otutree = otutree))
}

#' @keywords internal
.internal_import_dada2 <- function(seqtab, 
                                   taxatab=NULL, 
                                   reftree=NULL,
                                   sampleda=NULL, 
                                   ...){
    refseq <- colnames(seqtab)
    refseqnm <- paste0("OTU_", seq_len(length(refseq)))
    colnames(seqtab) <- refseqnm
    names(refseq) <- refseqnm
    if (!is.null(taxatab)){
        if (!identical(colnames(seqtab), rownames(taxatab))){
            taxatab <- taxatab[match(refseq, rownames(taxatab)), , drop=FALSE]
        }
        rownames(taxatab) <- refseqnm
        taxatab <- fillNAtax(taxatab) 
    }
    if (!is.null(reftree)){
        if (inherits(reftree, "phylo")){
            reftree <- .check_tree_tiplab(reftree, refseq, refseqnm)
        }else if(inherits(reftree, "character") && file.exists(reftree)){
            rlang::warn(c("The reftree is a tree file, it is parsing by read.tree function.",
                     "It is better to parse it with the function of treeio", 
                     "then the treedata or phylo result all can be supported."))
            reftree <- ape::read.tree(reftree)
            reftree <- .check_tree_tiplab(reftree, refseq, refseqnm)
        }
        reftree <- treeio::as.treedata(reftree)
    }
    if (!is.null(sampleda)){
        if (inherits(sampleda, "character") && file.exists(sampleda)){
           sampleda <- read_qiime_sample(sampleda)
        }else if (!inherits(sampleda, "data.frame")){
           rlang::abort("The sampleda should be a file or data.frame contained sample information!")
        }
    }
    refseq <- build_refseq(refseq)
    return (list(otutab=data.frame(t(seqtab),check.names=FALSE),
                sampleda=sampleda,
                taxatab=taxatab,
                refseq=refseq,
                otutree=reftree))
}


.check_tree_tiplab <- function(tree, oldnm, newnm){
    if (any(tree$tip.label %in% oldnm)){
        tree$tip.label <- newnm[match(tree$tip.label, oldnm)]
    }
    return(tree)
}


.build_mpse <- function(res){
    
    otuda <- res$otutab
    sampleda <- res$sampleda
    taxatab <- res$taxatab
    otutree <- res$otutree
    refseq <- res$refseq
    otu.metada <- res$otu.metada

    otuda <- otuda[rowSums(otuda) > 0, colSums(otuda) > 0, drop=FALSE]
    #otuda <- otuda[,colSums(otuda)>0, drop=FALSE]
    
    rownm <- rownames(otuda)
    colnm <- colnames(otuda)
    
    if (is.null(sampleda)){
        sampleda <- S4Vectors::DataFrame(row.names=colnames(otuda))
    }else{
        rnm <- rownames(sampleda)
        flagn <- intersect(rnm, colnm)
        if(!check_equal_len(list(rnm, flagn, colnm))){
            rlang::warn(c("The number of samples in otu table is not equal the number of samples in sample data.", 
                     "The same samples will be extract automatically !"))
        
            sampleda <- sampleda[match(flagn, rnm), , drop = FALSE]
            otuda <- otuda[, match(flagn, colnm), drop=FALSE]
        }else{
            sampleda <- sampleda[match(colnames(otuda), rownames(sampleda)),,drop=FALSE]
        }
    }

    if (!is.null(taxatab)){
        rnm <- rownames(taxatab)
        flagn <- intersect(rnm, rownm)
        if(! check_equal_len(list(rnm, flagn, rownm))){
            rlang::warn(c("The number of features in otu table is not equal the number of features in taxonomy annotation.",
                          "The same features will be extract automatically !"))
        
            taxatab <- taxatab[match(flagn, rnm), , drop = FALSE]
            otuda <- otuda[match(flagn, rownm), , drop = FALSE]
            if (!is.null(otu.metada)){
                otu.metada <- otu.metada[match(flagn, rownames(otu.metada)),,drop=FALSE]
            }
            rownm <- rownames(otuda)
        }
        taxatree <- try_convert_taxa(data.frame(taxatab))
    }else{
        taxatree <- NULL
    }

    if (!is.null(otutree)){
        tip.label <- otutree@phylo$tip.label
        flagn <- intersect(otutree@phylo$tip.label, rownm)
        if ( !check_equal_len(list(flagn, otutree@phylo$tip.label, rownm))){
            rlang::warn(c("The number of features in otu table is not equal the number of tips in otu tree.",
                          "The same features will be extract automatically !"))
            rmotus <- setdiff(tip.label, flagn)
        }else{
            rmotus <- NULL
        }
        if (length(rmotus) > 0){
            otutree <- treeio::drop.tip(otutree, tip=rmotus)
            otuda <- otuda[match(flagn, rownm), , drop = FALSE]
            rownm <- rownames(otuda)
            if (!is.null(taxatree)){
                taxatree <- treeio::drop.tip(taxatree, tip=setdiff(taxatree@phylo$tip.label, flagn), collapse.singles=FALSE)
            }else{
                taxatree <- NULL
            }
            if (!is.null(otu.metada)){
                otu.metada <- otu.metada[match(flagn, rownames(otu.metada)),,drop=FALSE]
            }
        }
    }

    if (!is.null(refseq)){
        rnm <- names(refseq)
        flagn <- intersect(rnm, rownm)
        if ( !check_equal_len(list(rnm, flagn, rownm))){
            rlang::warn(c("The number of features in otu table is not equal the number of sequence in refseq.",
                          "The same features will be extract automatically !"))
            rmotus <- setdiff(rownm, flagn)
        }else{
            rmotus <- NULL
        }
        refseq <- refseq[flagn]
        if (length(rmotus) > 0){
            otuda <- otuda[match(flagn, rownm), , drop = FALSE]
            if (!is.null(otutree) && length(rmotus) != treeio::Ntip(otutree)){
                otutree <- treeio::drop.tip(otutree, tip=rmotus)
            }else{
                otutree <- NULL
            }
            if (!is.null(taxatree) && length(rmotus) != treeio::Ntip(taxatree)){
                taxatree <- treeio::drop.tip(taxatree, tip=rmotus, collapse.singles=FALSE)
            }else{
                taxatree <- NULL
            }
            if (!is.null(otu.metada)){
                otu.metada <- otu.metada[match(flagn, rownames(otu.metada)),,drop=FALSE]
            }
        }
    }
    
    mpse <- MPSE(
                 assays = list(Abundance=otuda %>% as.matrix()),
                 colData = sampleda,
                 otutree = otutree,
                 taxatree =  taxatree,
                 refseq = refseq
            )

    if (!is.null(otu.metada) && ncol(otu.metada)>0){
        SummarizedExperiment::rowData(mpse) <- otu.metada
    }

    methods::validObject(mpse)
    return(mpse)
}

.build_ps <- function(res){
    if (!is.null(res$otutree)){
        reftree <- res$otutree@phylo
    }else{
        reftree <- NULL
    }
    ps <- phyloseq::phyloseq(
          phyloseq::otu_table(res$otutab,taxa_are_rows=TRUE),
          phyloseq::sample_data(res$sampleda, FALSE),
          phyloseq::tax_table(as.matrix(res$taxatab), FALSE),
          phyloseq::phy_tree(reftree, FALSE),
          phyloseq::refseq(res$refseq, FALSE)
      )
    return(ps)
}

read_qiime_otu <- function(otufilename){
    skipn <- guess_skip_nrow(otufilename)
    xx <- utils::read.table(otufilename, sep="\t", skip=skipn, header=TRUE, 
                            row.names=1, comment.char="", quote="")
    indy <- vapply(xx, function(i)is.numeric(i), logical(1)) %>% setNames(NULL)
    otuda <- xx[, indy, drop=FALSE]
    if (!all(indy)){
        taxada <- xx[, !indy, drop=FALSE]
        taxada <- apply(taxada, 2, function(i)gsub(";\\s+", ";", i)) %>% 
                  data.frame() %>%
                  split_str_to_list(sep=";") %>% 
                  fillNAtax() %>% 
                  setNames(taxaclass[seq_len(ncol(.))])
    }else{
        taxada <- NULL
    }
    return(list(otuda=otuda, taxada=taxada))
}

.internal_parse_biom <- function(biomobj){
    x <- data.frame(as.matrix(biomformat::biom_data(biomobj)), check.names=FALSE)
    taxflag <- unlist(lapply(biomobj$rows, function(i){length(i$metadata$taxonomy)}))
    sampleda <- lapply(biomobj$columns, function(i)c(Sample=i$id, i$metadata)) %>% 
                dplyr::bind_rows() %>% 
                tibble::column_to_rownames(var=colnames(.)[1])
    if (ncol(sampleda)==1){
        sampleda <- NULL
    }
    if (all(taxflag==0)){
        taxatab <- NULL
    }else if(all(taxflag == 1)){
        taxatab <- lapply(biomobj$rows, function(i)i$metadata$taxonomy) %>%
                   gsub(";\\s+", "@@",.) %>%
                   unlist() %>%
                   data.frame() %>%
                   split_str_to_list(sep="@@") %>%
                   apply(., 2, function(x)gsub("\\s+$", "", gsub("^\\s+", "", x))) %>%
                   data.frame() %>%
                   magrittr::set_colnames(taxaclass[seq_len(ncol(.))]) %>%
                   magrittr::set_rownames(rownames(x)) %>%
                   fillNAtax()
    }else{
        taxatab <- lapply(biomobj$rows, function(i)paste0(i$metadata$taxonomy, collapse="@@")) %>%
            unlist() %>%
            data.frame() %>%
            split_str_to_list(sep="@@") %>%
            dplyr::select(which(lapply(., 
                         function(i)!sum(grepl("[-]?[0-9]+[.]?[0-9]*|[-]?[0-9]+[L]?|[-]?[0-9]+[.]?[0-9]*[eE][0-9]+", i))> 0.5*nrow(.)) %>% 
                         unlist())) %>%
            dplyr::select(which(lapply(., function(i)!all(nchar(i)>60)) %>% unlist())) %>%
            apply(., 2, function(x)gsub("\\s+$", "", gsub("^\\s+", "", x))) %>%
            data.frame() %>%
            magrittr::set_colnames(taxaclass[seq_len(ncol(.))]) %>%
            magrittr::set_rownames(rownames(x)) %>%
            fillNAtax()
    }
    return(list(otutab = x, taxatab=taxatab, sampleda=sampleda))
}

# This function was refer to the biom of biomformat packages, but this packages
# was not maintained.
.internal_biom <- function(biomfilename){
    trash = try(silent = TRUE, expr = {
        x <- jsonlite::fromJSON(biomfilename, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
        }
    )
    if (inherits(trash, "try-error")) {
        trash = try(silent = TRUE, expr = {
            x <- biomformat::read_hdf5_biom(biomfilename)
        })
    }
    if (inherits(trash, "try-error")) {
        stop_wrap("Both attempts to read input file:\n", biomfilename,
              "either as JSON (BIOM-v1) or HDF5 (BIOM-v2). ", 
              "Check file path, file name, file itself, then try again.")
    }

    biom.type = c("OTU table", "Pathway table", "Function table",
                  "Ortholog table", "Gene table", "Metabolite table", "Taxon table")

    if (!x$type %in% biom.type){
       if (any(grepl(x$type, biom.type, ignore.case = TRUE))) {
           x$type <- biom.type[1]
       }
    }
    return(biomformat::biom(x))

}

taxaclass <- c("Kingdom", 
               "Phylum", 
               "Class", 
               "Order", 
               "Family", 
               "Genus", 
               "Speies", 
               "Rank8", 
               "Rank9", 
               "Rank10")

read_qiime_sample <- function(samplefile){
    sep <- guess_sep(samplefile)    
    sampleda <- utils::read.table(
                  samplefile, header=TRUE, sep=sep, row.names=1, 
                  comment.char="", quote=""
                )
    return(sampleda)
}

remove_na_taxonomy_rank <- function(x){
    x <- as.data.frame(x, check.names=FALSE)
    x <- remove_unclassfied(x) 
    indx <- vapply(x, 
                   function(i)all(is.na(i)),
                   FUN.VALUE=logical(1)) %>% 
            setNames(NULL)
    x <- x[, !indx, drop=FALSE]
    return(x)
}

build_refseq <- function(x){
    flag <- guess_rownames(x)
    refseq <- switch(flag,
                     DNA=Biostrings::DNAStringSet(x),
                     AA=Biostrings::AAStringSet(x),
                     RNA=Biostrings::RNAStringSet(x),
                     Other=NULL
                ) 
    return(refseq)
}

#' @keywords internal
#' this is refer to the seqmagick
guess_rownames <- function(string){
    a <- string
    a <- strsplit(toupper(a), split = "")[[1]]
    freq1 <- mean(a %in% c("A", "C", "G", "T", "X", "N", "-"))
    freq2 <- mean(a %in% c("A", "C", "G", "T", "N"))
    if (freq1 > 0.9 && freq2 > 0) {
        return("DNA")
    }
    freq3 <- mean(a %in% c("A", "C", "G", "U", "X", "N", "-"))
    freq4 <- mean(a %in% c("T"))
    freq5 <- mean(a %in% c("U"))
    freq6 <- mean(a %in% c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
    if (freq3 > 0.9 && freq4 == 0 && freq5 > 0) {
        return("RNA")
    }
    if (freq6 != 0){
        return("Other")
    }
    return("AA")
}

#guess_skipnum <- function(otufilename){
#    x <- readLines(otufilename, n=50)
#   n <- sum(substr(x, 1, 1) == "#") - 1
#    return(n)
#}

guess_skip_nrow <- function(filename){
    x <- readLines(filename, n=80)
    n <- x %>% 
         strsplit("\t") %>% 
         lapply(length) %>% 
         unlist() %>% 
         magrittr::equals(max(.)) %>% 
         magrittr::not() %>% 
         which()
    if (length(n)==0){
        return (0)
    }
    if (length(n)==1){
        return (n)
    }
    if (length(n)>1 &&       
        all(n %>% 
         diff() %>% 
         magrittr::equals(1))){ 
        n <- max(n)
    }else{
        rlang::abort(paste0("The number ", max(n), " line has different number elements!"))
    }
    return(n)
}

guess_sep <- function(samplefile){
    xx <- readLines(samplefile, n = 8)
    indx <- vapply(strsplit(xx,split=","), function(i)length(i), numeric(1))
    if (all(indx > 1) && var(indx)==0){
        return(",")
    }else{
        return("\t")
    }
}

check_equal_len <- function(x){
   xx <- vapply(x, function(i)length(i), numeric(1))
   return (var(xx) == 0)
}

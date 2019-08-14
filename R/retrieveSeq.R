#' @title Retriveing Sequencing from NCBI
#'
#' @description
#' Retriveing sequences from NCBI with the accession ids.
#' @param ids vector, the accession number or accession.
#' @param files character, the file name specified by a double-quoted string.
#' @param databases character, the name of databases to use, default is `protein`,
#' if nucleotide sequences to retrieve set nuccore,see \code{\link[rentrez]{entrez_fetch}}.
#' @param type character, the format in which to get data,such as fasta, xml ...,
#' see \code{\link[rentrez]{entrez_fetch}}.
#' @param times integer, the time of sleeping, default is 3, 
#' meaning 3 seconds.
#' @param checkids bool, whether check the sequence of ids has been retrieved.
#' default is FALSE.
#' @return the files of sequences downloaded by ids
#' @importFrom Biostrings readBStringSet
#' @importFrom rentrez entrez_fetch
#' @importFrom stringr str_trim
#' @author ShuangbinXu
#' @export
#' @examples
#' retrieveSeq(ids=c("ADM52729.1", "AAF82637.1"), 
#'             files="test.fasta",
#'             databases="protein",
#'             type="fasta",
#'             checkids=TRUE)
retrieveSeq <- function(ids, files, 
						databases="protein", 
						type="fasta", 
						times=3, checkids=FALSE){
	if (file.exists(files) && checkids){
		seqobj <- readBStringSet(files)
		tmpid <- names(seqobj)
		tmpid <- unlist(vapply(strsplit(tmpid, " "),function(x){x[1]},character(1)))
		ids <- setdiff(ids, tmpid)
	}
	if (length(ids)>400){
		stop("The length of ids vector should be shorter than 400!")
	}
	if (length(ids)==0){
		return(NA)
	}
	cat(ids)
	cat("\n")
	tryCatch({tmprecs <- entrez_fetch(db=databases, ids, rettype=type)
		    tmprecs <- gsub("\n\n", "\n", tmprecs)
			tmprecs <- str_trim(tmprecs)
	 	   write(tmprecs, 
			      files,
		          append=TRUE)},
		error=function(e){do.call("retrieveSeq", 
					     args=list(ids=ids,
							 files=files,
							 checkids=TRUE))})
	message(paste0("Sleeping ... ",times,"s"))
	Sys.sleep(times)
	
}

#' @title Retriveing Sequencing from NCBI By mapply
#'
#' @description
#' Retriveing sequences from NCBI with the accession ids.
#' @param idlist vector, the accession version.
#' @param files character, the file name specified by a double-quoted string.
#' @param databases character, the name of databases to use, default is `protein`,
#' if nucleotide sequences to retrieve set nuccore,see \code{\link[rentrez]{entrez_fetch}}.
#' @param type character, the format in which to get data,such as fasta, xml ...,
#' see \code{\link[rentrez]{entrez_fetch}}.
#' @param times integer, the time of sleeping, default is 3, 
#' meaning 3 seconds.
#' @param checkids bool, whether check the sequence of ids has been retrieved.
#' default is FALSE.
#' @return the files of sequences downloaded by ids
#' @seealso \code{\link[retrieveSeq]{MicrobiotaProcess}}
#' @author ShuangbinXu
#' @export
#' @examples
#' idslist <- list(c("ADM52729.1", "AAF82637.1"), c("CAA24729.1", "CAA83510.1"))
#' mapplyretrieveSeq(idlist=idslist,
#'                   files="test.fasta",
#'                   databases="protein",
#'                   type="fasta",
#'                    times=3,checkids=TRUE)
mapplyretrieveSeq <- function(idlist, files, 
							  databases="protein", 
							  type="fasta", 
				              times=3, checkids=TRUE){
	invisible(mapply(retrieveSeq, 
			   idlist, 
	  	 	   MoreArgs=list(files=files,
						  databases=databases,
						  type=type,
				          times=times,
                          checkids=checkids),
			   SIMPLIFY=FALSE))
}

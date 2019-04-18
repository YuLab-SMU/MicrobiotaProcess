#' @title Retriveing Sequencing from NCBI
#'
#' @description
#' Retriveing sequences from NCBI with the accession ids.
#' @param ids vector, the accession number or accession.
#' @param files character, the file name specified by a double-quoted string.
#' @param times integer, the time of sleeping, default is 3, 
#' meaning 3 seconds.
#' @param checkids bool, whether check the sequence of ids has been retrieved.
#' default is FALSE.
#' @importFrom ape read.FASTA read.GenBank
#' @author ShuangbinXu
#' @export

retrieveSeq <- function(ids, files, times=3, checkids=FALSE){
	if (file.exists(files) && checkids){
		seqobj <- read.FASTA(files)
		tmpid <- names(seqobj)
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
	tryCatch({tmpGenBanks <- read.GenBank(ids)
	 	   write.dna(tmpGenBanks, 
			      files,
			      format="fasta", 
			      nbcol=1,
		             append=TRUE, 
			      colw=100)},
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
#' @param ids list, the accession number or accession.
#' @param files character, the file name specified by a double-quoted string.
#' @param times integer, the time of sleeping, default is 3, 
#' meaning 3 seconds.
#' @param checkids bool, whether check the sequence of ids has been retrieved.
#' default is FALSE.
#' @importFrom ape read.FASTA read.GenBank
#' @seealso [retrieveSeq]
#' @author ShuangbinXu
#' @export

mapplyretrieveSeq <- function(idlist, files, 
				  times=3, checkids=TRUE){
	invisible(mapply(retrieveSeq, 
			   idlist, 
	  	 	   MoreArgs=list(files=files,
				          times=times,
                                      checkids=checkids),
			   SIMPLIFY=FALSE))
}


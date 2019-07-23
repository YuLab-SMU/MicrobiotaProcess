#' @author GhuangChuangYu
#' @importFrom grDevices colorRampPalette
#' @keywords internal
# this is from `ggtree`
getCols <- function (n){
     col <- c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072", "#80b1d3",
              "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd",
              "#ccebc5", "#ffed6f")
     col2 <- c("#1f78b4", "#ffff33", "#c2a5cf", "#ff7f00", "#810f7c",
               "#a6cee3", "#006d2c", "#4d4d4d", "#8c510a", "#d73027",
               "#78c679", "#7f0000", "#41b6c4", "#e7298a", "#54278f")
     col3 <- c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99",
               "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a",
               "#ffff99", "#b15928")
         colorRampPalette(col2)(n)
}

#' @keywords internal 
setfactorlevels <- function(data, factorlist){
        factornames <- intersect(colnames(data), names(factorlist))
        if (length(factornames)>0){
                for(i in factornames){
                        data[[i]] <- factor(data[[i]], 
                                   levels=as.vector(factorlist[[i]]))
                }
        }
        return(data)
}

#' @keywords internal
getotudata <- function(obj){
	if(taxa_are_rows(obj)){
		otudata <- data.frame(t(otu_table(obj)), check.names=FALSE)
	}else{
		otudata <- data.frame(otu_table(obj), check.names=FALSE)
	}
}

#' @keywords internal
checkotu <- function(obj){
	if (is.null(obj@otu_table)){
		stop("The otu table is empty!")
	}else{
		otuda <- getotudata(obj)
		return(otuda)
	}
}

#' @keywords internal
checksample <- function(obj){
	if (is.null(obj@sam_data)){
		stop("The sample_data is empty")
	}else{
		sampleda <- getsample(obj)
		return(sampleda)
	}
}

#' @keywords internal.
getsample <- function(obj){
	if (is.null(obj@sam_data)){
		sampleda <- NULL
	}else{
		sampleda <- sample_data(obj)
	}
	return(sampleda)
}

#' @keywords internal
taxlevel <- c("k", "p", "c", "o", "f", "g", "s")

#' @importFrom zoo na.locf
#' @keywords internal
filltaxname <- function(taxdf){
	tmprownames <- rownames(taxdf)
	indexmark <- apply(taxdf, 2, function(x){nchar(x, keepNA = TRUE)})==3
	taxdf[indexmark] <- NA
	indextmp <- apply(is.na(taxdf), 1, which)
	if(length(indextmp)==0){
		return(taxdf)
		break
	}
	taxdf <- apply(taxdf, 1, na.locf)
	taxdf <- lapply(seq_len(ncol(taxdf)), function(i) taxdf[,i])
	newtaxname <- function(x, y){
		y <- as.vector(y)
		x[y] <- paste(taxlevel[y], x[y], sep="__un_")
		x
	}
	taxdf <- data.frame(t(mapply(newtaxname, taxdf, indextmp)), 
						stringsAsFactors=FALSE)
	rownames(taxdf) <- tmprownames
	return(taxdf)
}

#' @keywords internal
addtaxlevel <- function(taxdf){
	#taxlevel <- c("k", "p", "c", "o", "f", "g", "s")
	paste(taxlevel, taxdf, sep="__")
}

#' @keywords internal
fillNAtax <- function(taxdf){
	if (!grepl(taxdf[1,1], "k_")){
		tmprownames <- rownames(taxdf)
		tmpcolnames <- colnames(taxdf)
		taxdf <- t(apply(taxdf, 1, as.character))
		taxdf[is.na(taxdf)] <- ""
		taxdf <- data.frame(t(apply(taxdf, 1, addtaxlevel)),
							stringsAsFactors=FALSE)
		rownames(taxdf) <- tmprownames
		colnames(taxdf) <- tmpcolnames
	}
	taxdf <- filltaxname(taxdf)
	return(taxdf)
}


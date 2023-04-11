#' @author GhuangChuangYu
#' @importFrom grDevices colorRampPalette
#' @keywords internal
# this is from `ggtree`
get_cols <- function (n){
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
                    data[[match(i,colnames(data))]] <- factor(data[[match(i, colnames(data))]], 
                               levels=as.vector(factorlist[[match(i,names(factorlist))]]))
            }
    }
    return(data)
}

#' @keywords internal
get_otudata <- function(obj){
    otudata <- obj@otu_table
    otudata <- data.frame(otudata, check.names=FALSE)
    if(obj@otu_table@taxa_are_rows){
        otudata <- data.frame(t(otudata), check.names=FALSE)
    }
    return (otudata)
}

#' @keywords internal
checkotu <- function(obj){
    if (is.null(obj@otu_table)){
    	stop("The otu table is empty!")
    }else{
    	otuda <- get_otudata(obj)
    	return(otuda)
    }
}

#' @keywords internal
checksample <- function(obj){
    if (is.null(obj@sam_data)){
    	stop("The sample_data is empty")
    }else{
    	sampleda <- get_sample(obj)
    	return(sampleda)
    }
}

#' @keywords internal.
get_sample <- function(obj){
    if (is.null(obj@sam_data)){
    	sampleda <- NULL
    }else{
    	sampleda <- data.frame(obj@sam_data, check.names=FALSE)
    }
    return(sampleda)
}

#' @keywords internal
taxlevelchar <- c("k", "p", "c", "o", "f", "g", "s", "st")

newtaxname <- function(x, y){
    y <- as.vector(y)
    x[y] <- paste(taxlevelchar[y], x[y], sep="__un_")
    x
}

#' @importFrom zoo na.locf
#' @keywords internal
filltaxname <- function(taxdf, type="species"){
    tmprownames <- rownames(taxdf)
    indexmark <- apply(taxdf, 2, function(x){nchar(x, keepNA = TRUE)})<=4
    taxdf[indexmark] <- NA
    if (any(is.na(taxdf[,1]))){
        if (type == "species"){
            prefix <- "k__"
        }else{
            prefix <- "d1__"
        }
        taxdf[is.na(taxdf[,1]), 1] <- paste0(prefix, "Unknown")
    }    
    indextmp <- apply(is.na(taxdf), 1, which)
    if(length(indextmp)==0){
        taxdf <- data.frame(taxdf, check.names=FALSE)
        return(taxdf)
    }
    taxdf <- apply(taxdf, 1, na.locf)
    taxdf <- lapply(seq_len(ncol(taxdf)), function(i) taxdf[,i])
    taxdf <- data.frame(t(mapply(newtaxname, taxdf, indextmp)), 
                        stringsAsFactors=FALSE)
    rownames(taxdf) <- tmprownames
    return(taxdf)
}

#' @keywords internal
addtaxlevel <- function(taxdf, taxlevelchar){#, type="species"){
    #if (type != "species"){
    #    taxlevelchar <- paste0("d", seq_len(ncol(taxdf)))
    #}else{
    #    taxlevelchar <- taxlevelchar[seq_len(ncol(taxdf))]
    #}
    taxlevelchar <- taxlevelchar[seq_len(length(taxdf))]
    paste(taxlevelchar, taxdf, sep="__")
}

#' @importFrom tibble column_to_rownames
#' @keywords internal
fillNAtax <- function(taxdf, type="species"){
    #taxdf <- remove_unclassfied(taxdf)
    taxdf <- remove_na_taxonomy_rank(taxdf)
    if (type!="species" || ncol(taxdf) >= 9){
        #assign("taxlevelchar", paste0("d", seq_len(ncol(taxdf))), envir = .GlobalEnv)
        taxlevelchar <- paste0('d', seq_len(ncol(taxdf)))
    }else{
        #assign("taxlevelchar", c("k", "p", "c", "o", "f", "g", "s", "st"), envir = .GlobalEnv)
        taxlevelchar <- c("k", "p", "c", "o", "f", "g", "s", "st")
    }
#    if (any(is.na(taxdf[,1]))){
#        if (type == "species"){
#            prefix <- "k__"
#        }else{
#            prefix <- "d1__"
#        }
#        taxdf[is.na(taxdf[,1]), 1] <- paste0(prefix, "Unknown")
#    }
    if (!(grepl("^k__", taxdf[1, 1], ignore.case = TRUE) || 
          grepl("^d1__", taxdf[1, 1], ignore.case = TRUE) || 
          grepl("^d__", taxdf[1, 1], ignore.case = TRUE)
           )
    ){
    	tmprownames <- rownames(taxdf)
    	tmpcolnames <- colnames(taxdf)
        taxdf <- t(apply(taxdf, 1, as.character))
        taxdf[is.na(taxdf)] <- ""
        taxdf <- data.frame(t(apply(taxdf, 1, addtaxlevel, taxlevelchar)),
                            stringsAsFactors=FALSE)
    	rownames(taxdf) <- tmprownames
    	colnames(taxdf) <- tmpcolnames
    }
    taxdf <- filltaxname(taxdf, type = type)
    taxdf <- repduplicatedtaxcheck(taxdf) #%>% column_to_rownames(var="rowname")
    attr(taxdf, "fillNAtax") <- TRUE 
    return(taxdf)
}

#' @keywords internal
remove_unclassfied <- function(taxdf){
    taxdf[grepl.data.frame("Unclassified|uncultured|Ambiguous|Unknown|unknown|metagenome|Unassig", taxdf, ignore.case=TRUE)] <- NA
    return(taxdf)
}

grepl.data.frame <- function(pattern, x, ...){
    y <- if (length(x)) {
             do.call("cbind", lapply(x, "grepl", pattern=pattern, ...))
         }else{
             matrix(FALSE, length(row.names(x)), 0)
         }
    if (.row_names_info(x) > 0L)
        rownames(y) <- row.names(x)
    y
}


#' @importFrom magrittr %>%
#' @importFrom tibble rownames_to_column column_to_rownames
#' @keywords internal
duplicatedtaxcheck <- function(taxdf){
    if (ncol(taxdf)==1){return(taxdf)}
    taxdf <- taxdf %>% rownames_to_column()
    for (i in ncol(taxdf):3){
    	tmp <- split(taxdf,taxdf[,i])
    	for (j in seq_len(length(tmp))){
    		flag <- length(unique(as.vector(tmp[[j]][,i-1])))
    		if (flag > 1){
    			tmp[[j]][,i] <- paste(tmp[[j]][,i],tmp[[j]][,i-1],sep="_")
    		}
    	}
    	taxdf <- do.call("rbind",c(tmp, make.row.names=FALSE)) 
    }
    return(taxdf)
}

#' @keywords internal
repduplicatedtaxcheck <- function(taxdf){
    for (i in seq_len(7)){
    	taxdf <- duplicatedtaxcheck(taxdf) %>% 
		column_to_rownames(var="rowname")
    }
    return(taxdf)
}

## #' @keywords internal
## ## reference https://rdrr.io/cran/stackoverflow/man/match.call.defaults.html
## match.call.defaults <- function(fun) {
##     if (!is.na(fun)){
##         print(args(diff_analysis.data.frame))
##         args(diff_analysis.data.frame)
##     }else{
##         call <- evalq(match.call(expand.dots=TRUE), parent.frame(1))
##         formals <- evalq(formals(), parent.frame(1))
##         for(i in setdiff(names(formals), c(names(call)))){
##             call[i] <- list(formals[[i]])
##         }
##         match.call(sys.function(sys.parent()), call)
##     }
## }

#' @keywords internal
extract_args <- function(obj, arg){
    if (!"someparams" %in% methods::slotNames(obj)){
        stop("The object don't have someparams slot!")
    }else{
        args <- obj@someparams
        argres <- args[[arg]]
        return(argres)
    }
}

# #' @importFrom utils globalVariables
# utils::globalVariables('taxlevelchar')

#' @importFrom stats sd 
#' @importFrom plyr ddply
# Adapted from Rmisc
summarySE <- function (data = NULL, measurevar, groupvars = NULL, na.rm = FALSE,
                       conf.interval = 0.95, .drop = TRUE){
    length2 <- function(x, na.rm = FALSE) {
        if (na.rm)
            sum(!is.na(x))
        else length(x)
    }
    datac <- ddply(data, 
                   groupvars, 
                   .drop = .drop, 
                   .fun = function(xx, col, na.rm) {
          c(N = length2(xx[, col], na.rm = na.rm), mean = mean(xx[, col], na.rm = na.rm), 
            sd = sd(xx[, col], na.rm = na.rm))
          }, measurevar, na.rm)
    datac %<>% dplyr::rename(!!measurevar:="mean")
    datac$se <- datac$sd/sqrt(datac$N)
    ciMult <- qt(conf.interval/2 + 0.5, datac$N - 1)
    datac$ci <- datac$se * ciMult
    return(datac)
}

#' @importFrom stats qt
# reference to Rmisc
CI <- function (x, ci = 0.95, na.rm=FALSE){
    a <- mean(x, na.rm = na.rm)
    s <- sd(x, na.rm = na.rm)
    if (na.rm){
        x <- x[!is.na(x)]
    }
    n <- length(x)
    error <- qt(ci + (1 - ci)/2, df = n - 1) * s/sqrt(n)
    return(c(upper = a + error, mean = a, lower = a - error))
}

.return_wrap <- function(...){
    msg <- paste(..., collapse = "", sep = "")
    wrapped <- strwrap(msg, width = getOption("width") - 2) %>%
               glue::glue_collapse(., "\n", last = "\n")
    wrapped
}

message_wrap <- function(...){
    msg <- .return_wrap(...)
    message(msg)
}

stop_wrap <- function(...){
    msg <- .return_wrap(...)
    stop(msg, call. = FALSE)
}

warning_wrap <- function(...){
    msg <- .return_wrap(...)
    warning(msg, call. = FALSE)
}

remove_MP_internal_res <- function(x){
    x <- x[,!vapply(x, function(i)is.list(i), logical(1))]
    index <- lapply(MP_internal_res, function(i)which(grepl(paste0("^",i), colnames(x)))) %>%
             unlist() %>%
             unique()
    if (length(index) > 0) x <- x[-index]
    return(x)
}

MP_internal_res <- c("Observe", "Chao1", "ACE", "Shannon", 
                     "Simpson", "Pielou", 'PC', 'PCo', 
                     'CA', 'NMDS', 'RDA', 'CCA', "PAE", "NRI", 
                     "NTI", "PD", 'IAC', "HAED", "EAED")

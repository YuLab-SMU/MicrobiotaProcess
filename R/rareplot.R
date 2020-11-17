#' @method ggrarecurve data.frame
#' @rdname ggrarecurve
#' @export
ggrarecurve.data.frame <- function(obj, sampleda, factorLevels, chunks=400, ...){
    obj <- get_rarecurve(obj=obj, 
                         sampleda=sampleda, 
                         chunks=chunks, 
                         factorLevels=factorLevels)
    p <- ggrarecurve.rarecurve(obj=obj, ...)
    return(p)
}

#' @method ggrarecurve rarecurve
#' @importFrom ggplot2 ggplot geom_ribbon aes_string geom_smooth facet_wrap scale_y_continuous
#' @importFrom dplyr filter
#' @importFrom rlang .data
#' @importFrom scales squish
#' @importFrom Rmisc summarySE
#' @rdname ggrarecurve
#' @export
ggrarecurve.rarecurve <- function(obj,
    #sampleda, 
    indexNames="Observe", 
    linesize=0.5, 
    facetnrow=1,
    shadow=TRUE, 
    #chunks=400, 
    factorNames, 
    #factorLevels, 
    se=FALSE,
    method="lm", 
    formula=y ~ log(x), ...){
    obj <- obj$data 
    #obj <- get_rarecurve(obj=obj, sampleda=sampleda, chunks=chunks, factorLevels=factorLevels)
    #obj <- stat_rare(data=obj, chunks=chunks, sampleda=sampleda, factorLevels=factorLevels, plotda=TRUE)
    mapping <- aes_string(x="readsNums", y="value", color="sample")
    if (!missing(factorNames)){
        if (shadow){
            obj <- summarySE(obj, measurevar="value", groupvars=c(factorNames, "readsNums", "Alpha"), na.rm=TRUE)
            obj$up <- obj$value - obj$sd
            obj$down <- obj$value + obj$sd
            mapping <- modifyList(mapping, aes_string(group=factorNames, color=factorNames, fill=factorNames, ymin="up", ymax="down"))
        }else{
            mapping <- modifyList(mapping, aes_string(group="sample", color=factorNames))
        }
    }
    if (!is.null(indexNames)){
        obj <- obj %>% filter(.data$Alpha %in% indexNames)
    }
    p <- ggplot(data=obj, mapping=mapping) #+
    if (!missing(factorNames) && shadow){
        #p <- p + geom_errorbar(alpha=0.5)
        p <- p + geom_ribbon(alpha=0.3, color=NA, show.legend=FALSE)
    }    
    message("The color has been set automatically, you can reset it manually by adding scale_color_manual(values=yourcolors)")
    p <- p + geom_smooth(se=se, method = method, size=linesize,formula = formula,...)+
         scale_y_continuous(limits=c(0,NA), oob=squish) +
         facet_wrap(~ Alpha, scales="free", nrow=facetnrow) +
         ylab("alpha metric")+xlab("number of reads")
    return(p)
}

#' @title mapping data of ggrarecurve
#' @description
#' generating the data of ggrarecurve.
#' @param data data.frame,(nrow sample * ncol taxonomy 
#' (feature) or and factor)
#' @param chunks integer, the number of subsample in a sample,
#' default is 400.
#' @param sampleda data.frame, (nrow sample * ncol factor)
#' @param factorLevels list, the levels of the factors, default is NULL,
#' if you want to order the levels of factor, you can set this.
#' @param plotda boolean, default is TRUE, whether build the data.frame for
#' `geom_bar` of `ggplot2`.
#' @return data.frame for ggrarecurve.
#' @author Shuangbin Xu
#' @importFrom dplyr bind_rows
#' @importFrom reshape melt
#' @importFrom magrittr %>%
#' @keywords internal
stat_rare <- function(data, 
    chunks=400, 
    sampleda,
    factorLevels,
    plotda=TRUE){
    tmpfeature <- colnames(data)[vapply(data,is.numeric,logical(1))]
    tmpfactor <- colnames(data)[!vapply(data,is.numeric,logical(1))]
    dat <- data[ , match(tmpfeature, colnames(data)), drop=FALSE]
    out <- apply(dat, 1, samplealpha, chunks=chunks) %>% 
    	bind_rows(,.id="sample")
    if (plotda){
    	if (!missing(sampleda)){
    		sampleda$sample <- rownames(sampleda)
    		out <- merge(out, sampleda)
    		out <- melt(out,id.vars=c(colnames(sampleda), "readsNums"),
    					variable_name="Alpha")
    	}
    	if (missing(sampleda) && length(tmpfactor) > 0){
    		tmpsample <- data[, tmpfactor, drop=FALSE]
    		tmpsample$sample <- rownames(tmpsample)
    		out <- merge(out, tmpsample)
    		out <- melt(out, id.vars=c("sample", "readsNums", tmpfactor),
    					variable_name="Alpha")
    	}
    	if (missing(sampleda)&&length(tmpfactor) == 0){
    		out <- melt(out, id.vars=c("sample", "readsNums"),
    			     variable_name="Alpha")
    	}		
    }else{
    	if (!missing(sampleda)){
    		sampleda$sample <- rownames(sampleda)
    		out <- merge(out, sampleda)
    	}
    	if (missing(sampleda) && length(tmpfactor) >0){
    		tmpsample <- data[,tmpfactor,drop=FALSE]
    		tmpsample$sample <- rownames(tmpsample)
    		out <- merge(out, tmpsample)
    	}
    }
    if (!missing(factorLevels)){
    	out <- setfactorlevels(out, factorLevels)
    }
    return(out)
}

#' @keywords internal
samplealpha <- function(data, chunks=200){
    sdepth <- sum(data)
    step <- trunc(sdepth/chunks)
    n <- seq(0, sdepth, by=step)[-1]
    n <- c(n, sdepth)
    out <- lapply(n, function(x){
                tmp <- get_alphaindex(data, mindepth=x)
    		#tmp <- tmp$indexs
    		tmp$readsNums <- x
    	    return(tmp)})
    out <- do.call("rbind", c(out, make.row.names=FALSE))
    out[is.na(out)] <- 0
    return (out)
}


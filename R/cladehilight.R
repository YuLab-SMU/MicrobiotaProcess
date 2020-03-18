#' @author Shuangbin Xu, GuangChuang Yu
#' @keywords internal 
getcladedf <- function(ggtree, node){
    data <- ggtree$data
    if ("nodeClass" %in% colnames(data)){
        data <- set_newlevels(data=data, 
                newlevels=taxlevel[seq_len(length(unique(data$nodeClass)))],
                factorNames="nodeClass")
        levelsnum <- length(levels(data$nodeClass)) + 1
        tmpnum <- levelsnum - as.numeric(data$nodeClass)
        data$extend <- get_extend(tmpnum)
    }
    df <-lapply(node, function(x)get_clade_position_(data=data, node=x))
    df <- do.call("rbind", df)
    df$node <- node
    df <- merge(df, data, by.x="node", by.y="node")
    df$xmax <- df$xmax + df$extend
    return(df)
}

#' @author Shuangbin Xu, GuangChuang Yu
#' @keywords internal
getlabeldf <- function(ggtree, node, angle="auto"){
    data <- ggtree$data
    if ("nodeClass" %in% colnames(data)){
        data <- set_newlevels(data=data,
                newlevels=taxlevel[seq_len(length(unique(data$nodeClass)))],
                factorNames="nodeClass")
    	levelsnum <- length(levels(data$nodeClass)) + 1
    	tmpnum <- levelsnum - as.numeric(data$nodeClass)
    	data$levelindex <- tmpnum
    	data$extend <- get_extend(tmpnum)*0.9
    }
    df <- lapply(node, function(x)getcladelabelposition(data, 
                 x, angle=angle, angleoff=90))
    df <- do.call("rbind", df)
    df$node <- node
    df$y <- as.numeric(apply(df[,match(c("y", "yend"), colnames(df))],
                             1, mean)) 
    data <- data %>% dplyr::select(-c("x", "y", "angle"))
    df <- merge(df, data, by.x="node", by.y="node")
    df$x <- df$x + df$extend
    return(df)
}

#' @author GuangChuang Yu, Shuangbin Xu
#' @keywords internal
getcladelabelposition <- function(data, 
                                  node, angle = "auto", 
                                  adjustRatio=1, angleoff=90, 
                                  extend = 0) { 
    if (length(extend) == 1) {
    	 extend <- rep(extend, 2)
    }
    sp <- tryCatch(offspringtbl_tree(data, node)$node, error=function(e) NULL)
    i <- match(node, data$node)
    if (is.null(sp)) {
    	sp.df <- data[i,]
    }else{
    	sp <- c(sp, node)
    	sp.df <- data[match(sp, data$node),]
    }
    y <- sp.df$y
    y <- y[!is.na(y)]
    mx <- max(sp.df$x, na.rm=TRUE)
    d <- data.frame(x=mx, y=min(y) - extend[2], yend=max(y) + extend[1])
    if (missing(angle)) 
    	return(d)
    if (angle == "auto") {
    	angletmp <- mean(range(sp.df$angle))
    	if(angletmp >180){
    		d$angle <- angletmp + 90
    	}else{
    		d$angle <- angletmp + 270
    	}
    }else{
    	d$angle <- angle
    }
    mx <- d$x
    mx <- mx * adjustRatio
    data.frame(x=mx, xend=mx, y=d$y, yend=d$yend, angle=d$angle)
}

#' @author Shuangbin Xu, GuangChuang Yu
#' @keywords internal
getannotlabel <- function(labeldf,classlevel=4){
    df <- labeldf[labeldf$levelindex <= classlevel, ]
    dat <- labeldf[labeldf$levelindex > classlevel, ]
    lett <- c(letters, toupper(letters))
    if (nrow(df)>52){
        lengthtmp <- round(nrow(df)/52,0) + 1
        lefttmp <- nrow(df) - 52
        annolabel <- c(lett,paste0(rep(lett,lengthtmp)[seq_len(lefttmp)], seq(1,lefttmp))) 
    }else{
        annolabel <- lett[seq_len(nrow(df))]
    }
    tmplabels <- paste(annolabel, df$label,sep=": ")
    df$label <- annolabel
    dat <- rbind(df, dat)
    df$label <- tmplabels
    return(list(labeldf=dat, annotdf=df))
}

###' @author GuangChuang Yu
###' @keywords internal
#get_clade_position_ <- function(data, node){
#    sp <- tryCatch(offspringtbl_tree(data, node)$node, error=function(e) NULL)
#    i <- match(node, data$node)
#    if (is.null(sp)) {
#    	sp.df <- data[i,]
#    }else{
#    	sp <- c(sp, node)
#    	sp.df <- data[match(sp, data$node),] 
#    }
#    x <- sp.df$x
#    y <- sp.df$y
#    if ("branch.length" %in% colnames(data)) {
#    	xmin <- min(x, na.rm=TRUE)-data[["branch.length"]][i]/2
#    }else {
#    	xmin <- min(sp.df$branch, na.rm=TRUE)
#    }
#    data.frame(xmin=xmin,
#    		   xmax=max(x, na.rm=TRUE),
#    		   ymin=min(y, na.rm=TRUE) - 0.5, 
#    		   ymax=max(y, na.rm=TRUE) + 0.5)
#}


#' @keywords internal
get_extend <- function(x) {x * 0.5}

#' @importFrom utils getFromNamespace
#' @keywords internal
offspringtbl_tree <- getFromNamespace("offspring.tbl_tree","tidytree")

#' @importFrom utils getFromNamespace
#' @keywords internal
get_clade_position_ <- getFromNamespace("get_clade_position_", "ggtree")

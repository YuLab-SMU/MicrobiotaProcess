#' @author Shuangbin Xu, GuangChuang Yu
#' @keywords internal 
getcladedf <- function(ggtree, node){
	data <- ggtree$data 
	if ("nodeClass" %in% colnames(data)){
		levelsnum <- length(levels(data$nodeClass)) + 1
		tmpnum <- levelsnum - as.numeric(data$nodeClass)
		data$extend <- get_extend(tmpnum)
	}
	df <-lapply(node, function(x)ggtree:::get_clade_position_(data, x))
	df <- do.call("rbind", df)
	df$node <- node
	df <- merge(df, data, by.x="node", by.y="node")
	df$xmax <- df$xmax + df$extend
	return(df)
}

#' @author Shuangbin Xu, GuangChuang Yu
#' @keywords internal
getlabeldf <- function(ggtree, node){
	data <- ggtree$data 
	if ("nodeClass" %in% colnames(data)){
		levelsnum <- length(levels(data$nodeClass)) + 1
		tmpnum <- levelsnum - as.numeric(data$nodeClass)
		data$levelindex <- tmpnum
		data$extend <- get_extend(tmpnum)*0.9
	}
	df <- lapply(node, function(x)getcladelabelposition(data, x, angle="auto", angleoff=90))
	df <- do.call("rbind", df)
	df$node <- node
	df$y <- as.numeric(apply(df[,match(c("y", "yend"), colnames(df))], 1, mean)) 
	data <- data %>% dplyr::select(-c("x", "y", "angle"))
	df <- merge(df, data, by.x="node", by.y="node")
	df$x <- df$x + df$extend
	return(df)
}

#' @author GuangChuang Yu, Shuangbin Xu
#' @keywords internal
getcladelabelposition <- function(data, node, angle = "auto", adjustRatio=1, angleoff=90, extend = 0) { 
	if (length(extend) == 1) {
		 extend = rep(extend, 2)
	}
	sp <- tryCatch(tidytree:::offspring.tbl_tree(data, node)$node, error=function(e) NULL)
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
		d$angle <- mean(range(sp.df$angle)) + angleoff
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
	if (nrow(df)>26){
		lengthtmp <- round(nrow(df)/26,0) + 1
		lefttmp <- nrow(df) - 26
		annolabel <- c(letters,paste0(rep(letters,lengthtmp)[1:lefttmp], seq(1,lefttmp))) 
	}else{
		annolabel <- letters[1:nrow(df)]
	}
	tmplabels <- paste(annolabel, df$label,sep=": ")
	df$label <- annolabel
	dat <- rbind(df, dat)
	df$label <- tmplabels
	return(list(labeldf=dat, annotdf=df))
}

#' @keywords internal
get_extend <- function(x) {x * 0.5}


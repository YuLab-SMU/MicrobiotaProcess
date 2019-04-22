#' @author GhuangChuangYu
#' @importFrom grDevices colorRampPalette
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

#' @keyword internal 
setfactorlevels <- function(data, factorlist){
        factornames <- intersect(colnames(dat), names(factorlist))
        if (length(factornames)>0){
                for(i in factornames){
                        data[[i]] <- factor(data[[i]], 
                                   levels=as.vector(factorlist[[i]]))
                }
        }
        return(data)
}


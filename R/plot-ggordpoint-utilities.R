#' @title add confidence ellipse to ordinary plot
#' @param data, default is NULL.
#' @param mapping, aes mapping, default is NULL.
#' @param ellipse_pro numeric, confidence value for the ellipse, default is 0.9
#' @param alpha numeric, alpha of ellipse, default is 0.2
#' @param show.legend, default is NA.
#' @param inherit.aes, default is TRUE
#' @param ... additional parameters
#' @return ggplot layer
#' @importFrom ggplot2 aes_
#' @importFrom ggplot2 layer
#' @importFrom utils modifyList
#' @author Guangchuang Yu
#' @references \url{http://lchblogs.netlify.com/post/2017-12-22-r-addconfellipselda/}
#' @keywords internal
#' @noRd
geom_ord_ellipse <- function(data=NULL, mapping = NULL, ellipse_pro = 0.9, alpha=0.3, show.legend=NA, inherit.aes = TRUE, ...) {
    #default_aes <- aes_(color = ~Groups, group = ~Groups)
    #if (is.null(mapping)) {
    #    mapping <- default_aes
    #} else {
    #    mapping <- modifyList(default_aes, mapping)
    #}
    layer(
        geom = "polygon",
        stat = StatOrdEllipse,
        mapping = mapping,
        position = 'identity',
        show.legend = show.legend,
        data = data,
        inherit.aes = inherit.aes,
        params = list(
            ellipse_pro = ellipse_pro,
            alpha = alpha,
            ...
        )
    )
}

#' @importFrom ggplot2 ggproto
#' @importFrom ggplot2 Stat
#' @importFrom grDevices chull
#' @keywords internal
StatOrdEllipse <- ggproto("StatOrdEllipse", Stat,
                          compute_group = function(self, data, scales, params, ellipse_pro) {
                              #names(data)[seq_len(2)] <- c('one', 'two')
                              theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
                              circle <- cbind(cos(theta), sin(theta))
                              ell <- ddply(data, .(group), function(x) {
                                  if(nrow(x) <= 2) {
                                      return(NULL)
                                  }
                                  sigma <- var(cbind(x$x, x$y))
                                  mu <- c(mean(x$x), mean(x$y))
                                  ed <- sqrt(qchisq(ellipse_pro, df = 2))
                                  data.frame(sweep(circle %*% chol(sigma) * ed, 2, mu, FUN = '+'))
                              })
                              names(ell)[2:3] <- c('x', 'y')
                              ell <- ddply(ell, .(group), function(x) x[chull(x$x, x$y), ])
                              names(ell) <- c('Groups', 'x', 'y')
                              return(ell)
                          },
                          required_aes = c("x", "y", "group")
                          )

#' @keywords internal
## . function was from plyr package
. <- function (..., .env = parent.frame()) {
    structure(as.list(match.call()[-1]), env = .env, class = "quoted")
}

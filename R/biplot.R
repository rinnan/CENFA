#' biplot
#'
#' Biplots of \code{cnfa} and \code{enfa} objects.
#'
#' @param x an object of class \code{cnfa} or \code{enfa} describing the
#'   occupied habitat
#' @param global an object of class \code{GLcnfa} describing the global reference
#'   habitat
#' @param xax the column number of the x-axis
#' @param yax the column number of the y-axis
#' @param p the proportion of observations to include in the calculations of
#'   the minimum convex polygons
#' @param n the number of projected variables to label
#' @param plot if \code{TRUE}, plot will be returned on function call
#' @param ... additional \code{plot} arguments
#'
#' @examples
#' mod1 <- cnfa(x = climdat.hist, s.dat = ABPR, field = "CODE")
#' glc <- GLcenfa(x = climdat.hist)
#' biplot(x = mod1, global = glc)
#'
#' @importFrom grDevices chull
#' @importFrom magrittr %>%
#' @importFrom stats ecdf sd
#' @importFrom graphics abline arrows legend par points polygon
#'
#' @export

# @rdname biplot
#if (!isGeneric("biplot")) {
  setGeneric("biplot", function(x, global, ...){
    standardGeneric("biplot")})
#}

#' @rdname biplot
setMethod("biplot",
          signature(x = "cnfa", global = "GLcenfa"),
          function(x, global, xax = 1, yax = 2, p = 0.99, n = 5, plot = TRUE, ...){

            s.ras <- raster(x)
            g.ras <- raster(global)

            good <- c(xax, yax) %in% 1:nlayers(s.ras)

            if(FALSE %in% good) {
              xax <- 1
              yax <- 2
              warning("selected axes not meaningful, using first two axes instead.")
            }

            co <- x@co[, c(xax, yax)]

            #if(canProcessInMemory(g.ras)) {
            g.dat <- na.omit(values(g.ras))
            g <- g.dat %*% co
            s <- na.omit(values(s.ras)[, c(xax, yax)])

            if(p < 1 & p > 0) {
              g.centroid  <- colMeans(g)
              gdists <- sweep(g, 2, g.centroid, "-")
              r1 <- ecdf(gdists[ ,1])
              r2 <- ecdf(gdists[, 2])
              outs <- sqrt(r1(gdists[ ,1])^2 + r2(gdists[ ,2])^2)
              g.qn <- which(outs < quantile(outs, p))
              g.in <- g[g.qn, ]
              g.ch <- chull(g.in)

              s.centroid  <- colMeans(s)
              sdists <- sweep(s, 2, s.centroid, "-")
              r1 <- ecdf(sdists[ ,1])
              r2 <- ecdf(sdists[, 2])
              outs <- sqrt(r1(sdists[ ,1])^2 + r2(sdists[ ,2])^2)
              s.qn <- which(outs < quantile(outs, p))
              s.in <- s[s.qn, ]
              s.ch <- chull(s.in)
            } else if(p == 1) {
              g.in <- g
              s.in <- s
              g.ch <- chull(g.in)
              s.ch <- chull(s.in)
            } else stop("percentage must be in the range [0, 100]")

            #} else {

            #   f1 <- function(p) as.numeric(p %*% co)
            #   g <- calc(g.ras, f1, forceapply = T)
            #   s <- subset(s.ras, c(xax, yax))
            #
            #   if(p < 1 & p > 0) {
            #
            #     gdists <- scale(g, center = T, scale = F)
            #     f2 <- function(p) sqrt(sum(p^2))
            #     gdists <- calc(gdists, f2, forceapply = T)
            #     gqn <- quantile(gdists, probs = p)
            #     g.qn <- which(values(gdists) < gqn) %>% values(g)[., ]
            #     g.ch <- chull(g.qn)
            #
            #     gqn <- quantile(gdists, probs = p)
            #     g.qn <- which(values(gdists) < gqn) %>% values(g)[., ]
            #     g.ch <- chull(g.qn)
            #
            #     sdists <- scale(s, center = T, scale = F)
            #     sdists <- calc(sdists, f2, forceapply = T)
            #     sqn <- quantile(sdists, probs = p)
            #     s.qn <- which(values(sdists) < sqn) %>% values(s)[., ]
            #     s.ch <- chull(s.qn)
            #
            #   } else if(p == 1) {
            #     g.qn <- na.omit(values(g))
            #     g.ch <- chull(g.qn)
            #     s.qn <- na.omit(values(s))
            #     s.ch <- chull(s.qn)
            #   } else error("percentage must be in the range [0, 100]")
            # }


              xmin <- min(g.in[ ,1], s.in[ ,1])
              xmax <- max(g.in[ ,1], s.in[ ,1])
              xfact <- xmax - xmin
              xmin <- xmin - xfact*.1
              xmax <- xmax + xfact*.1

              ymin <- min(g.in[ ,2], s.in[ ,2])
              ymax <- max(g.in[ ,2], s.in[ ,2])
              yfact <- ymax - ymin
              ymin <- ymin - yfact*.1
              ymax <- ymax + yfact*.1

            mags <- apply(co, 1, norm, "2") %>% order(decreasing = T) %>% .[1:n]

            par(mar = c(1, 1, 1, 1))
            plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax),
                 xlab = NA, ylab = NA, axes = F, ann = F, ...)
            abline(v = g.centroid[1], h = g.centroid[2])#, col = "grey70")
            polygon(g.in[g.ch, ])
            polygon(s.in[s.ch, ], col = "grey60", xpd = T)
            points(s.centroid[1], s.centroid[2], pch = 21, bg = "white")

            .adjust_arrows(x = co[, 1], y = co[, 2], xfact, yfact, xpd=T, length = .05)

            # x <- rep(xfact/40, length(co[ ,1]))
            # x[co[, 1] < 0] <- -x[co[, 1] < 0]
            # y <- rep(yfact/40, length(co[ ,2]))
            # y[co[, 2] < 0] <- -y[co[, 2] < 0]

            #text(co[mags, 1] * xfact + x[mags], co[mags, 2] * yfact + y[mags], labels = rownames(co[mags, ]), cex = .7, xpd = T)
            .adjust_labels(co[mags, 1], co[mags, 2], xfact, yfact, labels = rownames(co[mags, ]), cex = .7, xpd = T)
            legend("topright",
                   legend = c(as.expression(bquote(xax == .(names(s.ras)[xax]))),
                              as.expression(bquote(yax == .(names(s.ras)[yax]))),
                              as.expression(bquote(~~~p == .(p)))),
                   bty = "n")


            # beginCluster(n, exclude = "CENFA")
            #
            # dat <- clusterR(y, calc, args = list(fun=f), export = "co", progress = "text", style = 3, m = 4)
            # names(dat) <- colnames(co)
            #
            # endCluster()

          }
)

#' @keywords internal
.adjust_labels <- function(x, y, xfact, yfact, ...) {
  #xfact <- xmax - xmin
  #yfact <- ymax - ymin

  p <- rep(xfact/40, length(x))
  p[x < 0] <- -p[x < 0]
  q <- rep(yfact/40, length(y))
  q[y < 0] <- -q[y < 0]

  text(x * xfact + p, y * yfact + q, ...)

}

#' @keywords internal
.adjust_arrows <- function(x, y, xfact, yfact, ...) {
  #xfact <- xmax - xmin
  #yfact <- ymax - ymin

  arrows(x0 = 0, y0 = 0, x1 = x * xfact, y1 = y * yfact, ...)
}

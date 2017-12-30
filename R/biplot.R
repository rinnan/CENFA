#' @importFrom grDevices chull
#' @importFrom magrittr %>%
#' @importFrom graphics abline arrows legend points polygon
#' @importFrom stats sd

biplot <- function(x, global, xax = 1, yax = 2, percentage = .99, cores = 1){

  y <- raster(global)
  co <- as.matrix(x@co)[, c(xax, yax)]
  co[,1] <- co[,1] / sqrt(t(co[,1]) %*% co[,1])
  xmin <- min(gli.qn[,1])
  ymin <- min(co[,2])*1.1
  xmax <- max(gli.qn[,1])
  ymax <- max(co[,2])*1.1
  f <- function(p) as.numeric(p %*% co)
  #filename
  beginCluster(cores, exclude = "CENFA")

  dat <- clusterR(y, calc, args = list(fun=f), export = "co", progress = "text", style = 3, m = 4)
  names(dat) <- colnames(co)

  endCluster()

  gpres <- which(!is.na(values(max(dat))))
  gli <- values(dat)[gpres, ]
  g.sds <- apply(gli, 2, sd, na.rm = T)
  gli.c <- scale(gli, center = F, scale = g.sds)
  gcentroid  <- colMeans(gli, na.rm = T)
  gdists <- sweep(gli.c, 2, gcentroid, "-")
  gdists <- sqrt(rowSums(gdists^2))
  gqn <- quantile(gdists, probs = percentage)
  gli.qn <- gli[which(gdists < gqn),]
  gch <- chull(gli.qn)

  y <- raster(x)[[c(xax, yax)]]
  pres <- which(!is.na(values(max(y))))
  li <- values(y)[pres, ]
  li[,1] <- li[,1]/sqrt(t(x@mf) %*% x@mf)
  sds <- apply(li, 2, sd, na.rm = T)
  li.c <- scale(li, center = F, scale = sds)
  centroid  <- colMeans(li, na.rm = T)
  dists <- sweep(li.c, 2, centroid, "-")
  dists <- sqrt(rowSums(dists^2))
  qn <- quantile(dists, probs = percentage)
  li.qn <- li[which(dists < qn),]
  ch.qs <- chull(li.qn)

  mags <- apply(co, 1, norm, "2") %>% order(decreasing = T) %>% .[1:5]
  fact <- (ymax - ymin)/(xmax - xmin)
  xmin <- min(xmin, max(co[,1])/fact*1.05)

  plot(NA, xlim = c(xmin, xmax), ylim = c(ymin, ymax), xlab = names(dat)[1], ylab = names(dat)[2], axes = F, ann = F)
  #axis(1, at = c(-25,0,25,50,75,100), col = "grey70", col.ticks = "grey70", col.axis = "grey70", )
  #axis(2, las = 1, at = c(-1,-.5,0,.5,1), col = "grey70", col.ticks = "grey70", col.axis = "grey70")
  #mtext(names(dat)[1], side = 4, xpd = T, at = ymin, col = "grey70", las = 1, adj = 0, padj = 1)
  #mtext(names(dat)[2], side = 3, xpd = T, at = range(gli.qn[,1])[1], col = "grey70", padj = 0, adj = 1)
  polygon(gli.qn[gch, ])
  polygon(li.qn[ch.qs, ], col = "grey60", xpd = T)
  abline(v = gcentroid[1], h = gcentroid[2])#, col = "grey70")
  arrows(x0 = 0, y0 = 0, x1 = co[, 1]/fact, y1 = co[, 2], xpd=T, length = .05)
  points(centroid[1], centroid[2], pch = 21, bg = "white")
  text(co[mags,1]*1.05/fact, co[mags,2]*1.05, labels = rownames(co[mags,]), cex = .5, xpd = T)
  legend("topright",
         legend = c(as.expression(bquote(xax == .(names(dat)[1]))), as.expression(bquote(yax == .(names(dat)[2]))), as.expression(bquote(~~~p == .(percentage)))),
         bty = "n")
}


# biplot <- function (x, xax = 1, yax = 2, pts = FALSE, nc = TRUE, percent = 95,
#                     clabel = 1, side = c("top", "bottom", "none"), Adensity,
#                     Udensity, Aangle, Uangle, Aborder, Uborder, Acol, Ucol, Alty,
#                     Ulty, Abg, Ubg, Ainch, Uinch, ...) {
#   gpres <- which(!is.na(values(max(raster(x)))))
#   li <- values(raster(x))[gpres,]
#   pr <- rep(NA, ncell(raster(x)))
#   pr[gpres] <- 0
#   pr[x@present] <- 1
#   pr <- pr[!is.na(pr)]
#
#   side <- match.arg(side)
#   if (!inherits(x, "cnfa"))
#     stop("Object of class 'cnfa' expected")
#   old.par <- par(no.readonly = TRUE)
#   on.exit(par(old.par))
#   par(mar = c(0.1, 0.1, 0.1, 0.1))
#   x1 <- li[, xax]
#   x1 <- c(x1 - diff(range(x1)/50), x1 + diff(range(x1))/50)
#   xlim <- range(x1)
#   y1 <- li[, yax]
#   y1 <- c(y1 - diff(range(y1)/50), y1 + diff(range(y1))/50)
#   ylim <- range(y1)
#   pmar <- t(x@mf) %*% as.matrix(x@co[, c(xax, yax)])
#   scatterutil.base(dfxy = li[, c(xax, yax)], xax = 1, yax = 2,
#                    xlim = xlim, ylim = ylim, grid = F, addaxes = T,
#                    cgrid = 1, include.origin = TRUE, origin = c(0, 0), sub = "",
#                    csub = 1.25, possub = "bottomleft", pixmap = NULL, contour = NULL,
#                    area = NULL, add.plot = FALSE)
#   if (pts) {
#     if (missing(Acol))
#       Acol <- gray(0.8)
#     if (missing(Ucol))
#       Ucol <- "black"
#     if (missing(Abg))
#       Abg <- gray(0.8)
#     if (missing(Ubg))
#       Ubg <- "black"
#     if (missing(Ainch))
#       Ainch <- 0.03
#     if (missing(Uinch))
#       Uinch <- Ainch #* max(x$pr)
#     symbols(li[, c(xax, yax)], circles = rep(1, length(pr)),
#             fg = Acol, bg = Abg, inches = Ainch, add = TRUE)
#     symbols(li[pr > 0, c(xax, yax)], circles = pr[pr >
#                                                     0], fg = Ucol, bg = Ubg, inches = Uinch, add = TRUE)
#     abline(v = 0)
#     abline(h = 0)
#     if (nc)
#       symbols(pmar, circles = 1, fg = "black", bg = "white",
#               inches = Ainch * 2, add = TRUE)
#   }
#   else {
#     if (missing(Adensity))
#       Adensity <- NULL
#     if (missing(Udensity))
#       Udensity <- NULL
#     if (missing(Aangle))
#       Aangle <- 45
#     if (missing(Uangle))
#       Uangle <- 45
#     if (missing(Aborder))
#       Aborder <- NULL
#     if (missing(Uborder))
#       Uborder <- NULL
#     if (missing(Acol))
#       Acol <- gray(0.95)
#     if (missing(Ucol))
#       Ucol <- gray(0.6)
#     if (missing(Alty))
#       Alty <- NULL
#     if (missing(Ulty))
#       Ulty <- NULL
#     pcff <- function(xy) {
#       mo <- apply(xy, 2, mean)
#       dis <- apply(xy, 1, function(x) sum((x - mo)^2))
#       xy <- xy[dis < quantile(dis, percent/100), ]
#       return(xy[chull(xy[, 1], xy[, 2]), ])
#     }
#     mcpA <- pcff(li[, c(xax, yax)])
#     pres <- which(pr > 0)
#     mcpU <- pcff(li[pres, c(xax, yax)])
#     polygon(mcpA, density = Adensity, angle = Aangle, border = Aborder,
#             col = Acol, lty = Alty)
#     polygon(mcpU, density = Udensity, angle = Uangle, border = Uborder,
#             col = Ucol, lty = Ulty)
#     abline(v = 0)
#     abline(h = 0)
#     if (nc)
#       points(pmar, pch = 21, bg = "white", cex = 1.5)
#   }
#   dfarr <- x@co[, c(xax, yax)]
#   born <- par("usr")
#   k1 <- min(dfarr[, 1])/born[1]
#   k2 <- max(dfarr[, 1])/born[2]
#   k3 <- min(dfarr[, 2])/born[3]
#   k4 <- max(dfarr[, 2])/born[4]
#   k <- c(k1, k2, k3, k4)
#   dfarr <- 0.75 * dfarr/max(k)
#   s.arrow(dfarr, clabel = clabel, addaxes = FALSE, add.plot = TRUE)
#   if (xax == 1)
#     xax <- "mar"
#   else xax <- paste("sp", xax - 1)
#   if (yax == 1)
#     yax <- "mar"
#   else yax <- paste("sp", yax - 1)
#   if (side != "none") {
#     tra <- paste(" xax =", xax, "\n yax =", yax)
#     wt <- strwidth(tra, cex = 1)
#     ht <- strheight(tra, cex = 1) * 1.5
#     xl <- par("usr")[1]
#     yu <- par("usr")[4]
#     yd <- par("usr")[3]
#     if (side == "top") {
#       rect(xl, yu - ht, xl + wt, yu, col = "white", border = 0)
#       text(xl + wt/2, yu - ht/2, tra, cex = 1)
#     }
#     if (side == "bottom") {
#       rect(xl, yd + ht, xl + wt, yd, col = "white", border = 0)
#       text(xl + wt/2, yd + ht/2, tra, cex = 1)
#     }
#   }
#   box()
# }
#
# scatterutil.base <- function (dfxy, xax, yax, xlim, ylim, grid, addaxes, cgrid, include.origin,
#                               origin, sub, csub, possub, pixmap, contour, area, add.plot) {
#   df <- data.frame(dfxy)
#   if (!is.data.frame(df))
#     stop("Non convenient selection for df")
#   if ((xax < 1) || (xax > ncol(df)))
#     stop("Non convenient selection for xax")
#   if ((yax < 1) || (yax > ncol(df)))
#     stop("Non convenient selection for yax")
#   x <- df[, xax]
#   y <- df[, yax]
#   if (is.null(xlim)) {
#     x1 <- x
#     if (include.origin)
#       x1 <- c(x1, origin[1])
#     x1 <- c(x1 - diff(range(x1)/10), x1 + diff(range(x1))/10)
#     xlim <- range(x1)
#   }
#   if (is.null(ylim)) {
#     y1 <- y
#     if (include.origin)
#       y1 <- c(y1, origin[2])
#     y1 <- c(y1 - diff(range(y1)/10), y1 + diff(range(y1))/10)
#     ylim <- range(y1)
#   }
#   if (!is.null(pixmap)) {
#     if (is.null(class(pixmap)))
#       pixmap <- NULL
#     if (is.na(charmatch("pixmap", class(pixmap))))
#       pixmap <- NULL
#   }
#   if (!is.null(contour)) {
#     if (!is.data.frame(contour))
#       contour <- NULL
#     if (ncol(contour) != 4)
#       contour <- NULL
#   }
#   if (!is.null(area)) {
#     if (!is.data.frame(area))
#       area <- NULL
#     if (!is.factor(area[, 1]))
#       area <- NULL
#     if (ncol(area) < 3)
#       area <- NULL
#   }
#   if (!add.plot)
#     plot.default(0, 0, type = "n", xlab = "", ylab = "",
#                  xaxt = "n", yaxt = "n", xlim = xlim, ylim = ylim,
#                  xaxs = "i", yaxs = "i", frame.plot = FALSE)
#   if (!is.null(pixmap)) {
#     pixmap::plot(pixmap, add = TRUE)
#   }
#   if (!is.null(contour)) {
#     apply(contour, 1, function(x) segments(x[1], x[2], x[3],
#                                            x[4], lwd = 1))
#   }
#   if (grid & !add.plot)
#     scatterutil.grid(cgrid)
#   if (addaxes & !add.plot)
#     abline(h = 0, v = 0, lty = 1)
#   if (!is.null(area)) {
#     nlev <- nlevels(area[, 1])
#     x1 <- area[, 2]
#     x2 <- area[, 3]
#     for (i in 1:nlev) {
#       lev <- levels(area[, 1])[i]
#       a1 <- x1[area[, 1] == lev]
#       a2 <- x2[area[, 1] == lev]
#       polygon(a1, a2)
#     }
#   }
#   if (csub > 0)
#     scatterutil.sub(sub, csub, possub)
#   return(list(x = x, y = y))
# }
#
# scatterutil.sub <- function (cha, csub, possub = "bottomleft") {
#   cha <- as.character(cha)
#   if (length(cha) == 0)
#     return(invisible())
#   if (is.null(cha))
#     return(invisible())
#   if (is.na(cha))
#     return(invisible())
#   if (any(cha == ""))
#     return(invisible())
#   if (csub == 0)
#     return(invisible())
#   cex0 <- par("cex") * csub
#   cha <- paste(" ", cha, " ", sep = "")
#   xh <- strwidth(cha, cex = cex0)
#   yh <- strheight(cha, cex = cex0) * 5/3
#   if (possub == "bottomleft") {
#     x1 <- par("usr")[1]
#     y1 <- par("usr")[3]
#     rect(x1, y1, x1 + xh, y1 + yh, col = "white", border = 0)
#     text(x1 + xh/2, y1 + yh/2, cha, cex = cex0)
#   }
#   else if (possub == "topleft") {
#     x1 <- par("usr")[1]
#     y1 <- par("usr")[4]
#     rect(x1, y1, x1 + xh, y1 - yh, col = "white", border = 0)
#     text(x1 + xh/2, y1 - yh/2, cha, cex = cex0)
#   }
#   else if (possub == "bottomright") {
#     x1 <- par("usr")[2]
#     y1 <- par("usr")[3]
#     rect(x1, y1, x1 - xh, y1 + yh, col = "white", border = 0)
#     text(x1 - xh/2, y1 + yh/2, cha, cex = cex0)
#   }
#   else if (possub == "topright") {
#     x1 <- par("usr")[2]
#     y1 <- par("usr")[4]
#     rect(x1, y1, x1 - xh, y1 - yh, col = "white", border = 0)
#     text(x1 - xh/2, y1 - yh/2, cha, cex = cex0)
#   }
# }
#
# scatterutil.grid <- function (cgrid) {
#   col <- "lightgray"
#   lty <- 1
#   xaxp <- par("xaxp")
#   ax <- (xaxp[2] - xaxp[1])/xaxp[3]
#   yaxp <- par("yaxp")
#   ay <- (yaxp[2] - yaxp[1])/yaxp[3]
#   a <- min(ax, ay)
#   v0 <- seq(xaxp[1], xaxp[2], by = a)
#   h0 <- seq(yaxp[1], yaxp[2], by = a)
#   abline(v = v0, col = col, lty = lty)
#   abline(h = h0, col = col, lty = lty)
#   if (cgrid <= 0)
#     return(invisible())
#   cha <- paste(" d = ", a, " ", sep = "")
#   cex0 <- par("cex") * cgrid
#   xh <- strwidth(cha, cex = cex0)
#   yh <- strheight(cha, cex = cex0) * 5/3
#   x1 <- par("usr")[2]
#   y1 <- par("usr")[4]
#   rect(x1 - xh, y1 - yh, x1 + xh, y1 + yh, col = "white", border = 0)
#   text(x1 - xh/2, y1 - yh/2, cha, cex = cex0)
# }
#
# s.arrow <- function (dfxy, xax = 1, yax = 2, label = row.names(dfxy), clabel = 1,
#                      pch = 20, cpoint = 0, boxes = TRUE, edge = TRUE, origin = c(0,
#                                                                                  0), xlim = NULL, ylim = NULL, grid = TRUE, addaxes = TRUE,
#                      cgrid = 1, sub = "", csub = 1.25, possub = "bottomleft",
#                      pixmap = NULL, contour = NULL, area = NULL, add.plot = FALSE) {
#   arrow1 <- function(x0, y0, x1, y1, len = 0.1, ang = 15, lty = 1,
#                      edge) {
#     d0 <- sqrt((x0 - x1)^2 + (y0 - y1)^2)
#     if (d0 < 1e-07)
#       return(invisible())
#     segments(x0, y0, x1, y1, lty = lty)
#     h <- strheight("A", cex = par("cex"))
#     if (d0 > 2 * h) {
#       x0 <- x1 - h * (x1 - x0)/d0
#       y0 <- y1 - h * (y1 - y0)/d0
#       if (edge)
#         arrows(x0, y0, x1, y1, angle = ang, length = len,
#                lty = 1)
#     }
#   }
#   dfxy <- data.frame(dfxy)
#   opar <- par(mar = par("mar"))
#   on.exit(par(opar))
#   par(mar = c(0.1, 0.1, 0.1, 0.1))
#   coo <- scatterutil.base(dfxy = dfxy, xax = xax, yax = yax,
#                           xlim = xlim, ylim = ylim, grid = grid, addaxes = addaxes,
#                           cgrid = cgrid, include.origin = TRUE, origin = origin,
#                           sub = sub, csub = csub, possub = possub, pixmap = pixmap,
#                           contour = contour, area = area, add.plot = add.plot)
#   if (grid & !add.plot)
#     scatterutil.grid(cgrid)
#   if (addaxes & !add.plot)
#     abline(h = 0, v = 0, lty = 1)
#   if (cpoint > 0)
#     points(coo$x, coo$y, pch = pch, cex = par("cex") * cpoint)
#   for (i in 1:(length(coo$x))) arrow1(origin[1], origin[2],
#                                       coo$x[i], coo$y[i], edge = edge)
#   if (clabel > 0)
#     scatterutil.eti.circ(coo$x, coo$y, label, clabel, origin,
#                          boxes)
#   if (csub > 0)
#     scatterutil.sub(sub, csub, possub)
#   box()
#   invisible(match.call())
# }
#
# scatterutil.eti.circ <- function (x, y, label, clabel, origin = c(0, 0), boxes = TRUE) {
#   if (is.null(label))
#     return(invisible())
#   if (any(is.na(label)))
#     return(invisible())
#   if (any(label == ""))
#     return(invisible())
#   xref <- x - origin[1]
#   yref <- y - origin[2]
#   for (i in 1:(length(x))) {
#     cha <- as.character(label[i])
#     cha <- paste(" ", cha, " ", sep = "")
#     cex0 <- par("cex") * clabel
#     xh <- strwidth(cha, cex = cex0)
#     yh <- strheight(cha, cex = cex0) * 5/6
#     if ((xref[i] > yref[i]) & (xref[i] > -yref[i])) {
#       x1 <- x[i] + xh/2
#       y1 <- y[i]
#     }
#     else if ((xref[i] > yref[i]) & (xref[i] <= (-yref[i]))) {
#       x1 <- x[i]
#       y1 <- y[i] - yh
#     }
#     else if ((xref[i] <= yref[i]) & (xref[i] <= (-yref[i]))) {
#       x1 <- x[i] - xh/2
#       y1 <- y[i]
#     }
#     else if ((xref[i] <= yref[i]) & (xref[i] > (-yref[i]))) {
#       x1 <- x[i]
#       y1 <- y[i] + yh
#     }
#     if (boxes) {
#       rect(x1 - xh/2, y1 - yh, x1 + xh/2, y1 + yh, col = "white",
#            border = 1)
#     }
#     text(x1, y1, cha, cex = cex0)
#   }
# }

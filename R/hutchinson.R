remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

hutchinson<-function (x, xax = 1, yax = 2, pts = FALSE, nc = TRUE, percent = 95,
          clabel = 1, side = c("top", "bottom", "none"), Adensity,
          Udensity, Aangle, Uangle, Aborder, Uborder, Acol, Ucol, Alty,
          Ulty, Abg, Ubg, Ainch, Uinch, ...)
{
  side <- match.arg(side)
  if (!inherits(x, "enfa"))
    stop("Object of class 'enfa' expected")
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  par(mar = c(1, 1, 1, 1))
  x1 <- SpatialPoints(na.omit(values(x@global_ras)[,c(xax,yax)]))
  p1<-mcp(x1,percent=90)
  x2 <- SpatialPoints(na.omit(values(x@species_ras)[,c(xax,yax)]))
  p2<-mcp(x2,percent=90)
   y1 <- SpatialPoints(x@global_ras[[yax]])))
  dists<-sqrt(x1^2+y1^2)
  dists2<-remove_outliers(dists)
  temp<-which(is.na(dists2))
  x1<-x1[-temp]
  xlim <- range(x1)
  y1<-y1[-temp]
  ylim <- range(y1)
  x2 <- c(na.omit(values(x@species_ras[[xax]])))
  y2 <- c(na.omit(values(x@species_ras[[yax]])))
  dists<-sqrt(x2^2+y2^2)
  dists2<-remove_outliers(dists)
  temp<-which(is.na(dists2))
  x2<-x2[-temp]
  xlim <- range(x1)
  y2<-y2[-temp]
  ylim <- range(y1)

  #pmar <- t(x$mar * x$cw) %*% as.matrix(x$co[, c(xax, yax)])
  dat<-cbind(x1,y1)
  plot(NA,xlim=xlim,ylim=ylim,ann=F)
  plot(p1,add=T)
  plot(p2,add=T,col="darkgreen")
  mcp1<-chull(x1,y1)
  polygon(x1[mcp1],y1[mcp1])
  sdat<-cbind(x2,y2)
  mcp2<-chull(x2,y2)
  polygon(x2[mcp2],y2[mcp2],col="red")

  samp<-sample(nrow(sdat),1000)
  points(dat[samp,],pch=".")
  pcff <- function(xy) {
    mo <- apply(xy, 2, mean)
    dis <- apply(xy, 1, function(x) sum((x - mo)^2))
    xy <- xy[dis < quantile(dis, percent/100), ]
    return(xy[chull(xy[, 1], xy[, 2]), ])
  }
  mcpA <- pcff(x$li[, c(xax, yax)])
  mcpU <- pcff(x$li[rep(1:length(x$pr), x$pr), c(xax, yax)])
  polygon(mcpA, density = Adensity, angle = Aangle, border = Aborder,
          col = Acol, lty = Alty)
  ade4::scatterutil.base(dfxy = dat, xax = 1, yax = 2,
                   xlim = xlim, ylim = ylim, grid = TRUE, addaxes = FALSE,
                   cgrid = 1, include.origin = TRUE, origin = c(0, 0), sub = "",
                   csub = 1.25, possub = "bottomleft", pixmap = NULL, contour = NULL,
                   area = NULL, add.plot = FALSE)
  if (pts) {
    if (missing(Acol))
      Acol <- gray(0.8)
    if (missing(Ucol))
      Ucol <- "black"
    if (missing(Abg))
      Abg <- gray(0.8)
    if (missing(Ubg))
      Ubg <- "black"
    if (missing(Ainch))
      Ainch <- 0.03
    if (missing(Uinch))
      Uinch <- Ainch * max(x$pr)
    symbols(x$li[, c(xax, yax)], circles = rep(1, length(x$pr)),
            fg = Acol, bg = Abg, inches = Ainch, add = TRUE)
    symbols(x$li[x$pr > 0, c(xax, yax)], circles = x$pr[x$pr >
                                                          0], fg = Ucol, bg = Ubg, inches = Uinch, add = TRUE)
    abline(v = 0)
    abline(h = 0)
    if (nc)
      symbols(pmar, circles = 1, fg = "black", bg = "white",
              inches = Ainch * 2, add = TRUE)
  }
  else {
    if (missing(Adensity))
      Adensity <- NULL
    if (missing(Udensity))
      Udensity <- NULL
    if (missing(Aangle))
      Aangle <- 45
    if (missing(Uangle))
      Uangle <- 45
    if (missing(Aborder))
      Aborder <- NULL
    if (missing(Uborder))
      Uborder <- NULL
    if (missing(Acol))
      Acol <- gray(0.95)
    if (missing(Ucol))
      Ucol <- gray(0.6)
    if (missing(Alty))
      Alty <- NULL
    if (missing(Ulty))
      Ulty <- NULL
    pcff <- function(xy) {
      mo <- apply(xy, 2, mean)
      dis <- apply(xy, 1, function(x) sum((x - mo)^2))
      xy <- xy[dis < quantile(dis, percent/100), ]
      return(xy[chull(xy[, 1], xy[, 2]), ])
    }
    mcpA <- pcff(x$li[, c(xax, yax)])
    mcpU <- pcff(x$li[rep(1:length(x$pr), x$pr), c(xax, yax)])
    polygon(mcpA, density = Adensity, angle = Aangle, border = Aborder,
            col = Acol, lty = Alty)
    polygon(mcpU, density = Udensity, angle = Uangle, border = Uborder,
            col = Ucol, lty = Ulty)
    abline(v = 0)
    abline(h = 0)
    if (nc)
      points(pmar, pch = 21, bg = "white", cex = 1.5)
  }
  dfarr <- x$co[, c(xax, yax)]
  born <- par("usr")
  k1 <- min(dfarr[, 1])/born[1]
  k2 <- max(dfarr[, 1])/born[2]
  k3 <- min(dfarr[, 2])/born[3]
  k4 <- max(dfarr[, 2])/born[4]
  k <- c(k1, k2, k3, k4)
  dfarr <- 0.75 * dfarr/max(k)
  s.arrow(dfarr, clabel = clabel, addaxes = FALSE, add.plot = TRUE)
  if (xax == 1)
    xax <- "mar"
  else xax <- paste("sp", xax - 1)
  if (yax == 1)
    yax <- "mar"
  else yax <- paste("sp", yax - 1)
  if (side != "none") {
    tra <- paste(" xax =", xax, "\n yax =", yax)
    wt <- strwidth(tra, cex = 1)
    ht <- strheight(tra, cex = 1) * 1.5
    xl <- par("usr")[1]
    yu <- par("usr")[4]
    yd <- par("usr")[3]
    if (side == "top") {
      rect(xl, yu - ht, xl + wt, yu, col = "white", border = 0)
      text(xl + wt/2, yu - ht/2, tra, cex = 1)
    }
    if (side == "bottom") {
      rect(xl, yd + ht, xl + wt, yd, col = "white", border = 0)
      text(xl + wt/2, yd + ht/2, tra, cex = 1)
    }
  }
  box()
}

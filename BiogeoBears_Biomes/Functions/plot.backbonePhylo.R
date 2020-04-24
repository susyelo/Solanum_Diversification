scaleN<-function (x, k) 
{
  for (i in 1:length(x$tip.clade)) if (x$tip.clade[[i]]$N > 
                                       1) 
    x$tip.clade[[i]]$N <- x$tip.clade[[i]]$N * k
  x
}

plot.backbonePhylo<-function (x, ...) 
{
  if (!inherits(x, "backbonePhylo")) 
    stop("x not an object of class \"backbonePhylo\"")
  if (hasArg(vscale)) 
    vscale <- list(...)$vscale
  else vscale <- 1
  if (hasArg(col)) 
    col <- list(...)$col
  else col <- "grey"
  if (length(col) != Ntip(x)) {
    if (!is.null(names(col))) {
      tmp <- setNames(rep("grey", Ntip(x)), sapply(x$tip.clade, 
                                                   function(x) x$label))
      tmp[names(col)] <- col
      col <- tmp
    }
    else col <- setNames(rep(col[1], Ntip(x)), sapply(x$tip.clade, 
                                                      function(x) x$label))
  }
  if (is.null(names(col))) 
    names(col) <- sapply(x$tip.clade, function(x) x$label)
  x <- scaleN(x, vscale)
  tt <- backbone.toPhylo(x)
  n <- sum(sapply(x$tip.clade, function(x) x$N))
  cw <- reorder.backbonePhylo(x, "cladewise")
  y <- vector(length = length(cw$tip.clade) + cw$Nnode)
  ii <- order(cw$edge[, 2][cw$edge[, 2] <= Ntip(cw)])
  z <- c(0, cumsum(sapply(cw$tip.clade[order(ii)], function(x) x$N)))
  nn <- sapply(2:length(z), function(i, x) (x[i] - x[i - 1])/2 + 
                 x[i - 1], x = z)
  y[cw$edge[cw$edge[, 2] <= length(cw$tip.clade), 2]] <- nn[1:length(cw$tip.clade)]
  pw <- reorder.backbonePhylo(x, "pruningwise")
  nn <- unique(pw$edge[, 1])
  for (i in 1:length(nn)) {
    yy <- y[pw$edge[which(pw$edge[, 1] == nn[i]), 2]]
    y[nn[i]] <- mean(range(yy))
  }
  X <- nodeHeights(tt)
  plot.new()
  par(mar = rep(0.1, 4))
  pp <- par("pin")[1]
  sw <- par("cex") * (max(strwidth(sapply(cw$tip.clade, function(x) x$label), 
                                   units = "inches"))) + 1.37 * par("cex") * strwidth("W", 
                                                                                      units = "inches")
  alp <- optimize(function(a, H, sw, pp) (a * 1.04 * max(H) + 
                                            sw - pp)^2, H = X, sw = sw, pp = pp, interval = c(0, 
                                                                                              1e+06))$minimum
  xlim <- c(min(X), max(X) + sw/alp)
  plot.window(xlim = xlim, ylim = c(0, n))
  for (i in 1:nrow(X)) {
    if (cw$edge[i, 2] > length(cw$tip.clade)) 
      lines(X[i, ], rep(y[cw$edge[i, 2]], 2), lwd = 2, 
            lend = 2)
    else lines(c(X[i, 1], X[i, 2] - cw$tip.clade[[cw$edge[i, 
                                                          2]]]$depth), rep(y[cw$edge[i, 2]], 2), lwd = 2, lend = 2)
  }
  for (i in 1:x$Nnode + length(x$tip.clade)) lines(X[which(cw$edge[, 
                                                                   1] == i), 1], range(y[cw$edge[which(cw$edge[, 1] == i), 
                                                                                                 2]]), lwd = 2, lend = 2)
  for (i in 1:length(x$tip.clade)) {
    xx <- c(X[which(cw$edge[, 2] == i), 2] - cw$tip.clade[[i]]$depth, 
            X[which(cw$edge[, 2] == i), 2], X[which(cw$edge[, 
                                                            2] == i), 2])
    yy <- c(y[cw$edge[which(cw$edge[, 2] == i), 2]], y[cw$edge[which(cw$edge[, 
                                                                             2] == i), 2]] + cw$tip.clade[[i]]$N/2 - 0.5, y[cw$edge[which(cw$edge[, 
                                                                                                                                                  2] == i), 2]] - cw$tip.clade[[i]]$N/2 + 0.5)
    if (yy[2] < yy[3]) 
      yy[2] <- yy[3] <- yy[1]
    polygon(x = xx, y = yy, col = col[cw$tip.clade[[i]]$label], 
            lwd = 2)
  }
  for (i in 1:length(cw$tip.clade)) text(X[which(cw$edge[, 
                                                         2] == i), 2], y[i], cw$tip.clade[[i]]$label, pos = 4, 
                                         offset = 0.1)
  PP <- list(type = "phylogram", use.edge.length = TRUE, node.pos = 1, 
             show.tip.label = TRUE, show.node.label = FALSE, font = 1, 
             cex = 1, adj = 0, srt = 0, no.margin = FALSE, label.offset = 0.1, 
             x.lim = par()$usr[1:2], y.lim = par()$usr[3:4], direction = "rightwards", 
             tip.color = "black", Ntip = Ntip(cw), Nnode = cw$Nnode, 
             edge = cw$edge, xx = sapply(1:(Ntip(cw) + cw$Nnode), 
                                         function(x, y, z) y[match(x, z)], y = X, z = cw$edge), 
             yy = y)
  assign("last_plot.phylo", PP, envir = .PlotPhyloEnv)
  
}
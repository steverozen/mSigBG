
foo.catalog <- # Update of ICAMS::PlotCatalog.SBS96Catalog
function (catalog, plot.SBS12, cex = par("cex"), grid = TRUE, 
          upper = TRUE, xlabels = TRUE, ylabels = TRUE, ylim = NULL) 
{
  stopifnot(rownames(catalog) == ICAMS::catalog.row.order$SBS96)
  class.col <- c("#0000ff", "#000000", "#ff4040", "#838383", 
                 "#40ff40", "#ff667f")
  if (ncol(catalog) == 1) {
    cols <- rep(class.col, each = 16)
    to.plot <- catalog[, 1]
  } else {
    # stop("Can only handle 1 or 2 column catalogs")
    num.colors <- min(ncol(catalog), 8)
    cols <- RColorBrewer::brewer.pal(num.colors, "Dark2")
    to.plot <- t(catalog)
  }
  maj.class.names <- c("C>A", "C>G", "C>T", "T>A", "T>C", 
                       "T>G")
  num.classes <- 96
  catalog.type <- attributes(catalog)$catalog.type
  if (catalog.type == "density") {
    ylab <- "mut/million"
    to.plot <- 1e+06 * to.plot
    ymax <- 1e+06 * max(rowSums(catalog))
  }
  else if (catalog.type == "counts") {
    ymax <- 4 * ceiling(max(max(rowSums(catalog)), 10)/4)
    ylab <- "counts"
  }
  else if (catalog.type %in% c("counts.signature", "density.signature")) {
    ylab <- ifelse(catalog.type == "counts.signature", "counts proportion", 
                   "density proportion")
    ymax <- max(rowSums(catalog))
  }
  else {
    stop("Programming error, illegal catalog type ", catalog.type)
  }
  if (is.null(ylim)) {
    ylim <- c(0, ymax)
  }
  else {
    ymax <- ylim[2]
  }
  if (ylabels == FALSE) {
    bp <- barplot(to.plot, xaxt = "n", yaxt = "n", xaxs = "i", 
                  xlim = c(-1, 230), ylim = ylim, lwd = 3, space = 1.35, 
                  border = NA, col = cols, cex.lab = cex * par("cex.lab"))
  }
  else {
    bp <- barplot(to.plot, xaxt = "n", yaxt = "n", xaxs = "i", 
                  xlim = c(-1, 230), ylim = ylim, lwd = 3, space = 1.35, 
                  border = NA, col = cols, ylab = ylab, cex.lab = cex * 
                    par("cex.lab"))
  }
  segments(bp[1] - 0.5, 0, bp[num.classes] + 0.5, 0, col = "grey35", 
           lwd = 0.25)
  y.pos.sample.name <- ifelse(catalog.type == "counts", ymax * 
                                1.22, ymax * 1.08)
  text(bp[1], y.pos.sample.name, labels = colnames(catalog)[ncol(catalog)], 
       xpd = NA, cex = cex, font = 2, adj = c(0, 0))
  if (catalog.type == "counts") {
    count.cex <- cex
    for (i in 1:6) {
      j <- 16 + 16 * (i - 1)
      k <- 1 + 16 * (i - 1)
      text(bp[j], ymax * 1.15, labels = round(sum(abs(catalog[k:(16 * 
                                                                   i), ]))), adj = c(1, 1), xpd = NA, cex = count.cex)
    }
  }
  y.axis.values <- seq(0, ymax, ymax/4)
  if (catalog.type != "counts") {
    y.axis.labels <- format(round(y.axis.values, 2), nsmall = 2)
  }
  else {
    y.axis.labels <- y.axis.values
  }
  if (grid == TRUE) {
    segments(bp[1] - 0.5, seq(ymax/4, ymax, ymax/4), bp[num.classes] + 
               0.5, seq(ymax/4, ymax, ymax/4), col = "grey35", 
             lwd = 0.25)
    if (ylabels == TRUE) {
      text(-0.5, y.axis.values, labels = y.axis.labels, 
           las = 1, adj = 1, xpd = NA, cex = cex)
    }
  }
  else {
    if (ylabels == TRUE) {
      Axis(side = 2, at = y.axis.values, las = 1, cex.axis = cex, 
           labels = FALSE)
      text(-3.5, y.axis.values, labels = y.axis.labels, 
           cex = cex, las = 1, adj = 1, xpd = NA)
    }
  }
  if (xlabels == TRUE) {
    cex.xlabel <- cex
    xlabel.idx <- seq(1, 96, by = 4)
    label <- c("A", "C", "G", "T")
    text(bp[xlabel.idx], -ymax/7, labels = label, cex = cex.xlabel, 
         adj = 0.5, xpd = NA)
    x <- list(bp[xlabel.idx], bp[xlabel.idx + 1], bp[xlabel.idx + 
                                                       2], bp[xlabel.idx + 3])
    y <- c(-ymax/3.5, -ymax/2.8, -ymax/2.39, -ymax/2.1)
    for (i in 1:4) {
      text(x[[i]], y[i], labels = label[i], cex = cex.xlabel, 
           adj = 0.5, xpd = NA)
    }
    text(1.5, -ymax/7, labels = "preceded by 5'", pos = 2, 
         xpd = NA, cex = cex.xlabel)
    text(1.5, -ymax/3.5, labels = "followed by 3'", pos = 2, 
         xpd = NA, cex = cex.xlabel)
  }
  else {
    every.fourth <- seq(from = 1, to = length(bp), by = 4)
    every.fourth <- c(every.fourth, 96)
    Axis(at = bp[every.fourth], side = 1, labels = FALSE, 
         col = "grey35")
  }
  if (upper == TRUE) {
    x.left <- bp[seq(1, 81, 16)]
    x.right <- bp[seq(16, 96, 16)]
    y.pos.line.top <- ifelse(catalog.type == "counts", ymax * 
                               1.4, ymax * 1.3)
    y.pos.line.bottom <- ifelse(catalog.type == "counts", 
                                ymax * 1.38, ymax * 1.28)
    rect(xleft = x.left, ybottom = y.pos.line.bottom, xright = x.right, 
         ytop = y.pos.line.top, col = class.col, border = NA, 
         xpd = NA, adj = 0.5)
    y.pos.text <- ifelse(catalog.type == "counts", ymax * 
                           1.48, ymax * 1.38)
    text((x.left + x.right)/2, y.pos.text, labels = maj.class.names, 
         xpd = NA, cex = cex * 1.25)
  }
  invisible(list(plot.success = TRUE, plot.object = bp))
}

#' Plot a spectrum as a stacked bar chart with component signatures
#' 
#' @param sigs An SBS96 catalog.
#' 
#' @param exposure A vector of length `(ncol(sigs))`.
#' 
#' @export
#' 
#' @import graphics, RColorBrewer
#' 
# plot_stacked_sigs_by_exposure
plot_stacked_sigs_by_exposure <- function(sigs, exposures) {
  oo <- order(exposures, decreasing = T)
  sigs <- sigs[ , oo, drop = FALSE]
  exposures <- exposures[oo]
  
  mm <- rep(exposures, 96)
  
  scaled.sigs <- mm * sigs
  
  bp <- foo.catalog(scaled.sigs)$plot.object
  
  legend.txt <- rev(colnames(sigs))
  num.colors <- min(ncol(sigs), 8)
  legend.col <-     # stop("Can only handle 1 or 2 column catalogs")
    rev(RColorBrewer::brewer.pal(num.colors, "Dark2"))
  legend("topright", legend = legend.txt, fill = legend.col, cex = 0.8)

  return(invisible(bp))
}

if (FALSE) {
  vv <- cosmicsig::COSMIC_v3.2$signature$GRCh37$SBS96[ , 1:6]
  zz(vv, c(0.3, 0.25, 0.2, 0.15, 0.1, 0.05))
}

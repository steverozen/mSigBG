one.col <- function(col) {
  xx <- col / sum(col)
  return(xx)
}

yy <- apply(nitrosamines, MARGIN = 2, one.col)
attr(yy, "catalog.type") <- "counts.signature"

ICAMS::PlotCatalogToPdf(yy, "foo.pdf")

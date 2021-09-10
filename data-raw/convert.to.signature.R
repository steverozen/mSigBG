#import draw_pca and check_pca_quality function from test_pca.R
data.dir <- file.path('data-raw')
source(paste0(data.dir,'/test_pca.R'))

one.col <- function(col) {
  mean_count <- col / sum(col)
  return(mean_count)
}

sbs96.nitrosamine.meancount <- apply(nitrosamines, MARGIN = 2, one.col)
attr(sbs96.nitrosamine.meancount, "catalog.type") <- "counts.signature"
attr(sbs96.nitrosamine.meancount, "class") <- attr(nitrosamines, "class")

draw_pca(sbs96.nitrosamine.meancount)
check_pca_quality(sbs96.nitrosamine.meancount)
draw_pca(nitrosamines)
check_pca_quality(nitrosamines)

ICAMS::PlotCatalogToPdf(sbs96.nitrosamine.meancount, "sbs96.nitrosamine.meancount.pdf")

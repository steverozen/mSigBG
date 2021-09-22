library(PCAWG7)
library(tidyr)
library(philentropy)
# Example to generate one synthetic spectrum due to HepG2 background and
# SBS22 (AA, easy) SBS4 (smoking, hard example)

# 1. Generate partial spectrum for background

# We use info in mSigBG::background.info
set.seed(101010, kind = "L'Ecuyer-CMRG")

h2 <- mSigBG::background.info$HepG2
target <- "SBS22"
target.sig <- PCAWG7::signature$genome$SBS96[ , target, drop = FALSE]
target.exp <- PCAWG7::exposure$other.genome$SBS96[target,,drop = FALSE]
target.exp <- target.exp[target.exp > 0.5]
target.mean <- mean(target.exp)

placeholder <- target.sig
attr(placeholder,'catalog.type') <- 'counts'
plotlist <- placeholder#place holder, dont know how to generate empty catalog
diff.muts <- vector()
output.syn.spec <- placeholder

for(i in 1:5){
  num.muts <- rnbinom(1, size = h2$count.nbinom.size, mu = h2$count.nbinom.mu)
  syn.target.num.muts <- num.muts
  diff.muts <- c(diff.muts, num.muts)
  print(num.muts)
  ideal.bg.spec <- h2$background.sig * num.muts
  ideal.target.spec <- target.sig * syn.target.num.muts
  ideal.syn.spec <- ideal.bg.spec + ideal.target.spec
  attr(ideal.bg.spec, "catalog.type") <- "counts"
  attr(ideal.target.spec, "catalog.type") <- "counts"
  attr(ideal.target.spec, "class") <- attr(ideal.bg.spec, "class")
  attr(ideal.syn.spec, "catalog.type") <- "counts"
  attr(ideal.syn.spec, "class") <- attr(ideal.bg.spec, "class")
  colnames(ideal.bg.spec) <- paste0("HepG2.background",'.cl',i)
  colnames(ideal.target.spec) <- paste0("SBS22",'.cl',i)
  colnames(ideal.syn.spec) <- paste0("ideal.combined",'.cl',i)
  plotlist <- cbind(plotlist, 
                    ideal.bg.spec,ideal.target.spec,ideal.syn.spec)
  
  bg.partial.spectrum.w.noise <- unlist(lapply(ideal.bg.spec,  function(x) rnbinom(1, mu = x, size = 10)))
  target.partial.spectrum.w.noise  <- unlist(lapply(ideal.target.spec,  function(x) rnbinom(1, mu = x, size = 10)))
  syn.spec <- bg.partial.spectrum.w.noise + target.partial.spectrum.w.noise
  
  bg.partial.spectrum.w.noise <- as.matrix(bg.partial.spectrum.w.noise)
  row.names(bg.partial.spectrum.w.noise) <- rownames(ideal.bg.spec)
  target.partial.spectrum.w.noise <- as.matrix(target.partial.spectrum.w.noise)
  row.names(target.partial.spectrum.w.noise) <- rownames(ideal.target.spec)
  target.partial.spectrum.w.noise <- ICAMS::as.catalog(target.partial.spectrum.w.noise)
  syn.spec <- as.matrix(syn.spec)
  row.names(syn.spec) <- rownames(ideal.syn.spec)
  syn.spec <- ICAMS::as.catalog(syn.spec)
  colnames(bg.partial.spectrum.w.noise) <- paste0("HepG2.background.w/noise",'.cl',i)
  colnames(target.partial.spectrum.w.noise) <- paste0("SBS22.w/noise",'.cl',i)
  colnames(syn.spec) <- paste0("combined.w/noise",'.cl',i)
  attr(bg.partial.spectrum.w.noise, "catalog.type") <- "counts"
  attr(target.partial.spectrum.w.noise, "catalog.type") <- "counts"
  attr(syn.spec,'catalog.type') <- 'counts'
  plotlist <- cbind(plotlist, 
                    bg.partial.spectrum.w.noise, target.partial.spectrum.w.noise, syn.spec)
  output.syn.spec <- cbind(output.syn.spec, syn.spec)
  
}
ICAMS::PlotCatalogToPdf(plotlist[,-1],file='specs_for_5cls.pdf')
attr(output.syn.spec,'abundance') <- NULL
attr(output.syn.spec,'region') <- 'unknown'
ret<- 
  SeparateSignatureFromBackground(
    spectra     = output.syn.spec[,-1],
    bg.sig.info = h2,
    m.opts      = c(NULL,NULL,NULL,NULL,NULL),
    start.b.fraction = 0.5)

inferred.sig <- 
  ICAMS::as.catalog(
    object         = ret$inferred.target.sig, 
    catalog.type   = "counts.signature",
    region         = "genome",
    abundance      = attr(h2$background.sig,
                          "abundance"),
    infer.rownames = TRUE)
#draw inferred target signature spectrum
colnames(inferred.sig) <- paste0("Inferred-", target, "-sig")
two.sigs <- cbind(inferred.sig, h2$background.sig)
ICAMS::PlotCatalogToPdf(two.sigs,file='inferred_sig.pdf')

#draw inferred bg and target counted spectrum
inferred.bg.spectra <- 
  round(h2$background.sig %*% matrix(ret$exposures.to.bg.sig, nrow = 1))
inferred.bg.spectra <- ICAMS::as.catalog(inferred.bg.spectra,
                                         catalog.type   = "counts",
                                         region         = "genome",
                                         abundance = attr(ideal.syn.spec, "abundance"),
                                         infer.rownames = TRUE)
colnames(inferred.bg.spectra) <- paste0(sub('.*([1-9]).*','cl\\1',colnames(output.syn.spec[,-1])), "-inf-bg-spect")

inferred.target.spectra <- ideal.syn.spec[,rep(1,time =5)] - inferred.bg.spectra
colnames(inferred.target.spectra) <- paste0(sub('.*([1-9]).*','cl\\1',colnames(output.syn.spec[,-1])), "-inf-target-spect")
ICAMS::PlotCatalogToPdf(cbind(inferred.bg.spectra,inferred.target.spectra),file = 'inferred_count.pdf')

#draw scatter plots for count diff between inferred num of muts and ideal num of muts
inferred.muts <- cbind(ret$exposures.to.target.sig, ret$exposures.to.bg.sig, diff.muts, diff.muts)
colnames(inferred.muts) <- c('infer.target','infer.bg','assigned.target','assigned.bg')
inferred.muts <- data.frame(names = sub('.*([1-9]).*','cl\\1',rownames(inferred.muts)),inferred.muts)
rownames(inferred.muts) <- NULL;inferred.muts
inferred.muts$sum.inferred <- rowSums(cbind(inferred.muts$infer.target, inferred.muts$infer.bg))
inferred.muts$sum.ideal <- rowSums(cbind(inferred.muts$assigned.target, inferred.muts$assigned.bg))

forscatterplot = inferred.muts%>% #dataframe row to column transformation
  pivot_longer(cols = !1,
               names_to = "source",
               values_to = "count");

p <- ggplot(data = forscatterplot) +
  geom_point(mapping = aes(x=names,y=count,color=source,size = 1))+
  ggtitle('Diff between inferred count and ideal count')+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))+x
  guides(size='none')#hide legend 'size'
p

#cosine distance between target signatures
#cosine.dis.target <- lsa::cosine(inferred.sig[,1],target.sig[,1])
cosine.dis.target <- philentropy::cosine_dist(inferred.sig[,1], target.sig[,1],testNA=FALSE)
cosine.dis.bg <- philentropy::cosine_dist(h2$background.sig[,1], 
                             mSigBG::MeanOfSpectraAsSig(inferred.bg.spectra)[,1],testNA = FALSE)

#draw uncertainty
sub.target.sig.ests <- 
  ICAMS::TransformCatalog(inferred.target.spectra, 
                          target.catalog.type = "counts.signature",
                          target.abundance = attr(inferred.target.spectra,
                                                  "abundance",
                                                  exact = TRUE))

arrow.tops <- apply(sub.target.sig.ests, 1, max)
arrow.bottoms <- apply(sub.target.sig.ests, 1, min)
min.sub <- min(arrow.bottoms)
bp <- ICAMS::PlotCatalog(inferred.sig, ylim = c(min.sub, max(arrow.tops) + 0.005))

add.arrows <- function(bp, tops, bottoms) {
  oldopt <- getOption("warn")
  on.exit(options(warn = oldopt))
  options(warn = -1) # Does not seem to turn off warnings
  suppressWarnings(
    # Necessary because generates warnings for 0-length arrows
    arrows(
      x0     = bp,
      y0     = tops,    # location of up arrow tips
      y1     = bottoms, # location of down arrow tips
      angle  = 90,      # use "T" arrow heads
      code   = 3,       # use arrow heads at both ends
      length = 0.025    # width of arrow heads
    ))
  # Try this to get rid of the warnings:
  # assign("last.warning", NULL, envir = baseenv())
  
}

add.arrows(bp$plot.object, arrow.tops, arrow.bottoms)


#stop

if (F){# Generate a number of mutations to attribute to the HepG2 background
  #generate.syn.data <- function(cl.num.wanted, bg.catalog, target.catalog){
  
  #}
  set.seed(1)
  # Start with "ideal" spectra construction, use neg binorm, with size param at 10
  num.muts <- 1738
    #rnbinom(1, size = 10, mu = h2$count.nbinom.mu);num.muts
  syn.target.num.muts <- num.muts #rnbinom(1, size = 10, mu = target.mean);syn.target.num.muts
  ideal.bg.spec <- h2$background.sig * num.muts
  ideal.target.spec <- target.sig * syn.target.num.muts
  attr(ideal.bg.spec,'catalog.type') <- 'counts'
  attr(ideal.target.spec, "catalog.type") <- "counts"
  ideal.syn.spec <- ideal.bg.spec + ideal.target.spec
  attr(ideal.syn.spec, "$dimnames[[2]]") <- "ideal.combined"
  colnames(ideal.syn.spec) <- "ideal.combined"
  
  ideal.stacked.spec <- Plot1StackedSpectrum(
    background.spectrum = ideal.bg.spec,
    target.spectrum = ideal.target.spec
  )
  
  #ICAMS::PlotCatalogToPdf(cbind(ideal.bg.spec,ideal.target.spec,ideal.syn.spec),file='ideal_spec.pdf')
  
  #generate noises for partial bg and target spec 
  #following rnbinom
  
  bg.partial.spectrum.w.noise  <- unlist(lapply(ideal.bg.spec,  function(x) rnbinom(1, mu = x, size = 10)))
  target.partial.spectrum.w.noise  <- unlist(lapply(ideal.target.spec,  function(x) rnbinom(1, mu = x, size = 10)))
  syn.spec <- bg.partial.spectrum.w.noise + target.partial.spectrum.w.noise
  partial.stacked.spec <- Plot1StackedSpectrum(
    background.spectrum = bg.partial.spectrum.w.noise,
    target.spectrum = target.partial.spectrum.w.noise
  )
  bg.partial.spectrum.w.noise <- as.matrix(bg.partial.spectrum.w.noise)
  row.names(bg.partial.spectrum.w.noise) <- rownames(ideal.bg.spec)
  bg.partial.spectrum.w.noise <- ICAMS::as.catalog(bg.partial.spectrum.w.noise)
  target.partial.spectrum.w.noise <- as.matrix(target.partial.spectrum.w.noise)
  row.names(target.partial.spectrum.w.noise) <- rownames(ideal.target.spec)
  target.partial.spectrum.w.noise <- ICAMS::as.catalog(target.partial.spectrum.w.noise)
  syn.spec <- as.matrix(syn.spec)
  row.names(syn.spec) <- rownames(ideal.syn.spec)
  syn.spec <- ICAMS::as.catalog(syn.spec)
  colnames(bg.partial.spectrum.w.noise) <- paste0("HepG2.background.w/noise")
  colnames(target.partial.spectrum.w.noise) <- paste0("SBS22.w/noise")
  colnames(syn.spec) <- paste0("combined.w/noise")
  
  ICAMS::PlotCatalogToPdf(cbind(
    bg.partial.spectrum.w.noise,target.partial.spectrum.w.noise,syn.spec),file='noise_spec.pdf')
  ICAMS::PlotCatalogToPdf(cbind(ideal.bg.spec,ideal.target.spec,ideal.syn.spec,
                                bg.partial.spectrum.w.noise,target.partial.spectrum.w.noise,syn.spec),file='combined_spec.pdf')
  
  generate.bg.nbinom <- function(x){
    ff <- rnbinom(1,size=100,mu=ceiling(num.muts/96))
    return(ff)
  }
  
  # Partial spectrum with nbinomial noise
  bg.partial.spect <- apply(h2$background.sig, 1, generate.bg.nbinom);sum(bg.partial.spect)
  
  # 2 Generate partial spectrum for SBS22 (AA)
  
  # remotes::install_github("steverozen/PCAWG7")
  
  # Talk to Nanhai on what distribution to use, we will just use a normal
  # distribution
  target <- "SBS22"
  
  target.sig <- PCAWG7::signature$genome$SBS96[ , target, drop = FALSE]
  target.exp <- PCAWG7::exposure$other.genome$SBS96[target,,drop = FALSE]
  target.exp <- target.exp[target.exp > 0.5]
  target.mean <- mean(target.exp)
  target.sd   <- sd(target.exp)
  target.size <- length(target.exp)
  
  syn.target.num.muts <- rnbinom(1, size = 10, mu = target.mean);syn.target.num.muts
  
  generate.t.nbinom <- function(x){
    ff <- rnbinom(1,size=100,mu=ceiling(syn.target.num.muts/96))
    return(ff)
  }
  target.partial.spect <- apply(target.sig, 1, generate.t.nbinom);sum(target.partial.spect)
  
  syn.spectrum <- bg.partial.spect + target.partial.spect
  a <- h2$background.sig*100
  b<-target.sig *100
  attr(a,'catalog.type') <- 'counts'
  attr(syn.spectrum, "catalog.type") <- "counts"
  attr(b, "catalog.type") <- "counts"
  
  
  Plot1StackedSpectrum(
    background.spectrum = bg.partial.spect,
    target.spectrum = target.partial.spect
  )
  ICAMS::PlotCatalogToPdf(cbind(a,b,syn.spectrum),file='test.pdf')
}

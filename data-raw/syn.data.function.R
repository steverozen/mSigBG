library(PCAWG7)

set.seed(101010, kind = "L'Ecuyer-CMRG")

#target: 'SBSXX' form
#bg: 'mSigBG::background.info$HepG2', spectra form
#ratio: num of target muts/num of muts in bg, must be <= 1, default = 1
#num: number of synthesized cell lines

generate.syn.spec <- function(num, ratio=1, target, bg){
  target.sig <- PCAWG7::signature$genome$SBS96[ , target, drop = FALSE]
  name.of.bg = strsplit(colnames(bg$background.sig), '[.]')[[1]][1]
  
  #intialize empty variables
  placeholder <- target.sig
  attr(placeholder,'catalog.type') <- 'counts'
  plotlist <- placeholder#place holder, dont know how to generate empty catalog
  output.syn.spec <- placeholder
  diff.muts.in.bg <- vector()
  diff.muts.in.target <- vector()
  
  for(i in 1:num){
    num.muts <- rnbinom(1, size = bg$count.nbinom.size, mu = bg$count.nbinom.mu)
    syn.target.num.muts <- round(ratio * num.muts)
    diff.muts.in.bg <- c(diff.muts.in.bg, num.muts)
    diff.muts.in.target <- c(diff.muts.in.target, syn.target.num.muts)
    print(c(num.muts, syn.target.num.muts))
    
    #generate ideal bg and target spectra
    ideal.bg.spec <- bg$background.sig * num.muts
    ideal.target.spec <- target.sig * syn.target.num.muts
    ideal.syn.spec <- ideal.bg.spec + ideal.target.spec
    attr(ideal.bg.spec, "catalog.type") <- "counts"
    attr(ideal.target.spec, "catalog.type") <- "counts"
    attr(ideal.target.spec, "class") <- attr(ideal.bg.spec, "class")
    attr(ideal.syn.spec, "catalog.type") <- "counts"
    attr(ideal.syn.spec, "class") <- attr(ideal.bg.spec, "class")
    colnames(ideal.bg.spec) <- paste0(name.of.bg,'.cl',i)
    colnames(ideal.target.spec) <- paste0(target,'.cl',i)
    colnames(ideal.syn.spec) <- paste0("ideal.combined",'.cl',i)
    plotlist <- cbind(plotlist, 
                      ideal.bg.spec,ideal.target.spec,ideal.syn.spec)
    
    #add nbinorm noise to each column
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
    colnames(bg.partial.spectrum.w.noise) <- paste0(name.of.bg,".w/noise",'.cl',i)
    colnames(target.partial.spectrum.w.noise) <- paste0(target, ".w/noise",'.cl',i)
    colnames(syn.spec) <- paste0("combined.w/noise",'.cl',i)
    attr(bg.partial.spectrum.w.noise, "catalog.type") <- "counts"
    attr(target.partial.spectrum.w.noise, "catalog.type") <- "counts"
    attr(syn.spec,'catalog.type') <- 'counts'
    plotlist <- cbind(plotlist, 
                      bg.partial.spectrum.w.noise, target.partial.spectrum.w.noise, syn.spec)
    output.syn.spec <- cbind(output.syn.spec, syn.spec)
    attr(output.syn.spec,'abundance') <- NULL
    attr(output.syn.spec,'region') <- 'unknown'
    
  }
  ICAMS::PlotCatalogToPdf(plotlist[,-1],file=paste0('specs_for_', num,'cls_for_', ratio, target, '_in_', name.of.bg,'.pdf'))
  
  return(list(output.syn.spec, diff.muts.in.bg, diff.muts.in.target, target.sig))
}

#output: pdf of catalogs for different cell lines;
#output.syn.spec: only contain catalogs of synthetic spectra with noise
#diff.muts.in.bg, diff.muts.in.target: 
  #total number of mutations used to intialize ideal spectra of bg and target
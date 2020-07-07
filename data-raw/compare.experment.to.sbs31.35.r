
inferred.spectra.list <- lapply(all.separated, `[[`, "inferred.target.spectra")
all.target.spectra <- do.call(cbind, inferred.spectra.list)

all.target.sig.list <- lapply(all.separated, `[[`, "inferred.target.sig.as.catalog")

all.target.sigs <- do.call(cbind, all.target.sig.list)

grist <- cbind(
  ICAMS::TransformCatalog(all.target.spectra, 
                          target.catalog.type = "counts.signature"), 
  all.target.sigs)

cosine.sim <- philentropy::distance(t(grist), method = "cosine")
rownames(cosine.sim) <- colnames(grist)
colnames(cosine.sim) <- colnames(grist)
gplots::heatmap.2(x = cosine.sim,
                  dendrogram = "column",
                  margins = c(9, 9),
                  cex.axis = 0.5,
                  symm = TRUE,
                  trace = "none")



pc <- prcomp(t(grist), center = FALSE, scale = FALSE, retx = TRUE)

factoextra::fviz_pca_ind(pc, repel = TRUE)

factoextra::fviz_contrib(
  pc, choice = "var", 
  font.xtickslab = 6, 
  xtickslab.rt = 90)


sbs31 <- PCAWG7::signature$genome$SBS96[ , "SBS31"]
sbs35 <- PCAWG7::signature$genome$SBS96[ , "SBS35"]

one.linear.decomp <- function(sig.name) {
  my.sig <- all.target.sigs[ , sig.name]
  lm <- glm(my.sig ~ sbs31 + sbs35)
  coef <- lm$coefficients[-1]
  recon <- sbs31 * coef["sbs31"] + sbs35 * coef["sbs35"]
  cos.dist <- philentropy::distance(rbind(recon, my.sig), method = "cosine")
  return(list(sig = sig.name, vect = c(cos.dist = cos.dist), coef))
}

foo <- lapply(colnames(all.target.sigs), one.linear.decomp)
lapply(foo, `[[`, c("sig"))



one.opt <- function(sig, method = "cosine", invert = -1) {
  
  sig <- as.vector(sig)
  
  my.obj.fn <- function(coef) {
    recon <- coef[1]*sbs31 + coef[2]*sbs35 
    return(
      suppressMessages(
        invert * philentropy::distance(
          rbind(recon, sig), method = method)))
  }
  
  g_ineq <- function(coef) {
    abs(1 - sum(coef))  
  }
  
  retval <- nloptr::nloptr(
    x0          = c(0.5, 0.5),
    eval_f      = my.obj.fn,
    eval_g_ineq = g_ineq,
    lb          = c(0, 0),
    ub          = c(1, 1),
    opt         = list(
      algorithm = "NLOPT_LN_COBYLA",
      maxeval   = 2000,
      xtol_rel  = 1e-6,
      xtol_abs  = 1e-7)
  )
  
  retval <- list(similarity = invert * retval$objective, 
              coef = retval$solution / sum(retval$solution))
  
  return(retval)
  
}

retval <- one.opt(all.target.sigs[ , 1])

t( 
  as.data.frame(
    lapply(
      all.target.sig.list, 
      function(x) { ret <- one.opt(x)
      names(ret$coef) <- c("sbs31", "sbs35")
      return(c(cos.sim = ret$similarity, ret$coef))})
  ))



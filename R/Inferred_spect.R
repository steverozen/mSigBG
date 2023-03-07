#### Two Ways of getting Inferred spectra 
##### Method 1: obs.spectra - (bg.exposure * the background signature) 
# Inferred spectra obtained from observed spectra minus background spectra 

infer_spec_M1 <- function(obs.spectra, spectra.sep) {
  bg.exposure <- spectra.sep$exposures.to.bg.sig
  bg.spectra <- spectra.sep$inferred.bg.spectra
  bg.sig <- spectra.sep$inferred.bg.spectra/colSums(bg.spectra)
  
  target.spectra <- obs.spectra[,1, drop=FALSE]  - bg.exposure[1] * bg.sig[, 1, drop=FALSE]
  if(length(bg.exposure)>1){
    for (i in 2:length(bg.exposure)){
      target.spectra <- cbind(target.spectra, obs.spectra[,i, drop=FALSE]  - bg.exposure[i] * bg.sig[, i, drop=FALSE])
    }
  }
  #target.sig.method1 <- target.spectra/colSums(target.spectra) # convert clonal target spectra to signature 
  
  return(target.spectra)
  #return(target.sig.method1)
}

###############################################################################################
##### Method 2: (obs.counts - bg.exposure) * sig.to.return 
# Inferred spectra obtained from inferred target signature multiply # of target mutations (observed -bg) 

# No. of mutation due to bg in each sample
infer_spec_M2 <- function(obs.spectra, spectra.sep){
  bg.exposures <- spectra.sep$exposures.to.bg.sig
  names(bg.exposures) <- sub("-inf.bg.spect", "", colnames(spectra.sep$inferred.bg.spectra), fixed = TRUE)
  
  target.sig <-spectra.sep$inferred.target.sig.as.catalog
  
  inferred.target.count <- sum(obs.spectra[,1]) - bg.exposures[1]
  target.spec.method2 <- inferred.target.count * target.sig
  
  if(length(bg.exposures>1)){
    for (j in 2:length(bg.exposures)){
      inferred.target.count <- sum(obs.spectra[,j]) - bg.exposures[j]
      temp.inferred.spectra <- inferred.target.count * target.sig
      target.spec.method2 <- cbind(target.spec.method2,temp.inferred.spectra)
    }
  }
  return(target.spec.method2)
}


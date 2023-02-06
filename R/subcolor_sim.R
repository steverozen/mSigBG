subcolor_sim <- function(spectra1, spectra2){
  
  # A function to obtain cosine similarity between all subcolor of two input spectra/signature
  # The two input spectra should be a one-dimentional vector 
  
  all.subcolor.cossim <- c()
  weight <- c()
  
  for (j in seq(from = 1, to = 96, by = 16)) {
    s1 <- spectra1[j:(j+15), 1, drop = FALSE]
    s2 <- spectra2[j:(j+15), 1, drop = FALSE]
    
    
    all.subcolor.cossim <- c(all.subcolor.cossim,
                             philentropy::cosine_dist(s1, s2, testNA = TRUE))
    weight <- c(weight, sum(s1, s2))
    
  }
  
  names(all.subcolor.cossim) <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  
  return(list(by.color     = all.subcolor.cossim, 
              weighted.avg = weighted.mean(all.subcolor.cossim, weight)))
}

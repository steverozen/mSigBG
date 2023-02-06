subcolor_sim <- function(spectra1, spectra2){
  
  # A function to obtain cosine similarity between all subcolor of two input spectra/signature
  # The two input spectra should be a one-dimentional vector 
  
  all.subcolor.cossim <- c()
  
  for (j in seq(from = 1, to = 96, by = 16)){
    this.cossim <- 
      philentropy::cosine_dist(spectra1[j:(j+15), 1, drop = FALSE], 
                               spectra2[j:(j+15), 1, drop = FALSE],
                               testNA = TRUE)
    

    all.subcolor.cossim <- c(all.subcolor.cossim, this.cossim)
  }
  
  names(all.subcolor.cossim) <- c("C>A","C>G","C>T","T>A","T>C","T>G")
  
  return(all.subcolor.cossim)
}

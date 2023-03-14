convert_base <- function(x){
  if (x == "A")
    ret <- "T"
  else if (x == "C")
    ret <- "G"
  else if (x == "G")
    ret <- "C"
  else if (x == "T")
    ret <- "A"
  else{
    message("Incorrent input")
    ret <- 0
  }
  
  return(ret)
}


Jill_to_SBS96 <- function(input_table){
  subs <- input_table[,c("pre_context","Ref","rear_context","Alt")]
  for (i in 1:nrow(input_table)){
    if (subs[i,2] == "A" | subs[i,2] =="G"){
      tmp1 <- convert_base(subs[i,1])
      tmp2 <- convert_base(subs[i,2])
      tmp3 <- convert_base(subs[i,3])
      tmp4 <- convert_base(subs[i,4])
      subs[i,1] <- tmp3
      subs[i,2] <- tmp2
      subs[i,3] <- tmp1
      subs[i,4] <- tmp4
    }
  }
  subs.pasted <- apply(subs, MARGIN = 1, function(x){paste0(x[1], x[2], x[3], x[4])})
  SBS96.template <- c("ACAA", "ACCA", "ACGA", "ACTA", "CCAA", "CCCA", "CCGA", "CCTA", 
                      "GCAA", "GCCA", "GCGA", "GCTA", "TCAA", "TCCA", "TCGA", "TCTA", 
                      "ACAG", "ACCG", "ACGG", "ACTG", "CCAG", "CCCG", "CCGG", "CCTG", 
                      "GCAG", "GCCG", "GCGG", "GCTG", "TCAG", "TCCG", "TCGG", "TCTG",
                      "ACAT", "ACCT", "ACGT", "ACTT", "CCAT", "CCCT", "CCGT", "CCTT", 
                      "GCAT", "GCCT", "GCGT", "GCTT", "TCAT", "TCCT", "TCGT", "TCTT", 
                      "ATAA", "ATCA", "ATGA", "ATTA", "CTAA", "CTCA", "CTGA", "CTTA", 
                      "GTAA", "GTCA", "GTGA", "GTTA", "TTAA", "TTCA", "TTGA", "TTTA", 
                      "ATAC", "ATCC", "ATGC", "ATTC", "CTAC", "CTCC", "CTGC", "CTTC", 
                      "GTAC", "GTCC", "GTGC", "GTTC", "TTAC", "TTCC", "TTGC", "TTTC", 
                      "ATAG", "ATCG", "ATGG", "ATTG", "CTAG", "CTCG", "CTGG", "CTTG", 
                      "GTAG", "GTCG", "GTGG", "GTTG", "TTAG", "TTCG", "TTGG", "TTTG") 
  subs.spectra <- rep(0,96)
  names(subs.spectra) <- SBS96.template
  
  # count occurrence 
  for (i in subs.pasted){
    subs.spectra[i] <- subs.spectra[i] + 1
  }
  
  subs.cat <- as.catalog(as.matrix(subs.spectra), ref.genome = "GRCh37", region = "genome", catalog.type = "counts" , infer.rownames = TRUE)
  
  return(subs.cat)
}
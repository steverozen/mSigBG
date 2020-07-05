

# Run this to create package variables HepG2.background.spectra
# and MCF10A.background.spectra
HepG2.background.spectra <- ICAMS::ReadCatalog(
  "data-raw/package.variable.source.data/Background_HepG2.csv",
  region = "genome",
  catalog.type = "counts")
attr(HepG2.background.spectra, "abundance") <- 
  ICAMS::all.abundance$BSgenome.Hsapiens.1000genomes.hs37d5$genome$"96"
usethis::use_data(HepG2.background.spectra, overwrite = TRUE)


MCF10A.background.spectra <- ICAMS::ReadCatalog(
  "data-raw/package.variable.source.data/Background_MCF10A.csv",
  region = "genome",
  catalog.type = "counts")
attr(MCF10A.background.spectra, "abundance") <- 
  ICAMS::all.abundance$BSgenome.Hsapiens.1000genomes.hs37d5$genome$"96"
usethis::use_data(MCF10A.background.spectra, overwrite = TRUE)
                          

# Run this code to create the package variables HepG2.background
# and MCF10A.background.
background.info <- list()
background.info[["HepG2"]] <- 
  MakeBackgroundInfo(mSigBG::HepG2.background.spectra, title = "HepG2.background")
background.info[["MCF10A"]] <- 
  MakeBackgroundInfo(mSigBG::MCF10A.background.spectra, title = "MCF10A.background")
usethis::use_data(background.info, overwrite = TRUE)

vcfs <- list.files("data-raw/example.vcfs/HepG2_Cis/", full.names = TRUE)
hepg2.cis <- ICAMS::StrelkaSBSVCFFilesToCatalog(
  files      = vcfs,
  ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
  region     = "genome")$catSBS96
attr(hepg2.cis, "ref.genome") <- NULL

vcfs <- list.files("example.vcfs/MCF10A_Cis/", full.names = TRUE)
mcf10a.cis <- ICAMS::StrelkaSBSVCFFilesToCatalog(
  files      = vcfs,
  ref.genome = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
  region     = "genome")$catSBS96
attr(mcf10a.cis, "ref.genome") <- NULL
example.spectra <- list()
example.spectra[["HepG2.cisplatin"]] <- hepg2.cis
example.spectra[["MCF10A.cisplatin"]] <- mcf10a.cis
usethis::use_data(example.spectra, overwrite = TRUE)


OLD.MakeMCF10HepG2BackgroundVars <- function() {
  data.dir <- file.path("data-raw",
                        "background.signature.spectra", 
                        "backgroundVCFS-2020-02-04")
  sbs.files <- dir(path = data.dir, pattern = "SNV", full.names = TRUE)
  names <- sub(".*((MCF10A|HepG2)_SC._cl.).*", "\\1", sbs.files, perl = TRUE)
  
  # No substantial indel background, so will not read these VCFs
  # id.files  <- dir(path = data.dir, pattern = "INDEL", full.names = TRUE)
  sbs.cat <- 
    ICAMS::StrelkaSBSVCFFilesToCatalogAndPlotToPdf(
      files         = sbs.files, 
      ref.genome    = BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5,
      region        = "genome",
      names.of.VCFs = names,
      output.file = file.path(data.dir, "MCF10A-and-HepG2-bg"))
  ICAMS::WriteCatalog(sbs.cat$catSBS96, 
                      file.path(data.dir, "MCF10A-and-HepG2-bg-SBS96-bg.csv"))
  ICAMS::WriteCatalog(sbs.cat$catSBS192, 
                      file.path(data.dir, "MCF10A-and-HepG2-bg-SBS192-bg.csv"))
  
  # The ref genomes have an arbitrary file path encoded in them
  attr(sbs.cat$catSBS96, "ref.genome") <- NULL
  
  HepG2.background.spectra  <- sbs.cat$catSBS96[ , 1:3]
  MCF10A.background.spectra <- sbs.cat$catSBS96[ , 4:6]
  usethis::use_data(HepG2.background.spectra, overwrite = TRUE)
  usethis::use_data(MCF10A.background.spectra, overwrite = TRUE)
}

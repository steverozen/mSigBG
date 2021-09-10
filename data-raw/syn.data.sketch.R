# Example to generate one synthetic spectrum due to HepG2 background and
# SBS22 (AA, easy) SBS4 (smoking, hard example)

# 1. Generate partial spectrum for background

# We use info in mSigBG::background.info

h2 <- mSigBG::background.info$HepG2

# Generate a number of mutations to attribute to the HepG2 backgroun

num.muts <- 
  rnbinom(1, size = h2$count.nbinom.size, mu = h2$count.nbinom.mu)


# Partial spectrum with multinomial noise

pg.partial.spect <- rmultinom(1, num.muts, h2$background.sig)

# 2 Generate partial spectrum for SBS22 (AA)

library(PCAWG7)
# remotes::install_github("steverozen/PCAWG7")

# Talk to Nanhai on what distribution to use, we will just use a normal
# distribution
target <- "SBS22"

target.sig <- PCAWG7::signature$genome$SBS96[ , target, drop = FALSE]

target.exp <- PCAWG7::exposure$other.genome$SBS96[target, ]
target.exp <- target.exp[target.exp > 0.5]
target.mean <- mean(target.exp)
target.sd   <- sd(target.exp)

syn.target.num.muts <- 2000 # rnorm(1, mean = target.mean, sd = target.sd)

target.partial.spect <- rmultinom(1, syn.target.num.muts, target.sig)

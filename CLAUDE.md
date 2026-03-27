# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Package Overview

mSigBG (Mutational Signature BackGround) separates a background mutational signature from observed mutational spectra. Designed for delineating signatures from cell cultures exposed to mutagens.

The core approach uses maximum likelihood estimation with negative binomial distributions to decompose each observed spectrum into a background component (from untreated cells) and a target signature (from mutagen exposure). Optimization uses `nloptr` with the COBYLA algorithm.

## Build and Test Commands

```bash
# Run all tests
Rscript -e 'devtools::test()'

# Run a single test file
Rscript -e 'testthat::test_file("tests/testthat/test-SeparateSignatureFromBackground.R")'

# R CMD check
Rscript -e 'devtools::check()'

# Build and install
R CMD INSTALL .

# Regenerate documentation from roxygen2 comments
Rscript -e 'devtools::document()'
```

## Architecture

### Core Workflow

1. **`MakeBackgroundInfo()`** — Creates a background signature info structure from untreated cell spectra. Returns a list with `background.sig`, negative binomial parameters (`sig.nbinom.size`, `count.nbinom.mu`, `count.nbinom.size`), and `input.spectra`.

2. **`SeparateSignatureFromBackground()`** — The main optimization function. Takes exposed-cell spectra and background info, estimates the target signature profile and background mutation counts per spectrum. Uses constrained optimization (signature probabilities sum to 1, background exposure ≤ total count).

3. **`SeparateSignatureAndSpectra()`** — Wrapper that calls `SeparateSignatureFromBackground()` and also returns the decomposed spectra (background + target components).

### Supporting Functions

- **`LLHSpectrumNegBinom()`** — Calculates log-likelihood of an observed spectrum under the negative binomial model.
- **`MeanOfSpectraAsSig()`** — Converts spectra to signatures and computes their mean.
- **`ObjFn1()` / `NegLLHOfSignature()`** — Internal objective and likelihood functions for the optimizer.

### Visualization

- **`Plot1StackedSpectrum()`** — Stacked bar chart showing background vs target components.
- **`PlotSpectraAsSigsWithUncertainty()`** — Mean signature with uncertainty bars.
- **`plot_stacked_sigs_by_exposure()`** — Stacked signatures scaled by exposure.

## Key Dependencies

- **ICAMS** — Handles mutational signature catalogs (SBS96 format). All spectra and signatures are ICAMS catalog objects.
- **nloptr** — Constrained nonlinear optimization.

## Package Data

- `background.info` — Pre-computed background info for HepG2 and MCF10A cell lines.
- `example.spectra` — Example cisplatin-exposed spectra (MCF10A.cisplatin).
- `MCF10A.background.spectra` / `HepG2.background.spectra` — Raw background spectra catalogs.

Data regeneration script: `data-raw/PrepareBackgroundSignaturePackageData.R`.

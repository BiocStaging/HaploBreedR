# HaploBreedR <img src="man/figures/logo.png" align="right" height="139" alt="HaploBreedR logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/aksruchi/HaploBreedR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/aksruchi/HaploBreedR/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/aksruchi/HaploBreedR/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/aksruchi/HaploBreedR/actions/workflows/pkgdown.yaml)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**Superior haplotype discovery and haplotype-guided genomic selection for crop breeding.**

HaploBreedR is an end-to-end quantitative genomics and molecular
breeding suite. It identifies superior functional haplotypes from
nonsynonymous CDS and UTR variants under false discovery rate control,
fits two-kernel spectral REML genomic prediction models with the
pipeline's superior haplotypes as a feature kernel, performs
haplotype-guided cross-usefulness prediction and linkage-aware
pyramiding optimization, and selects diverse parent sets with a
relatedness ceiling.

📊 **[One-page infographic](https://aksruchi.github.io/HaploBreedR/infographic.html)** — the pipeline at a glance.
📊 **[One-page infographic](https://aksruchi.github.io/HaploBreedR/infographic.html)** ([print TIFF](https://aksruchi.github.io/HaploBreedR/reference/figures/HaploBreedR_infographic.tiff)) — the pipeline at a glance.
---

## Installation

```r
# install.packages("devtools")
devtools::install_github("aksruchi/HaploBreedR")
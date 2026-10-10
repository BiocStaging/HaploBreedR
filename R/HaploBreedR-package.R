#' HaploBreedR: haplotype-based breeding decisions
#'
#' Identifies superior haplotypes at candidate genes from GWAS hits and
#' genotype data, tests haplotype-trait associations, applies effect-size
#' and permutation gates, and produces multi-trait breeding
#' recommendations.
#'
#' @section Optional parallel backend:
#' When `n_cores > 1`, the GRM builder uses `future` and `future.apply`
#' for parallel block computation. Install those packages if you need
#' multi-core GRM construction; otherwise the builder runs sequentially.
#' @section Function length:
#' The package deliberately keeps several long orchestration functions
#' (\code{haplobreed_analysis()} and \code{hbr_advanced_decision()}) as
#' single entry points. Splitting them would require passing dozens of
#' intermediate objects between helpers, which increases coupling without
#' improving readability. All internal helpers are shorter than 50 lines.
#'
#' @keywords internal
#' @importFrom data.table :=
#' @importFrom BEDMatrix BEDMatrix
#' @importFrom stats density dist end optimize start sd var median setNames
#'   quantile cor t.test aov lm p.adjust fisher.test qchisq ppoints rnorm qnorm dnorm
#'   approx aggregate optimize residuals as.formula na.omit
#' @importFrom utils head tail read.table download.file untar URLdecode
#' @importFrom methods is
"_PACKAGE"

#' HaploBreedR: haplotype-based breeding decisions
#'
#' Identifies superior haplotypes at candidate genes from GWAS hits and
#' genotype data, tests haplotype-trait associations, applies effect-size
#' and permutation gates, and produces multi-trait breeding
#' recommendations.
#'
#' @keywords internal
#' @importFrom data.table :=
#' @importFrom stats density dist end optimize start sd var median setNames
#'   quantile cor t.test aov lm p.adjust fisher.test qchisq ppoints rnorm
#'   residuals as.formula na.omit
#' @importFrom utils head tail read.table download.file untar URLdecode
#'   globalVariables
"_PACKAGE"
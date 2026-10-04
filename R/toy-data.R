#' Tiny demonstration dataset for HaploBreedR
#'
#' A minimal toy dataset shipped with the package so that examples,
#' tests, and the vignette can run without network access. It is not
#' biologically meaningful; it is sized only to exercise the code paths.
#'
#' Files live in `system.file("extdata", "toy", package = "HaploBreedR")`:
#'
#' \describe{
#'   \item{toy.raw, toy.bim, toy.fam}{PLINK binary input for 10
#'     individuals and 50 SNPs on one chromosome.}
#'   \item{toy_pheno.csv}{Two columns, `ID` and `trait` (a continuous
#'     response).}
#'   \item{toy_gwas.csv}{Three GWAS hits on chromosome 1.}
#'   \item{toy.gff3}{A three-kilobase mock gene model with one
#'     transcript, CDS, and 5'-/3'-UTRs.}
#'   \item{toy.fa}{A placeholder genome FASTA to satisfy the reference
#'     validator.}
#' }
#'
#' @format Files under `inst/extdata/toy/`. See `inst/scripts/make_toy_data.R`
#'   for the generation script.
#' @source Simulated.
#' @return `NULL`, invisibly. The object is documentation-only; the files
#'   are read with `system.file("extdata", "toy", package = "HaploBreedR")`.
#' @keywords internal
#' @name toy_data
NULL
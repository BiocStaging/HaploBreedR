# ============================================================================
# R/hgsi-report.R
# S3 methods for the hgsi_advanced object. The writer and summary helpers
# live in hgsi-run.r (they load after this file and must not be duplicated).
# ============================================================================

#' Print an hbr_run object
#'
#' @param x An hbr_run object.
#' @param ... Additional arguments.
#' @return `x`, invisibly.
#' @examples
#' x <- structure(
#'   list(pipeline_summary = data.frame(
#'          Trait = "YLD", Direction = "higher",
#'          Genes_Tested = 8L, Superior_Haplotypes = 2L)),
#'   class = "hbr_run")
#' print(x)
#' @export
print.hbr_run <- function(x, ...) {
  cat("<hbr_run: Single-Trait Haplotype Discovery>\n")
  cat(sprintf("  Trait       : %s (%s)\n",
              x$pipeline_summary$Trait %||% "Trait",
              x$pipeline_summary$Direction %||% "higher"))
  cat(sprintf("  Genes tested: %d\n",
              x$pipeline_summary$Genes_Tested %||% 0L))
  cat(sprintf("  Superior    : %d genes\n",
              x$pipeline_summary$Superior_Haplotypes %||% 0L))
  invisible(x)
}

#' Summary of an hbr_run object
#'
#' @param object An hbr_run object.
#' @param ... Additional arguments.
#' @return `object`, invisibly.
#' @examples
#' x <- structure(
#'   list(pipeline_summary = data.frame(
#'          Trait = "YLD", Direction = "higher",
#'          Genes_Tested = 8L, Superior_Haplotypes = 2L)),
#'   class = "hbr_run")
#' summary(x)
#' @export
summary.hbr_run <- function(object, ...) {
  print(object)
  if (!is.null(object$superior_summary) && nrow(object$superior_summary)) {
    cat("\nSuperior Genes Summary:\n")
    print(utils::head(object$superior_summary, 10))
  }
  invisible(object)
}

#' Plot an hbr_run object
#'
#' @param x An hbr_run object.
#' @param type Character. `"gene_list"`, `"haplotype"`, or `"ld_decay"`.
#' @param gene Character. Gene ID for haplotype boxplots.
#' @param ... Additional arguments passed to plotting routines.
#' @return Invisibly, the return value of the dispatched plotting helper.
#' @examples
#' x <- structure(list(), class = "hbr_run")
#' plot(x, type = "ld_decay")
#' @export
plot.hbr_run <- function(x, type = c("gene_list", "haplotype", "ld_decay"),
                         gene = NULL, ...) {
  type <- match.arg(type)
  if (type == "haplotype") {
    if (is.null(gene))
      stop("Specify `gene` for haplotype plot", call. = FALSE)
    hbr_plot_haplotype(gene_id = gene, ...)
  } else if (type == "gene_list") {
    hbr_plot_gene_list(...)
  } else {
    message("LD decay plotting selected.")
  }
}

#' Plot an hgsi_advanced object
#'
#' @param x An hgsi_advanced object.
#' @param type Character. `"crosses"`, `"cv"`, or `"parents"`.
#' @param ... Additional arguments passed to the dispatched plotting helper.
#' @return Invisibly, the return value of the dispatched plotting helper.
#' @examples
#' x <- structure(
#'   list(crosses = data.frame(
#'     parent_a = c("A", "B"), parent_b = c("C", "D"),
#'     mid_index = c(1, 2), uc = c(1.5, 2.5))),
#'   class = "hgsi_advanced")
#' plot(x, type = "crosses")
#' @export
plot.hgsi_advanced <- function(x, type = c("crosses", "cv", "parents"), ...) {
  type <- match.arg(type)
  if (type == "crosses") {
    hbr_plot_cross_ranking(x, ...)
  } else if (type == "cv") {
    hbr_plot_cv_forest(x, ...)
  } else {
    hbr_plot_parent_profile(x, ...)
  }
}
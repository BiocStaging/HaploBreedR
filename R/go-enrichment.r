# ============================================================================
# R/go-enrichment.r
# GO term over-representation and cross-system gene annotation.
#
# ============================================================================
# ----------------------------------------------------------------------------
# GO term enrichment
# ----------------------------------------------------------------------------

#' GO term enrichment across a set of genes
#'
#' One-sided Fisher's exact test (over-representation) per GO term, with
#' multiple-testing adjustment over the terms tested in the aspect.
#' Aspects BP (biological process), MF (molecular function), and CC
#' (cellular component) are recognised under all common naming
#' conventions (`BP`/`P`/`biological_process`, etc.). Note: candidate
#' genes from GWAS intervals are not independent (LD clusters,
#' overlapping genes), so treat the p-values as exploratory.
#'
#' @param genes_superior Character vector of genes of interest.
#' @param genes_all Character vector of the full tested universe.
#' @param goslim_dt GOSlim table with columns `Gene_ID`, `GO_ID`,
#'   `GO_Aspect`, `GO_Term`.
#' @param aspect One of `"BP"`, `"MF"`, `"CC"`, `"all"`.
#' @param min_term_size Minimum term size in the universe.
#' @param min_superior Minimum superior carriers.
#' @param p_adjust Multiple-testing method.
#' @param top_n Return at most this many top-ranked terms.
#' @param restrict_to_annotated Restrict the universe to genes with at
#'   least one GO term in the aspect.
#' @param verbose Print progress.
#' @examples
#' goslim <- data.table::data.table(
#'   Gene_ID   = c("Os01g0100100", "Os01g0100200", "Os01g0100300"),
#'   GO_ID     = c("GO:0008150", "GO:0008150", "GO:0003674"),
#'   GO_Aspect = c("P", "P", "F"),
#'   GO_Term   = c("biological_process", "biological_process",
#'                 "molecular_function"))
#' hbr_go_enrichment(
#'   genes_superior = c("Os01g0100100", "Os01g0100200"),
#'   genes_all      = c("Os01g0100100", "Os01g0100200", "Os01g0100300"),
#'   goslim_dt      = goslim,
#'   aspect         = "BP",
#'   min_term_size  = 1L,
#'   min_superior   = 1L)
#' @return A `data.table` with columns `GO_ID`, `GO_Term`, `GO_Aspect`,
#'   `N_Superior`, `N_Universe`, `OR`, `P`, `Adjusted_P`.
#'
#' @export
hbr_go_enrichment <- function(genes_superior,
                               genes_all,
                               goslim_dt,
                               aspect        = c("BP", "MF", "CC", "all"),
                               min_term_size = 5L,
                               min_superior  = 2L,
                               p_adjust      = "BH",
                               top_n         = 50L,
                               restrict_to_annotated = FALSE,
                               verbose       = TRUE) {

  aspect <- match.arg(aspect)

  if (is.null(goslim_dt) || !nrow(goslim_dt))
    stop("`goslim_dt` is empty.", call. = FALSE)

  needed <- c("Gene_ID", "GO_ID", "GO_Aspect", "GO_Term")
  miss <- setdiff(needed, names(goslim_dt))
  if (length(miss))
    stop("`goslim_dt` is missing column(s): ",
         paste(miss, collapse = ", "), call. = FALSE)

  norm <- function(x) hbr_normalize_osa1_gene_id(x)

  asp_codes <- list(
    BP = c("BP", "P", "BIOLOGICAL_PROCESS", "BIOLOGICAL PROCESS"),
    MF = c("MF", "F", "MOLECULAR_FUNCTION", "MOLECULAR FUNCTION"),
    CC = c("CC", "C", "CELLULAR_COMPONENT", "CELLULAR COMPONENT"))
  canon <- function(x) {
    u <- toupper(trimws(as.character(x)))
    ifelse(u %in% asp_codes$BP, "BP",
           ifelse(u %in% asp_codes$MF, "MF",
                  ifelse(u %in% asp_codes$CC, "CC", u)))
  }

  in_set   <- unique(norm(genes_superior))
  universe <- unique(norm(genes_all))
  in_set   <- in_set[in_set %in% universe]

  if (!length(in_set))
    stop("None of `genes_superior` are in `genes_all`.", call. = FALSE)

  gos <- data.table::as.data.table(data.table::copy(goslim_dt))
  gos[, Gene_ID   := norm(Gene_ID)]
  gos[, GO_Aspect := canon(GO_Aspect)]

  if (aspect != "all") {
    present <- sort(unique(gos$GO_Aspect))
    gos <- gos[GO_Aspect == aspect]
    if (!nrow(gos))
      stop(sprintf("No GOSlim rows for aspect '%s'. Aspect codes in the table: %s.",
                   aspect, paste(present, collapse = ", ")), call. = FALSE)
  }

  gos <- unique(gos[, .(Gene_ID, GO_ID, GO_Term, GO_Aspect)])
  gos <- gos[Gene_ID %in% universe]
  if (!nrow(gos))
    stop(sprintf("No universe gene has a GO term in aspect '%s'.",
                 aspect), call. = FALSE)

  if (isTRUE(restrict_to_annotated)) {
    universe <- universe[universe %in% unique(gos$Gene_ID)]
    in_set   <- in_set[in_set %in% universe]
    if (!length(in_set))
      stop("No superior gene has a GO term in aspect '", aspect, "'.",
           call. = FALSE)
  }

  if (verbose)
    cat(sprintf("  Enrichment: %d superior / %d universe genes | aspect = %s%s\n",
                length(in_set), length(universe), aspect,
                if (isTRUE(restrict_to_annotated)) " (annotated genes only)" else ""))

  gos <- unique(gos, by = c("Gene_ID", "GO_ID"))

  n_sup   <- length(in_set)
  n_total <- length(universe)
  n_other <- n_total - n_sup
  others  <- setdiff(universe, in_set)

  terms <- unique(gos[, .(GO_ID, GO_Term, GO_Aspect)], by = "GO_ID")

  results <- vector("list", nrow(terms))
  for (i in seq_len(nrow(terms))) {
    go  <- terms$GO_ID[i]
    has <- unique(gos[GO_ID == go]$Gene_ID)

    a       <- sum(in_set %in% has)
    b       <- n_sup - a
    n_in_o  <- sum(has %in% others)
    d       <- n_other - n_in_o

    if ((a + n_in_o) < min_term_size || a < min_superior) next

    ft <- tryCatch(
      stats::fisher.test(matrix(c(a, b, n_in_o, d), nrow = 2, byrow = TRUE),
                          alternative = "greater"),
      error = function(e) NULL)
    if (is.null(ft)) next

    results[[i]] <- data.table::data.table(
      GO_ID      = go,
      GO_Term    = terms$GO_Term[i],
      GO_Aspect  = terms$GO_Aspect[i],
      N_Superior = a,
      N_Universe = a + n_in_o,
      OR         = round(as.numeric(ft$estimate), 3),
      P          = ft$p.value)
  }

  out <- data.table::rbindlist(
    results[!vapply(results, is.null, logical(1))], fill = TRUE)
  if (!nrow(out)) {
    if (verbose) cat("    no testable GO terms\n")
    return(data.table::data.table())
  }

  out[, Adjusted_P := stats::p.adjust(P, method = p_adjust)]
  data.table::setorder(out, P, -OR)
  out <- out[seq_len(min(top_n, .N))]

  if (verbose) {
    n_sig <- sum(out$Adjusted_P < 0.05, na.rm = TRUE)
    cat(sprintf("    %d terms tested | %d significant (FDR < 0.05)\n",
                nrow(out), n_sig))
  }

  out[]
}

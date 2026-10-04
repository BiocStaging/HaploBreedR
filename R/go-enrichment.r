# ============================================================================
# R/go-enrichment.r
# GO term over-representation and cross-system gene annotation.
#
# ============================================================================


# ----------------------------------------------------------------------------
# Cross-annotation (MSU <-> RAP-DB) with GO terms
# ----------------------------------------------------------------------------

#' Attach the paired gene ID and GO terms from the opposite rice annotation
#'
#' Uses `org.Osativa.eg.db` to map each gene ID in `genes[[id_col]]` to
#' the corresponding ID in the opposite annotation system (MSU to RAP-DB
#' or vice versa) and optionally attaches GO IDs with BP / MF / CC
#' splits. Returns `genes` with the newly attached columns appended.
#'
#' @param genes A `data.table` containing a gene-ID column.
#' @param id_col Name of the gene-ID column in `genes`.
#' @param from_type One of `"MSU"`, `"RAP"`.
#' @param to_type One of `"RAP"`, `"MSU"`.
#' @param prefix Column-name prefix for attached columns.
#' @param attach_go Attach GO ID lists and BP/MF/CC split columns.
#' @param verbose Print progress.
#'
#' @return `genes` with the attached columns.
#'
#' @noRd
hbr_add_cross_annotation <- function(
    genes,
    id_col    = "Gene_ID",
    from_type = c("MSU", "RAP"),
    to_type   = c("RAP", "MSU"),
    prefix    = NULL,
    attach_go = TRUE,
    verbose   = TRUE) {

  from_type <- match.arg(from_type)
  to_type   <- match.arg(to_type)
  if (is.null(prefix))
    prefix <- if (to_type == "RAP") "RAPDB" else "MSU"

  genes <- data.table::as.data.table(data.table::copy(genes))
  if (!id_col %in% names(genes))
    stop("hbr_add_cross_annotation(): '", id_col,
         "' not found in `genes`.", call. = FALSE)

  if (!requireNamespace("org.Osativa.eg.db", quietly = TRUE))
    stop("Package `org.Osativa.eg.db` is required. Install with:\n",
         "  remotes::install_github(\"orgz/org.Osativa.eg.db\")",
         call. = FALSE)
  if (!requireNamespace("AnnotationDbi", quietly = TRUE))
    stop("Package `AnnotationDbi` is required.", call. = FALSE)

  org <- get("org.Osativa.eg.db",
             envir = asNamespace("org.Osativa.eg.db"))

  keytype_from <- if (from_type == "MSU") "GID" else "RAP"
  keytype_to   <- if (to_type   == "MSU") "GID" else "RAP"

  normalize_id <- function(x) {
    x <- trimws(as.character(x))
    x <- sub("^gene:", "", x, ignore.case = TRUE)
    x <- sub("\\..*$", "", x)
    x <- sub("\\|.*$", "", x)
    x <- sub("-[0-9]+$", "", x)
    x
  }

  ids_raw    <- as.character(genes[[id_col]])
  ids_norm   <- normalize_id(ids_raw)
  ids_unique <- unique(ids_norm[!is.na(ids_norm) & nzchar(ids_norm)])

  if (verbose)
    cat(sprintf("    org.Osativa.eg.db: %d unique %s ID(s) for lookup\n",
                length(ids_unique), from_type))

  map <- tryCatch(
    AnnotationDbi::select(org,
                          keys    = ids_unique,
                          columns = keytype_to,
                          keytype = keytype_from),
    error = function(e) {
      warning("org.Osativa.eg.db select() failed for ", keytype_to,
              ": ", conditionMessage(e), call. = FALSE)
      NULL
    })

  out_id <- paste0(prefix, "_Gene_ID")
  genes[, (out_id) := NA_character_]

  if (!is.null(map) && nrow(map)) {
    map_dt <- data.table::as.data.table(map)
    map_dt <- map_dt[!is.na(get(keytype_from)) & get(keytype_from) != ""]
    map_dt <- unique(map_dt, by = keytype_from)

    idx <- match(ids_norm, map_dt[[keytype_from]])
    genes[, (out_id) := map_dt[[keytype_to]][idx]]
  }

  out_go    <- paste0(prefix, "_GO_IDs")
  out_go_bp <- paste0(prefix, "_GO_BP")
  out_go_mf <- paste0(prefix, "_GO_MF")
  out_go_cc <- paste0(prefix, "_GO_CC")
  genes[, (out_go)    := NA_character_]
  genes[, (out_go_bp) := NA_character_]
  genes[, (out_go_mf) := NA_character_]
  genes[, (out_go_cc) := NA_character_]

  if (isTRUE(attach_go)) {
    map_go <- tryCatch(
      AnnotationDbi::select(org,
                            keys    = ids_unique,
                            columns = c("GO", "ONTOLOGY"),
                            keytype = keytype_from),
      error = function(e) {
        warning("GO lookup failed: ", conditionMessage(e), call. = FALSE)
        NULL
      })

    if (!is.null(map_go) && nrow(map_go) && "GO" %in% names(map_go)) {
      go_dt <- data.table::as.data.table(map_go)
      go_dt <- go_dt[!is.na(get(keytype_from)) & !is.na(GO) & GO != ""]

      go_all <- go_dt[, .(GO_IDs = paste(unique(GO), collapse = ";")),
                      by = keytype_from]

      go_split <- function(onto_label) {
        sub_dt <- go_dt[ONTOLOGY == onto_label]
        if (!nrow(sub_dt))
          return(data.table::data.table(k = character(0),
                                        v = character(0)))
        out <- sub_dt[, .(v = paste(unique(GO), collapse = ";")),
                      by = keytype_from]
        data.table::setnames(out, keytype_from, "k")
        out
      }
      bp <- go_split("BP")
      mf <- go_split("MF")
      cc <- go_split("CC")

      idx_all <- match(ids_norm, go_all[[keytype_from]])
      genes[, (out_go) := go_all$GO_IDs[idx_all]]

      if (nrow(bp)) {
        idx_bp <- match(ids_norm, bp$k)
        genes[, (out_go_bp) := bp$v[idx_bp]]
      }
      if (nrow(mf)) {
        idx_mf <- match(ids_norm, mf$k)
        genes[, (out_go_mf) := mf$v[idx_mf]]
      }
      if (nrow(cc)) {
        idx_cc <- match(ids_norm, cc$k)
        genes[, (out_go_cc) := cc$v[idx_cc]]
      }
    }
  }

  if (verbose) {
    n_id <- sum(!is.na(genes[[out_id]]))
    n_go <- sum(!is.na(genes[[out_go]]))
    cat(sprintf(
      "    %s: %d/%d IDs mapped, %d/%d with GO terms\n",
      prefix, n_id, nrow(genes), n_go, nrow(genes)))
  }

  genes[]
}


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

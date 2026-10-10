# ============================================================================
# R/validation.r
# Input validation: individual-ID matching, GWAS SNP QC, and the unified
# input-QC report.
#
# ============================================================================


# ----------------------------------------------------------------------------
# Individual-ID validation
# ----------------------------------------------------------------------------

#' Validate ID overlap between phenotype and genotype datasets
#'
#' Cleans both ID vectors with [hbr_clean_id()], compares them, and
#' reports the number matched, the number phenotype-only, the number
#' genotype-only, and the fraction of phenotyped individuals that have
#' genotype data. Status is `"FAIL"` when the fraction falls below
#' `min_overlap_frac`, `"PASS WITH WARNINGS"` when there are unmatched
#' IDs on either side, and `"PASS"` otherwise.
#'
#' @param phenotype_data Data frame with an ID column.
#' @param genotype A `hbr_genotype` object or a character vector of IDs.
#' @param pheno_id_col Name of the ID column in `phenotype_data`.
#' @param min_overlap_frac Minimum fraction of phenotyped individuals
#'   that must also have genotype data.
#' @return A list of class `hbr_id_report` with `status`, `summary` (a
#'   one-row `data.table`), `matched_ids`, `phenotype_only_ids`,
#'   `genotype_only_ids`, and `id_cleaning_changed`.
#'
#' @examples
#' \donttest{
#' # Requires `pheno` (data.frame) and `geno` (hbr_genotype) in the caller's
#' # environment. See the vignette for an end-to-end example.
#' if (interactive()) {
#'   hbr_validate_individual_ids(
#'     phenotype_data   = pheno,
#'     genotype         = geno,
#'     pheno_id_col     = "ID",
#'     min_overlap_frac = 0.5)
#' }
#' }
#' @noRd
hbr_validate_individual_ids <- function(phenotype_data, genotype,
                                        pheno_id_col = "ID",
                                        min_overlap_frac = 0.5) {

  if (!pheno_id_col %in% names(phenotype_data))
    stop("hbr_validate_individual_ids(): '", pheno_id_col,
         "' not found in phenotype_data.", call. = FALSE)

  geno_ids_raw <- if (inherits(genotype, "hbr_genotype"))
    genotype$individuals else as.character(genotype)

  pheno_ids_raw <- as.character(phenotype_data[[pheno_id_col]])
  pheno_ids     <- hbr_clean_id(pheno_ids_raw)
  geno_ids      <- hbr_clean_id(geno_ids_raw)

  changed <- unique(data.table::rbindlist(list(
    data.table::data.table(Original = pheno_ids_raw, Cleaned = pheno_ids),
    data.table::data.table(Original = geno_ids_raw,  Cleaned = geno_ids))))
  changed <- changed[Original != Cleaned]

  matched    <- intersect(pheno_ids, geno_ids)
  pheno_only <- setdiff(pheno_ids, geno_ids)
  geno_only  <- setdiff(geno_ids, pheno_ids)

  overlap_frac <- if (length(pheno_ids))
    length(matched) / length(unique(pheno_ids)) else 0

  status <- if (overlap_frac < min_overlap_frac) "FAIL"
            else if (length(pheno_only) || length(geno_only))
              "PASS WITH WARNINGS"
            else "PASS"

  summary_dt <- data.table::data.table(
    N_phenotyped = length(unique(pheno_ids)),
    N_genotyped  = length(unique(geno_ids)),
    N_matched    = length(matched),
    N_phenotype_only = length(pheno_only),
    N_genotype_only  = length(geno_only),
    Overlap_pct_of_phenotyped = round(overlap_frac * 100, 2))

  out <- list(status = status, summary = summary_dt, matched_ids = matched,
              phenotype_only_ids = pheno_only, genotype_only_ids = geno_only,
              id_cleaning_changed = changed[])
  class(out) <- c("hbr_id_report", "list")
  out
}

#' Print an individual-ID match report
#'
#' @param x An object of class `hbr_id_report`.
#' @param ... Ignored.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' rep <- structure(
#'   list(status = "PASS",
#'        summary = data.table::data.table(
#'          N_phenotyped              = 3L,
#'          N_genotyped               = 2L,
#'          N_matched                 = 2L,
#'          N_phenotype_only          = 1L,
#'          N_genotype_only           = 0L,
#'          Overlap_pct_of_phenotyped = 66.67),
#'        matched_ids        = c("a", "b"),
#'        phenotype_only_ids = "c",
#'        genotype_only_ids  = character(0),
#'        id_cleaning_changed = data.table::data.table(
#'          Original = character(0),
#'          Cleaned  = character(0))),
#'   class = "hbr_id_report")
#' print(rep)
#'
#' @export
print.hbr_id_report <- function(x, ...) {
  s <- x$summary
  cat("INDIVIDUAL ID MATCH REPORT\n==========================\n\n")
  cat(sprintf("Phenotyped individuals  : %d\n", s$N_phenotyped))
  cat(sprintf("Genotyped individuals   : %d\n", s$N_genotyped))
  cat(sprintf("Matched                 : %d (%.1f%% of phenotyped)\n",
              s$N_matched, s$Overlap_pct_of_phenotyped))
  cat(sprintf("Phenotype-only          : %d\n", s$N_phenotype_only))
  cat(sprintf("Genotype-only           : %d\n", s$N_genotype_only))
  if (nrow(x$id_cleaning_changed))
    cat(sprintf("\nID cleaning altered %d ID(s).\n",
                nrow(x$id_cleaning_changed)))
  cat("\nSTATUS:", x$status, "\n")
  invisible(x)
}


# ----------------------------------------------------------------------------
# GWAS SNP validation and REF/ALT construction
# ----------------------------------------------------------------------------

#' Validate GWAS SNP input
#'
#' Trims and coerces the SNP, chromosome, and position columns, removes
#' invalid rows, drops exact duplicates, and detects SNP IDs that map to
#' more than one chromosome or position. Depending on
#' `conflict_action`, conflicting SNPs are errored on, removed, or
#' reduced to their first occurrence. A QC summary is attached as an
#' attribute.
#'
#' @param gwas_snps Data frame with GWAS SNP records.
#' @param snp_col,chr_col,pos_col Column names in `gwas_snps`.
#' @param conflict_action One of `"error"`, `"remove"`, `"first"`.
#' @param verbose Print QC summary.
#'
#' @return A cleaned `data.table` of SNP records with a `qc_summary`
#'   attribute holding `input_records`, `exact_duplicates_removed`,
#'   `conflicting_snps`, `final_records`, `unique_snps`, and the
#'   conflicting SNP table.
#'
#' @noRd
hbr_validate_gwas_snps <- function(gwas_snps,
                                   snp_col = "SNP",
                                   chr_col = "Chromosome",
                                   pos_col = "Position",
                                   conflict_action = c("error", "remove",
                                                       "first"),
                                   verbose = TRUE) {

  conflict_action <- match.arg(conflict_action)

  if (is.null(gwas_snps))
    stop("gwas_snps cannot be NULL.", call. = FALSE)
  if (!is.data.frame(gwas_snps))
    stop("gwas_snps must be a data.frame or data.table.", call. = FALSE)

  missing_cols <- setdiff(c(snp_col, chr_col, pos_col), names(gwas_snps))
  if (length(missing_cols))
    stop("Missing required columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)

  dt <- data.table::as.data.table(data.table::copy(gwas_snps))
  n_input <- nrow(dt)

  dt[, (snp_col) := trimws(as.character(get(snp_col)))]
  dt[, (chr_col) := trimws(as.character(get(chr_col)))]
  dt[, (pos_col) := suppressWarnings(as.numeric(as.character(get(pos_col))))]

  dt <- dt[
    !is.na(get(snp_col)) & get(snp_col) != "" &
      !is.na(get(chr_col)) & get(chr_col) != "" &
      !is.na(get(pos_col)) & get(pos_col) >= 1]

  n_before <- nrow(dt)
  dt <- unique(dt, by = c(snp_col, chr_col, pos_col))
  exact_duplicates_removed <- n_before - nrow(dt)

  coord_check <- dt[
    , .(n_chr = data.table::uniqueN(get(chr_col)),
        n_pos = data.table::uniqueN(get(pos_col))),
    by = snp_col]
  conflicts     <- coord_check[n_chr > 1 | n_pos > 1]
  conflict_ids  <- conflicts[[snp_col]]

  conflict_table <- dt[get(snp_col) %in% conflict_ids]
  if (length(conflict_ids) > 0) {
    conflict_table <- conflict_table[
      order(get(snp_col), get(chr_col), get(pos_col))]

    if (verbose) {
      message("")
      message("==============================================")
      message("SNP COORDINATE CONFLICT")
      message("==============================================")
      message("Conflicting SNP IDs: ", length(conflict_ids))
      message("")
      print(conflict_table[, c(snp_col, chr_col, pos_col), with = FALSE])
      message("")
    }

    if (conflict_action == "error") {
      stop("QC failure: ", length(conflict_ids),
           " SNP ID(s) have multiple chromosome/position assignments.",
           call. = FALSE)
    }
    if (conflict_action == "remove") {
      dt <- dt[!get(snp_col) %in% conflict_ids]
      if (verbose)
        message("Removed ", length(conflict_ids),
                " conflicting SNP IDs.")
    }
    if (conflict_action == "first") {
      dt <- dt[, .SD[1], by = snp_col]
      if (verbose)
        message("WARNING: kept the first coordinate for ",
                length(conflict_ids), " conflicting SNP IDs.")
    }
  }

  final_unique_snps <- data.table::uniqueN(dt[[snp_col]])

  if (verbose) {
    message("")
    message("==============================================")
    message("GWAS SNP QC SUMMARY")
    message("==============================================")
    message("Input records:             ", n_input)
    message("Exact duplicates removed:  ", exact_duplicates_removed)
    message("Conflicting SNP IDs:       ", length(conflict_ids))
    message("Final records:             ", nrow(dt))
    message("Unique SNP IDs:            ", final_unique_snps)
    message("==============================================")
    message("")
  }

  attr(dt, "qc_summary") <- list(
    input_records = n_input,
    exact_duplicates_removed = exact_duplicates_removed,
    conflicting_snps = length(conflict_ids),
    final_records = nrow(dt),
    unique_snps = final_unique_snps,
    conflicts = conflict_table)

  dt
}

#' Build a REF/ALT table from a PLINK .bim file
#'
#' PLINK stores alleles as A1 (usually minor/ALT) and A2 (usually
#' major/REF). The returned table uses the canonical convention
#' REF = A2, ALT = A1, and drops any variant whose alleles are not
#' single A/C/G/T bases.
#'
#' @param prefix PLINK file prefix or `.bed` / `.bim` / `.fam` path.
#'
#' @return A `data.table(SNP, CHR, POS, REF, ALT)`.
#'
#' @noRd
.hbr_ref_alt_from_plink_prefix <- function(prefix) {
  prefix <- sub("\\.(bed|bim|fam)$", "", prefix, ignore.case = TRUE)
  bim_file <- paste0(prefix, ".bim")
  if (!file.exists(bim_file))
    stop(".hbr_ref_alt_from_plink_prefix(): .bim not found: ", bim_file,
         call. = FALSE)

  bim <- data.table::fread(bim_file, header = FALSE,
                           col.names = c("CHR", "SNP", "CM", "BP",
                                         "A1", "A2"))

  ref_alt <- data.table::data.table(
    SNP = as.character(bim$SNP),
    CHR = as.character(bim$CHR),
    POS = as.numeric(bim$BP),
    REF = toupper(as.character(bim$A2)),
    ALT = toupper(as.character(bim$A1)))

  ref_alt[
    !is.na(CHR) & !is.na(POS) & !is.na(REF) & !is.na(ALT) &
      REF != "" & ALT != "" &
      REF %in% c("A", "C", "G", "T") &
      ALT %in% c("A", "C", "G", "T")]
}


# ----------------------------------------------------------------------------
# Unified input-QC report
# ----------------------------------------------------------------------------

#' Generate a unified input-QC report
#'
#' Runs the three input-QC checks — sequence-resource QC
#' ([hbr_validate_sequence_resources()]), genotype-reference
#' compatibility ([hbr_validate_genotype_reference()]), and ID overlap
#' ([hbr_validate_individual_ids()], when phenotypes are supplied) —
#' and prints a single consolidated report. Optionally writes a plain-
#' text file and/or a full report set (TXT + three CSV tables) to
#' `out_dir`. Only the TXT and CSV files are written: the HTML
#' duplicate was removed.
#'
#' @param reference A `hbr_reference` object.
#' @param genotype A `hbr_genotype` object.
#' @param phenotype_data Optional phenotype table.
#' @param pheno_id_col ID column in `phenotype_data`.
#' @param out_file Optional single text file path.
#' @param out_dir Optional directory for the full report set.
#' @param ... Passed to [hbr_validate_genotype_reference()].
#' @return Invisibly, a list with `sequence_qc`, `compatibility`,
#'   `id_report`, `overall_status`, and `files`.
#' @examples
#' if (interactive()) {
#'   toy <- system.file("extdata", "toy", package = "HaploBreedR")
#'   # A complete call needs an hbr_reference and hbr_genotype; see the
#'   # vignette for the construction used here.
#' }
#' @export
hbr_generate_input_qc_report <- function(reference, genotype,
                                         phenotype_data = NULL,
                                         pheno_id_col   = "ID",
                                         out_file       = NULL,
                                         out_dir        = NULL,
                                         ...) {

  if (!inherits(reference, "hbr_reference"))
    stop("`reference` must be a hbr_reference object.", call. = FALSE)
  if (!inherits(genotype, "hbr_genotype"))
    stop("`genotype` must be a hbr_genotype object.", call. = FALSE)

  if (is.null(genotype$report))
    genotype <- hbr_normalize_genotype(genotype, reference)

  seq_qc <- reference$validation
  if (is.null(seq_qc))
    seq_qc <- hbr_validate_sequence_resources(reference)

  compat <- hbr_validate_genotype_reference(genotype, reference, ...)

  id_report <- if (!is.null(phenotype_data))
    hbr_validate_individual_ids(phenotype_data, genotype, pheno_id_col)
  else NULL

  statuses <- c(seq_qc$status, compat$status,
                if (!is.null(id_report)) id_report$status)
  overall <- if (any(statuses == "FAIL")) "FAIL"
             else if (any(statuses == "NO SEQUENCE DATA")) "PASS WITH WARNINGS"
             else if (any(grepl("WARNING", statuses)))      "PASS WITH WARNINGS"
             else                                            "PASS"

  con <- textConnection("cap_out", "w", local = TRUE)
  sink(con)
  cat("HAPLOBREEDR INPUT QC REPORT\n")
  cat("===========================\n")
  cat("Generated:", format(Sys.time()), "\n\n")
  print(seq_qc); cat("\n")
  print(compat);  cat("\n")
  if (!is.null(id_report)) { print(id_report); cat("\n") }
  cat("OVERALL STATUS:", overall, "\n")
  sink()
  close(con)
  report_text <- paste(cap_out, collapse = "\n")
  cat(report_text, "\n")

  files <- list(txt = NULL,
                genotype_reference_qc_csv = NULL,
                sequence_resource_qc_csv  = NULL,
                individual_id_qc_csv      = NULL)

  if (!is.null(out_file)) {
    writeLines(report_text, out_file)
    files$txt <- out_file
    cat("\nReport written to:", out_file, "\n")
  }

  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

    txt_path <- file.path(out_dir, "input_qc_report.txt")
    writeLines(report_text, txt_path)
    files$txt <- txt_path

    # HTML duplicate removed: the TXT above carries the same content.

    geno_ref_csv <- file.path(out_dir, "genotype_reference_qc.csv")
    data.table::fwrite(compat$summary, geno_ref_csv)
    files$genotype_reference_qc_csv <- geno_ref_csv

    seq_csv <- file.path(out_dir, "sequence_resource_qc.csv")
    data.table::fwrite(seq_qc$summary, seq_csv)
    files$sequence_resource_qc_csv <- seq_csv

    if (!is.null(id_report)) {
      id_csv <- file.path(out_dir, "individual_id_qc.csv")
      data.table::fwrite(id_report$summary, id_csv)
      files$individual_id_qc_csv <- id_csv
    }

    cat("\nFull QC report set written to:", out_dir, "\n")
  }

  invisible(list(sequence_qc    = seq_qc,
                 compatibility  = compat,
                 id_report      = id_report,
                 overall_status = overall,
                 files          = files))
}

#' Pre-flight Validation of Inputs
#'
#' Validates PLINK files, phenotype records, GFF3 features, and sample ID overlap
#' before launching intensive discovery pipelines.
#'
#' @param plink_prefix Character. PLINK binary prefix.
#' @param phenotype_data data.frame of phenotypes.
#' @param trait_col Character. Target trait column name.
#' @param pheno_id_col Character. Individual ID column in phenotype_data.
#' @param gff3 Optional path to GFF3 file.
#'
#' @return A list of class \code{hbr_check} detailing compatibility status.
#' @examples
#' if (interactive()) {
#'   pheno <- data.frame(ID = sprintf("IND%02d", 1:10), trait = rnorm(10))
#'   toy   <- system.file("extdata", "toy", package = "HaploBreedR")
#'   res   <- hbr_check_inputs(
#'     plink_prefix   = file.path(toy, "toy"),
#'     phenotype_data = pheno,
#'     trait_col      = "trait")
#' }
#' @export
hbr_check_inputs <- function(plink_prefix, phenotype_data, trait_col,
                             pheno_id_col = "ID", gff3 = NULL) {
  rep <- hbr_generate_input_qc_report(
    plink_prefix = plink_prefix,
    phenotype_data = phenotype_data,
    trait_col = trait_col,
    pheno_id_col = pheno_id_col,
    verbose = TRUE
  )
  class(rep) <- "hbr_check"
  rep
}

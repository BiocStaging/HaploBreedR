# ============================================================================
# R/genotype.r
# Genotype class, PLINK reading (in-memory and BEDMatrix-streaming),
# normalization, reference validation, strand harmonization, per-gene QC,
# PLINK writing, and packing helpers.

# ----------------------------------------------------------------------------
# S3 print method
# ----------------------------------------------------------------------------

#' Print an hbr_genotype object
#'
#' @param x An object of class `hbr_genotype`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @examples
#' g <- structure(
#'   list(genotypes      = matrix(c(0, 1, 2), nrow = 1),
#'        variants       = data.frame(SNP = "s1", CHR = "1",
#'                                    POS = 100, REF = "A", ALT = "T"),
#'        individuals    = "ind1",
#'        source_format  = "PLINK",
#'        source_path    = NA_character_,
#'        declared_build = NA_character_,
#'        n_individuals  = 1L,
#'        n_variants     = 1L,
#'        streaming      = FALSE),
#'   class = "hbr_genotype")
#' print(g)
#'
#' @export
print.hbr_genotype <- function(x, ...) {
  mode <- if (isTRUE(x$streaming)) "BEDMatrix (streaming)" else "in-memory"
  cat("<hbr_genotype>", x$source_format, "|", mode, "\n")
  cat("  Individuals:", x$n_individuals,
      " Variants:", x$n_variants, "\n")
  if (!is.null(x$source_path))
    cat("  Source:", x$source_path, "\n")
  if (!is.null(x$row_keep) && length(x$row_keep) != x$n_individuals)
    cat("  Kept rows:", length(x$row_keep), "\n")
  if (length(x$swap_cols))
    cat("  Strand-swapped columns:", length(x$swap_cols), "\n")
  invisible(x)
}


# ----------------------------------------------------------------------------
# Internal low-level helpers
# ----------------------------------------------------------------------------

#' Number of rows in an hbr_genotype's genotype storage
#' @noRd
.hbr_geno_nrow <- function(g) {
  if (inherits(g, "hbr_genotype")) g$n_individuals else nrow(g)
}

#' Number of columns in an hbr_genotype's genotype storage
#' @noRd
.hbr_geno_ncol <- function(g) {
  if (inherits(g, "hbr_genotype")) g$n_variants else ncol(g)
}

#' Row names of an hbr_genotype
#' @noRd
.hbr_geno_rownames <- function(g) {
  if (inherits(g, "hbr_genotype")) g$row_names else rownames(g)
}

#' Column names of an hbr_genotype
#' @noRd
.hbr_geno_colnames <- function(g) {
  if (inherits(g, "hbr_genotype")) g$col_names else colnames(g)
}

#' Pull a subset of genotype columns (and rows) as a numeric matrix
#'
#' Reads only the requested columns from the underlying storage (a
#' BEDMatrix handle or a plain matrix), applies any strand-swap
#' correction recorded in `g$swap_cols`, and subsets rows with
#' `g$row_keep` (or the supplied `row_idx`).
#'
#' @param g An `hbr_genotype` or a plain numeric matrix.
#' @param col_idx Integer vector of column indices.
#' @param row_idx Optional integer vector of row indices.
#'
#' @return A numeric matrix with `NA` for missing calls.
#'
#' @noRd
.hbr_pull_genotypes <- function(g, col_idx, row_idx = NULL) {

  if (inherits(g, "hbr_genotype")) {
    if (is.null(row_idx)) row_idx <- g$row_keep
    swap_cols <- g$swap_cols
    handle    <- g$genotypes
  } else {
    swap_cols <- integer(0)
    handle    <- g
  }

  x <- handle[, col_idx, drop = FALSE]
  storage.mode(x) <- "double"

  # Apply 0<->2 flip on swapped columns
  if (length(swap_cols)) {
    flip <- which(col_idx %in% swap_cols)
    for (j in flip) {
      v <- x[, j]
      v[v == 0L] <- -1L
      v[v == 2L] <-  0L
      v[v == -1L] <- 2L
      x[, j] <- v
    }
  }

  if (!is.null(row_idx) && length(row_idx) != nrow(x))
    x <- x[row_idx, , drop = FALSE]
  
  if (inherits(g, "hbr_genotype")) {
    rn <- g$individuals
    if (length(rn) == nrow(x)) rownames(x) <- rn
  } else if (!is.null(row_idx) && length(row_idx) == nrow(x)) {
    rn <- rownames(g)
    if (!is.null(rn)) rownames(x) <- rn[row_idx]
  }
  
  x
}

#' Iterate over chunks of columns of an hbr_genotype
#'
#' Calls `FUN(chunk_matrix, start_col, end_col)` once per chunk of
#' columns, with the swap correction applied. Used by the PLINK writer
#' and the RDS packer.
#'
#' @noRd
.hbr_iterate_chunks <- function(g, chunk_cols = 20000L, FUN) {
  n_cols <- .hbr_geno_ncol(g)
  out <- vector("list", ceiling(n_cols / chunk_cols))
  k <- 0L
  for (s in seq.int(1L, n_cols, by = chunk_cols)) {
    e   <- min(n_cols, s + chunk_cols - 1L)
    idx <- s:e
    chunk <- .hbr_pull_genotypes(g, idx)
    k <- k + 1L
    out[[k]] <- FUN(chunk, s, e)
  }
  out[seq_len(k)]
}


# ----------------------------------------------------------------------------
# PLINK readers
# ----------------------------------------------------------------------------

#' Read PLINK files into an hbr_genotype (in-memory, legacy path)
#'
#' @param prefix PLINK file prefix.
#' @param plink_path Path to the PLINK executable.
#' @param verbose Print progress.
#'
#' @return A list with `genotypes`, `variants`, `individuals`.
#'
#' @noRd
.hbr_read_genotype_plink <- function(prefix, plink_path, verbose = TRUE) {

  vcat <- function(...) if (isTRUE(verbose)) cat(...)

  bim_file <- paste0(prefix, ".bim")
  fam_file <- paste0(prefix, ".fam")
  raw_file <- paste0(prefix, ".raw")

  if (!file.exists(bim_file)) stop("PLINK .bim not found: ", bim_file)
  if (!file.exists(fam_file)) stop("PLINK .fam not found: ", fam_file)

  if (!file.exists(raw_file)) {
    if (!file.exists(paste0(prefix, ".bed"))) {
      stop(".raw not found and no .bed present to generate it from: ",
           prefix)
    }
    vcat("  .raw not found; generating via PLINK --recode A...\n")
    .hbr_run_plink(plink_path,
                   c("--bfile", shQuote(prefix),
                     "--recode", "A", "--out", shQuote(prefix)))
    if (!file.exists(raw_file))
      stop("PLINK reported success but .raw was not created: ", raw_file)
  }

  bim <- data.table::fread(bim_file, header = FALSE)
  data.table::setnames(bim, c("CHR", "SNP", "CM", "BP", "A1", "A2"))

  fam <- data.table::fread(fam_file, header = FALSE)
  data.table::setnames(fam, c("FID", "IID", "FATHER", "MOTHER",
                              "SEX", "PHENOTYPE"))

  vcat("  Reading .raw file: ", basename(raw_file), "...\n")

  raw <- tryCatch(
    data.table::fread(
      raw_file, sep = " ", header = TRUE, fill = TRUE,
      strip.white = TRUE, blank.lines.skip = TRUE, check.names = FALSE),
    error = function(e) NULL)

  if (is.null(raw) || nrow(raw) == 0L) {
    vcat("  fread failed or returned 0 rows; falling back to read.table...\n")
    raw <- data.table::as.data.table(
      read.table(raw_file, header = TRUE, sep = "", comment.char = "",
                 quote = "", stringsAsFactors = FALSE, check.names = FALSE))
  }

  if (nrow(raw) == 0L)
    stop("Failed to read ", raw_file,
         " with both fread() and read.table().", call. = FALSE)

  ind_ids   <- as.character(raw$IID)
  meta_cols <- c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE")
  geno_cols <- setdiff(names(raw), meta_cols)
  n_ind     <- nrow(raw)
  n_snp     <- length(geno_cols)

  vcat(sprintf("  Converting to matrix: %d individuals x %d SNPs (%.1f MB)\n",
               n_ind, n_snp, n_ind * n_snp * 8 / 1e6))

  chunk_size <- 100000L
  n_chunks   <- ceiling(n_snp / chunk_size)

  geno_mat <- matrix(NA_real_, nrow = n_ind, ncol = n_snp)

  for (chunk in seq_len(n_chunks)) {
    start_col  <- (chunk - 1L) * chunk_size + 1L
    end_col    <- min(chunk * chunk_size, n_snp)
    cols_chunk <- geno_cols[seq.int(start_col, end_col)]

    sub_dt <- raw[, ..cols_chunk]
    geno_mat[, seq.int(start_col, end_col)] <- as.matrix(sub_dt)
    rm(sub_dt); invisible(gc(verbose = FALSE))
  }

  rownames(geno_mat) <- ind_ids
  colnames(geno_mat) <- sub("_[^_]+$", "", geno_cols)

  rm(raw); invisible(gc(verbose = FALSE))

  variants <- bim[match(colnames(geno_mat), SNP),
                  .(SNP = as.character(SNP), CHR = as.character(CHR),
                    POS = as.numeric(BP),
                    REF = as.character(A2), ALT = as.character(A1))]

  list(genotypes = geno_mat, variants = variants, individuals = ind_ids)
}

#' Open PLINK files as a streaming (BEDMatrix) hbr_genotype
#'
#' Returns an `hbr_genotype` object whose `$genotypes` slot is a
#' `BEDMatrix::BEDMatrix` handle when the package is installed, or a
#' fully-materialised numeric matrix otherwise. The object also carries
#' the fields the pipeline needs downstream: `variants`, `individuals`,
#' `swap_cols` (strand-swap column indices), `row_keep` (subset rows),
#' `canonical_map` (canonical SNP name -> column index), `col_names`,
#' `row_names`, and `streaming` flag.
#'
#' @param prefix PLINK file prefix.
#' @param plink_path Path to the PLINK executable. Only used when
#'   BEDMatrix is unavailable and the `.raw` file must be generated.
#' @param verbose Print a one-line mode report.
#'
#' @return An `hbr_genotype` object.
#'
#' @noRd
.hbr_open_genotype_plink <- function(prefix, plink_path = "plink",
                                      verbose = TRUE) {

  vcat <- function(...) if (isTRUE(verbose)) cat(...)

  bim_file <- paste0(prefix, ".bim")
  fam_file <- paste0(prefix, ".fam")
  bed_file <- paste0(prefix, ".bed")

  if (!file.exists(bim_file)) stop("PLINK .bim not found: ", bim_file)
  if (!file.exists(fam_file)) stop("PLINK .fam not found: ", fam_file)
  if (!file.exists(bed_file)) stop("PLINK .bed not found: ", bed_file)

  bim <- data.table::fread(bim_file, header = FALSE)
  data.table::setnames(bim, c("CHR", "SNP", "CM", "BP", "A1", "A2"))

  fam <- data.table::fread(fam_file, header = FALSE)
  data.table::setnames(fam, c("FID", "IID", "FATHER", "MOTHER",
                              "SEX", "PHENOTYPE"))

  inds <- as.character(fam$IID)
  variants <- data.table::data.table(
    SNP = as.character(bim$SNP),
    CHR = as.character(bim$CHR),
    POS = as.numeric(bim$BP),
    REF = as.character(bim$A2),
    ALT = as.character(bim$A1))

  n_ind <- length(inds)
  n_snp <- nrow(variants)

  if (requireNamespace("BEDMatrix", quietly = TRUE)) {
    vcat(sprintf("  Streaming: %d individuals x %d SNPs (BEDMatrix)\n",
                 n_ind, n_snp))
    handle <- suppressMessages(
      BEDMatrix::BEDMatrix(prefix, simple_names = TRUE))

    return(structure(list(
      genotypes      = handle,
      individuals    = inds,
      variants       = variants,
      n_individuals  = n_ind,
      n_variants     = n_snp,
      swap_cols      = integer(0),
      row_keep       = seq_len(n_ind),
      canonical_map  = NULL,
      col_names      = as.character(bim$SNP),
      row_names      = inds,
      streaming      = TRUE,
      source_format  = "PLINK",
      source_path    = prefix,
      declared_build = NA_character_,
      report         = NULL),
      class = c("hbr_genotype", "list")))
  }

  vcat(sprintf("  In-memory: %d individuals x %d SNPs (BEDMatrix not installed)\n",
               n_ind, n_snp))
  gm <- .hbr_read_genotype_plink(prefix, plink_path, verbose = FALSE)

  structure(list(
    genotypes      = gm$genotypes,
    individuals    = gm$individuals,
    variants       = variants,
    n_individuals  = nrow(gm$genotypes),
    n_variants     = ncol(gm$genotypes),
    swap_cols      = integer(0),
    row_keep       = seq_len(nrow(gm$genotypes)),
    canonical_map  = NULL,
    col_names      = colnames(gm$genotypes),
    row_names      = rownames(gm$genotypes),
    streaming      = FALSE,
    source_format  = "PLINK",
    source_path    = prefix,
    declared_build = NA_character_,
    report         = NULL),
    class = c("hbr_genotype", "list"))
}


# ----------------------------------------------------------------------------
# Normalization
# ----------------------------------------------------------------------------

#' Normalize a genotype object's chromosome names, alleles, and variant keys
#'
#' @param genotype A `hbr_genotype` object.
#' @param reference Optional `hbr_reference` object.
#' @param verbose Print the normalization summary.
#'
#' @return The input `hbr_genotype` with normalized `$variants` and an
#'   added `$report` element.
#'
#' @noRd
hbr_normalize_genotype <- function(genotype, reference = NULL,
                                    verbose = TRUE) {

  if (!inherits(genotype, "hbr_genotype"))
    stop("hbr_normalize_genotype(): `genotype` must be a hbr_genotype.",
         call. = FALSE)
  if (!is.null(reference) && !inherits(reference, "hbr_reference")) {
    stop("hbr_normalize_genotype(): `reference` must be a hbr_reference ",
         "object or NULL.", call. = FALSE)
  }

  v <- data.table::copy(genotype$variants)
  v[, CHR_original := CHR]
  v[, REF_original := REF]
  v[, ALT_original := ALT]

  norm_chr <- function(x) {
    x <- trimws(as.character(x))
    x <- sub("^(?i)chr(omosome)?[_ ]?", "", x, perl = TRUE)
    x <- sub("^0+(?=[0-9])", "", x, perl = TRUE)
    x
  }
  v[, CHR := norm_chr(CHR_original)]

  chr_map <- unique(v[, .(CHR_original, CHR_normalized = CHR)])
  n_chr_normalized <- sum(chr_map$CHR_original != chr_map$CHR_normalized)

  unrecognized_chromosomes <- character(0)
  if (!is.null(reference)) {
    ref_chrs <- as.character(reference$chromosomes$Chr)
    unrecognized_chromosomes <- sort(unique(v$CHR[!v$CHR %in% ref_chrs]))
  }

  v[, REF := toupper(trimws(REF_original))]
  v[, ALT := toupper(trimws(ALT_original))]
  is_clean_base <- function(x) !is.na(x) & grepl("^[ACGT]$", x)
  v[, is_biallelic_snp := is_clean_base(REF) & is_clean_base(ALT)]
  n_non_biallelic <- sum(!v$is_biallelic_snp)

  v[, POS := suppressWarnings(as.numeric(POS))]
  v[, variant_key := paste(CHR, POS, REF, ALT, sep = ":")]
  dup_keys <- v$variant_key[duplicated(v$variant_key) |
                              duplicated(v$variant_key, fromLast = TRUE)]
  n_duplicate <- length(unique(dup_keys))

  genotype$variants <- v[]
  genotype$report <- list(
    n_chr_normalized         = n_chr_normalized,
    chr_map                  = chr_map[],
    unrecognized_chromosomes = unrecognized_chromosomes,
    n_non_biallelic_snp      = n_non_biallelic,
    n_duplicate_variant_key  = n_duplicate)

  if (isTRUE(verbose))
    cat("Normalized genotype variants: ",
        n_chr_normalized, " chromosome name(s) changed, ",
        length(unrecognized_chromosomes), " unrecognized, ",
        n_non_biallelic, " non-biallelic-SNP allele(s) flagged, ",
        n_duplicate, " duplicate coordinate(s)\n", sep = "")

  genotype
}


# ----------------------------------------------------------------------------
# Genotype-reference validation
# ----------------------------------------------------------------------------

#' Validate a genotype dataset's compatibility with a reference genome
#'
#' @param genotype A `hbr_genotype` object.
#' @param reference A `hbr_reference` object.
#' @param min_chr_match_frac Minimum fraction of variants whose
#'   chromosome must be present in the reference.
#' @param min_ref_match_frac Minimum REF-allele match fraction.
#' @param max_unknown_chr_frac Maximum tolerated fraction of unrecognized
#'   chromosomes.
#' @param verbose Print step-by-step progress.
#'
#' @return A list of class `hbr_compat_report`.
#'
#' @noRd
hbr_validate_genotype_reference <- function(genotype, reference,
                                            min_chr_match_frac   = 0.95,
                                            min_ref_match_frac   = 0.90,
                                            max_unknown_chr_frac = 0.05,
                                            verbose = TRUE) {

  if (!inherits(genotype, "hbr_genotype"))
    stop("hbr_validate_genotype_reference(): `genotype` must be a ",
         "hbr_genotype object.")
  if (!inherits(reference, "hbr_reference"))
    stop("hbr_validate_genotype_reference(): `reference` must be a ",
         "hbr_reference object.")

  t_step <- function(t0, label)
    if (verbose)
      cat(sprintf("    [validate] %-28s %.2fs\n", label,
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))

  t0 <- Sys.time()
  if (is.null(genotype$report)) {
    if (verbose)
      cat("    Genotype not pre-normalized; running ",
          "hbr_normalize_genotype()...\n", sep = "")
    genotype <- hbr_normalize_genotype(genotype, reference,
                                        verbose = verbose)
  }
  t_step(t0, "normalize")

  v <- data.table::as.data.table(genotype$variants)
  v <- v[, .(SNP = as.character(SNP),
             CHR = as.character(CHR),
             POS = as.numeric(POS),
             REF = toupper(as.character(REF)),
             is_biallelic_snp)]
  n_total  <- nrow(v)
  ref_chrs <- as.character(reference$chromosomes$Chr)

  t0 <- Sys.time()
  chrs_seen    <- unique(v$CHR)
  chrs_matched <- intersect(chrs_seen, ref_chrs)
  in_ref       <- v$CHR %in% ref_chrs
  n_chr_matched_variants <- sum(in_ref)
  chr_match_frac <- if (n_total) n_chr_matched_variants / n_total else 0
  t_step(t0, "chr match")

  t0 <- Sys.time()
  chr_len <- stats::setNames(as.numeric(reference$chromosomes$Length_Approx),
                             as.character(reference$chromosomes$Chr))
  chr_length_v   <- chr_len[v$CHR]
  position_valid <- !is.na(v$POS) & v$POS > 0 &
    (is.na(chr_length_v) | v$POS <= chr_length_v)
  bad_pos        <- !is.na(chr_length_v) & !position_valid
  position_out_of_range <- v[bad_pos,
                             .(SNP, CHR, POS,
                               chr_length = chr_length_v[bad_pos])]
  valid_pos_frac <- if (n_total) mean(position_valid, na.rm = TRUE) else 0
  t_step(t0, "position range")

  n_multiallelic <- sum(!v$is_biallelic_snp, na.rm = TRUE)
  n_duplicate    <- genotype$report$n_duplicate_variant_key

  ref_mismatches <- data.table::data.table(
    SNP = character(), CHR = character(), POS = numeric(),
    REF_stated = character(), REF_reference = character())
  n_ref_checked <- 0L
  n_ref_match   <- 0L
  coverage_mode <- "unavailable"

  checkable <- v[!is.na(CHR) & !is.na(POS) &
                   CHR %in% names(reference$genome_seqs)]

  if (isTRUE(reference$has_genome_fasta) && nrow(checkable)) {
    coverage_mode <- "genome_fasta"
    t0 <- Sys.time()

    ref_base <- rep(NA_character_, nrow(checkable))
    chr_names_ok <- intersect(unique(checkable$CHR),
                              names(reference$genome_seqs))

    for (chr in chr_names_ok) {
      chr_seq <- reference$genome_seqs[[chr]]
      if (is.null(chr_seq)) next
      idx <- which(checkable$CHR == chr)
      pos       <- as.integer(round(checkable$POS[idx]))
      chr_chars <- as.character(chr_seq)
      chr_len_s <- nchar(chr_chars)
      ok <- !is.na(pos) & pos >= 1L & pos <= chr_len_s
      if (!any(ok)) { rm(chr_chars); next }
      ref_base[idx[ok]] <- substring(chr_chars, pos[ok], pos[ok])
      rm(chr_chars)
    }

    checkable[, REF_reference := ref_base]
    n_ref_checked <- sum(!is.na(checkable$REF_reference))
    match_vec     <- !is.na(checkable$REF_reference) &
      checkable$REF == checkable$REF_reference
    n_ref_match   <- sum(match_vec, na.rm = TRUE)
    ref_mismatches <- checkable[!is.na(REF_reference) &
                                  REF != REF_reference,
                                .(SNP, CHR, POS,
                                  REF_stated = REF, REF_reference)]

    if (n_ref_checked == 0L)
      warning("REF check ran but matched 0 variants. ",
              "Chromosome names in reference$genome_seqs do not match ",
              "genotype$variants$CHR.", call. = FALSE)
    t_step(t0, sprintf("REF lookup (%d var)", nrow(checkable)))

  } else if (!is.null(reference$cds) && isTRUE(reference$has_cds)) {
    coverage_mode <- "cds_overlap_only"
    cds_plus <- data.table::as.data.table(reference$cds)[strand == "+"]

    if (nrow(cds_plus)) {
      t0 <- Sys.time()
      cds_dt <- cds_plus[, .(Chr = as.character(Chr),
                             cds_start = as.numeric(start),
                             cds_end = as.numeric(end),
                             transcript_id = as.character(transcript_id))]
      vv <- v[, .(SNP, CHR, POS, REF)]

      hits <- cds_dt[vv,
                     on = .(Chr = CHR, cds_start <= POS, cds_end >= POS),
                     nomatch = NULL, allow.cartesian = TRUE,
                     .(SNP = i.SNP, CHR = i.CHR, POS = i.POS, REF = i.REF,
                       transcript_id = x.transcript_id,
                       exon_offset = i.POS - x.cds_start)]
      t_step(t0, sprintf("CDS overlap join (%d hits)", nrow(hits)))

      if (nrow(hits) && !is.null(reference$cds_seqs)) {
        hits <- unique(hits, by = "SNP")
        hits <- hits[transcript_id %in% names(reference$cds_seqs)]

        if (nrow(hits)) {
          t0 <- Sys.time()
          ref_base <- rep(NA_character_, nrow(hits))

          for (tx in unique(hits$transcript_id)) {
            seq <- reference$cds_seqs[[tx]]
            if (is.null(seq)) next
            idx       <- which(hits$transcript_id == tx)
            pos       <- as.integer(round(hits$exon_offset[idx])) + 1L
            seq_chars <- as.character(seq)
            seq_len_s <- nchar(seq_chars)
            ok <- !is.na(pos) & pos >= 1L & pos <= seq_len_s
            if (!any(ok)) { rm(seq_chars); next }
            ref_base[idx[ok]] <- substring(seq_chars, pos[ok], pos[ok])
            rm(seq_chars)
          }

          hits[, REF_reference := ref_base]
          n_ref_checked <- sum(!is.na(hits$REF_reference))
          n_ref_match   <- sum(!is.na(hits$REF_reference) &
                                 hits$REF == hits$REF_reference,
                               na.rm = TRUE)
          ref_mismatches <- hits[!is.na(REF_reference) &
                                   REF != REF_reference,
                                 .(SNP, CHR, POS,
                                   REF_stated = REF, REF_reference)]
          t_step(t0, sprintf("CDS REF lookup (%d hits)", nrow(hits)))
        }
      }
    }
  }

  ref_match_frac <- if (n_ref_checked) n_ref_match / n_ref_checked else NA_real_

  unknown_chr_frac <- 1 - chr_match_frac
  fail <- (chr_match_frac < min_chr_match_frac) ||
    (!is.na(ref_match_frac) && ref_match_frac < min_ref_match_frac) ||
    (unknown_chr_frac > max_unknown_chr_frac)

  any_warning <- n_multiallelic > 0 || n_duplicate > 0 ||
    nrow(position_out_of_range) > 0 ||
    length(genotype$report$unrecognized_chromosomes) > 0 ||
    nrow(ref_mismatches) > 0

  status <- if (fail) "FAIL" else
    if (any_warning) "PASS WITH WARNINGS" else "PASS"

  summary_dt <- data.table::data.table(
    Reference_genome            = reference$genome,
    Reference_annotation        = reference$annotation,
    Declared_genotype_build     = genotype$declared_build,
    Individuals                 = genotype$n_individuals,
    Total_variants              = n_total,
    Chromosomes_matched         = length(chrs_matched),
    Chromosomes_total_reference = length(ref_chrs),
    Chromosome_match_pct        = round(chr_match_frac * 100, 2),
    Valid_position_pct          = round(valid_pos_frac * 100, 2),
    REF_check_coverage          = coverage_mode,
    REF_checked_n               = n_ref_checked,
    REF_match_pct               = if (is.na(ref_match_frac)) NA_real_ else
      round(ref_match_frac * 100, 2),
    Duplicate_positions         = n_duplicate,
    Multiallelic_variants       = n_multiallelic,
    Reference_mismatches        = nrow(ref_mismatches),
    Unknown_chromosomes         = length(genotype$report$unrecognized_chromosomes))

  out <- list(status                   = status,
              summary                  = summary_dt,
              ref_check_coverage       = coverage_mode,
              position_out_of_range    = position_out_of_range,
              ref_mismatches           = ref_mismatches,
              unrecognized_chromosomes =
                genotype$report$unrecognized_chromosomes)
  class(out) <- c("hbr_compat_report", "list")
  out
}

#' Print a genotype-reference compatibility report
#'
#' @param x An object of class `hbr_compat_report`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @examples
#' x <- structure(
#'   list(status = "PASS",
#'        summary = data.table::data.table(
#'          Reference_genome = "TOY", Reference_annotation = "CUSTOM",
#'          Declared_genotype_build = NA_character_,
#'          Individuals = 1L, Total_variants = 1L,
#'          Chromosomes_matched = 1L, Chromosomes_total_reference = 1L,
#'          Chromosome_match_pct = 100, Valid_position_pct = 100,
#'          REF_check_coverage = "genome_fasta",
#'          REF_checked_n = 1L, REF_match_pct = 100,
#'          Duplicate_positions = 0L, Multiallelic_variants = 0L,
#'          Reference_mismatches = 0L, Unknown_chromosomes = 0L)),
#'   class = "hbr_compat_report")
#' print(x)
#'
#' @export
print.hbr_compat_report <- function(x, ...) {
  cat("GENOTYPE-REFERENCE COMPATIBILITY\n")
  cat("================================\n\n")
  s <- x$summary
  cat(sprintf("Reference genome       : %s (%s)\n",
              s$Reference_genome, s$Reference_annotation))
  cat(sprintf("Declared genotype build: %s\n",
              s$Declared_genotype_build))
  cat(sprintf("Individuals            : %d\n", s$Individuals))
  cat(sprintf("Total variants         : %d\n", s$Total_variants))
  cat(sprintf("Chromosome match       : %d / %d (%.2f%%)\n",
              s$Chromosomes_matched, s$Chromosomes_total_reference,
              s$Chromosome_match_pct))
  cat(sprintf("Valid positions        : %.2f%%\n", s$Valid_position_pct))
  cat(sprintf("REF check mode         : %s\n", s$REF_check_coverage))
  cat(sprintf("REF match              : %d checked, %.2f%%\n",
              s$REF_checked_n,
              if (is.na(s$REF_match_pct)) NA else s$REF_match_pct))
  cat(sprintf("Duplicate positions    : %d\n", s$Duplicate_positions))
  cat(sprintf("Multiallelic variants  : %d\n", s$Multiallelic_variants))
  cat(sprintf("Reference mismatches   : %d\n", s$Reference_mismatches))
  cat(sprintf("Unknown chromosomes    : %d\n\n", s$Unknown_chromosomes))
  cat("STATUS:", x$status, "\n")
  invisible(x)
}


# ----------------------------------------------------------------------------
# Strand harmonization
# ----------------------------------------------------------------------------

#' Classify REF/ALT mismatches against the reference genome
#'
#' @param ref_alt A data.frame with `SNP`, `CHR`, `POS`, `REF`, `ALT`.
#' @param reference A `hbr_reference` object with `genome_seqs`.
#' @param verbose Print the per-class tally.
#'
#' @return A `data.table` copy of `ref_alt` with `REF_ref` and `class`.
#'
#' @noRd
.hbr_classify_ref_mismatches <- function(ref_alt, reference,
                                          verbose = TRUE) {

  if (is.null(reference$genome_seqs))
    stop(".hbr_classify_ref_mismatches(): reference has no genome_seqs.")

  dt <- data.table::as.data.table(data.table::copy(ref_alt))
  dt[, SNP := as.character(SNP)]
  dt[, CHR := as.character(CHR)]
  dt[, POS := suppressWarnings(as.numeric(POS))]
  dt[, REF := toupper(as.character(REF))]
  dt[, ALT := toupper(sub(",.*$", "", as.character(ALT)))]

  dt[, REF_ref := NA_character_]
  chr_names <- names(reference$genome_seqs)

  for (chr in intersect(unique(dt$CHR), chr_names)) {
    idx <- which(dt$CHR == chr)
    if (!length(idx)) next
    seq_chars <- as.character(reference$genome_seqs[[chr]])
    seq_len   <- nchar(seq_chars)
    pos       <- as.integer(round(dt$POS[idx]))
    ok <- !is.na(pos) & pos >= 1L & pos <= seq_len
    if (!any(ok)) next
    dt$REF_ref[idx[ok]] <- toupper(substring(seq_chars, pos[ok], pos[ok]))
    rm(seq_chars)
  }

  comp <- c(A = "T", T = "A", C = "G", G = "C")

  dt[, class := data.table::fcase(
    is.na(REF_ref),                                             "no_ref",
    REF == REF_ref,                                             "match",
    (REF %in% c("A","T") & ALT %in% c("A","T")) |
      (REF %in% c("C","G") & ALT %in% c("G","C")) |
      (REF %in% c("A","T") & ALT %in% c("T","A")) |
      (REF %in% c("C","G") & ALT %in% c("C","G")),               "ambiguous",
    !is.na(comp[REF]) & comp[REF] == REF_ref & REF != ALT,       "strand_flip",
    ALT == REF_ref & REF != ALT,                                 "allele_swap",
    !is.na(comp[ALT]) & comp[ALT] == REF_ref & REF != ALT,       "flip_and_swap",
    default = "unresolvable")]

  if (isTRUE(verbose)) {
    tab <- dt[, .N, by = class][order(-N)]
    tab[, pct := round(100 * N / sum(N), 2)]
    cat("    Strand classification:\n")
    for (i in seq_len(nrow(tab)))
      cat(sprintf("      %-16s %10s (%.2f%%)\n",
                  tab$class[i],
                  format(tab$N[i], big.mark = ","),
                  tab$pct[i]))
  }
  dt[]
}

#' Harmonize REF/ALT in a hbr_genotype against the reference genome
#'
#' Updates REF/ALT labels and records the indices of columns whose
#' dosage should be flipped (0<->2) at read time in `genotype$swap_cols`.
#' The underlying genotype matrix (or BEDMatrix handle) is never
#' mutated, so peak memory is unchanged. When `drop_ambiguous = TRUE`
#' on an in-memory genotype, the ambiguous columns are physically
#' dropped; on a streaming genotype the call falls back to marking the
#' columns as dropped via `genotype$drop_cols`.
#'
#' @param genotype A `hbr_genotype` object.
#' @param reference A `hbr_reference` object with `genome_seqs`.
#' @param drop_ambiguous Drop A/T and C/G variants after harmonization.
#' @param verbose Print progress.
#'
#' @return The genotype object with harmonized variants and updated
#'   `swap_cols` / `drop_cols`.
#'
#' @noRd
hbr_harmonize_strand <- function(genotype, reference,
                                  drop_ambiguous = FALSE,
                                  verbose = TRUE) {

  if (!inherits(genotype, "hbr_genotype"))
    stop("`genotype` must be a hbr_genotype object.")
  if (is.null(reference$genome_seqs))
    stop("`reference` has no genome_seqs; cannot harmonize.")

  cls <- .hbr_classify_ref_mismatches(genotype$variants, reference,
                                       verbose = verbose)

  comp <- c(A = "T", T = "A", C = "G", G = "C")

  idx_flip <- which(cls$class %in% c("strand_flip", "flip_and_swap"))
  idx_swap <- which(cls$class %in% c("allele_swap", "flip_and_swap"))

  if (length(idx_flip)) {
    genotype$variants$REF[idx_flip] <- comp[genotype$variants$REF[idx_flip]]
    genotype$variants$ALT[idx_flip] <- comp[genotype$variants$ALT[idx_flip]]
  }
  if (length(idx_swap)) {
    tmp_ref <- genotype$variants$REF[idx_swap]
    genotype$variants$REF[idx_swap] <- genotype$variants$ALT[idx_swap]
    genotype$variants$ALT[idx_swap] <- tmp_ref
    # Record for read-time flip; do NOT touch the matrix
    genotype$swap_cols <- sort(unique(c(genotype$swap_cols, idx_swap)))
  }

  genotype$strand_report <- cls

  if (isTRUE(verbose)) {
    cat(sprintf("    Corrections applied: %d flip, %d swap, %d flip+swap\n",
                sum(cls$class == "strand_flip"),
                sum(cls$class == "allele_swap"),
                sum(cls$class == "flip_and_swap")))
  }

  if (isTRUE(drop_ambiguous)) {
    idx_amb <- which(cls$class == "ambiguous")
    if (length(idx_amb)) {
      if (isTRUE(genotype$streaming)) {
        warning("drop_ambiguous = TRUE with streaming genotype: ambiguous ",
                "columns will be dropped logically (no matrix copy).",
                call. = FALSE)
        genotype$drop_cols <- sort(unique(c(genotype$drop_cols %||%
                                              integer(0), idx_amb)))
      } else {
        genotype$genotypes  <- genotype$genotypes[, -idx_amb, drop = FALSE]
        genotype$row_names  <- rownames(genotype$genotypes)
        genotype$col_names  <- colnames(genotype$genotypes)
        genotype$n_variants <- length(genotype$col_names)
        # Reindex swap_cols
        if (length(genotype$swap_cols)) {
          keep_mask <- rep(TRUE, nrow(genotype$variants))
          keep_mask[idx_amb] <- FALSE
          old_to_new <- cumsum(keep_mask)
          genotype$swap_cols <- old_to_new[genotype$swap_cols]
          genotype$swap_cols <- genotype$swap_cols[
            genotype$swap_cols > 0]
        }
      }
      genotype$variants <- genotype$variants[-idx_amb]
      if (!genotype$streaming)
        genotype$n_variants <- nrow(genotype$variants)
    }
  }

  genotype
}


# ----------------------------------------------------------------------------
# Per-gene genotype QC
# ----------------------------------------------------------------------------

#' Extract and filter a genotype sub-matrix for a gene's functional SNPs
#'
#' Given a `functional_snps` data.table and either a plain numeric
#' matrix or an `hbr_genotype` streaming source, extracts the relevant
#' columns and applies SNP-level and individual-level QC. In streaming
#' mode only the gene's columns are ever read from disk.
#'
#' @param functional_snps A `data.table` with `SNP` and `Genomic_Order`.
#' @param geno_mat Either a numeric genotype matrix with canonical
#'   column names, or an `hbr_genotype` source (streaming or in-memory).
#' @param max_miss_snp Maximum per-SNP missingness.
#' @param maf_min Minimum minor-allele frequency.
#' @param max_het_ind Maximum per-individual heterozygote fraction.
#' @param max_miss_ind Maximum per-individual missingness.
#' @param min_individuals Minimum individuals after QC.
#' @param gene_label Optional gene label for QC-diagnostic prints.
#' @param verbose Print the `[QC-diag]` filter diagnostics.
#'
#' @return A list with `geno`, `functional_snps`, `n_snps`,
#'   `n_individuals`, and per-filter counts.
#'
#' @noRd
hbr_extract_genotypes <- function(functional_snps, geno_mat,
                                  max_miss_snp = 0.20, maf_min = 0.05,
                                  max_het_ind = 0.10, max_miss_ind = 0.10,
                                  min_individuals = 10,
                                  gene_label = NULL,
                                  verbose = TRUE) {

  msg <- function(...)
    if (isTRUE(verbose) && !is.null(gene_label))
      cat("    [QC-diag]", gene_label, ..., "\n")

  snps <- functional_snps[order(Genomic_Order)]

  if (inherits(geno_mat, "hbr_genotype")) {
    gs <- geno_mat
    if (is.null(gs$canonical_map))
      stop("hbr_extract_genotypes(): hbr_genotype has no canonical_map. ",
           "Did the pipeline canonicalize the genotype?",
           call. = FALSE)
    bed_cols <- gs$canonical_map[as.character(snps$SNP)]
    ok <- !is.na(bed_cols)
    if (!any(ok)) {
      msg("no functional SNPs present in the canonical map"); return(NULL)
    }
    G <- .hbr_pull_genotypes(gs, unname(bed_cols[ok]))
    colnames(G) <- as.character(snps$SNP[ok])
    snps <- snps[ok]
  } else {
    col_idx <- which(colnames(geno_mat) %in% snps$SNP)
    if (!length(col_idx)) {
      msg("no SNP IDs found in genotype matrix"); return(NULL)
    }
    G <- geno_mat[, col_idx, drop = FALSE]
    colnames(G) <- colnames(geno_mat)[col_idx]
  }

  miss_snp  <- apply(G, 2, function(x) sum(is.na(x)) / length(x))
  n_rm_miss <- sum(miss_snp > max_miss_snp)
  G <- G[, miss_snp <= max_miss_snp, drop = FALSE]
  if (!ncol(G)) { msg("all SNPs removed by missingness"); return(NULL) }

  maf_v    <- apply(G, 2, hbr_calc_maf)
  n_rm_maf <- sum(maf_v < maf_min)
  G <- G[, maf_v >= maf_min, drop = FALSE]
  if (!ncol(G)) { msg("all SNPs removed by MAF <", maf_min); return(NULL) }

  mono1     <- apply(G, 2, function(x)
    length(unique(stats::na.omit(x))) == 1)
  n_rm_mono <- sum(mono1)
  G <- G[, !mono1, drop = FALSE]
  if (!ncol(G)) { msg("all SNPs monomorphic"); return(NULL) }

  ind_het  <- apply(G, 1, hbr_calc_het)
  keep_het <- ind_het <= max_het_ind
  G[G == 1] <- NA
  mono2     <- apply(G, 2, function(x)
    length(unique(stats::na.omit(x))) == 1)
  n_rm_mono <- n_rm_mono + sum(mono2)
  G <- G[, !mono2, drop = FALSE]
  if (!ncol(G)) {
    msg("all SNPs monomorphic after heterozygote removal"); return(NULL)
  }

  miss_ind  <- apply(G, 1, function(x) sum(is.na(x)) / length(x))
  keep_miss <- miss_ind <= max_miss_ind
  keep_ind  <- keep_het & keep_miss
  G <- G[keep_ind, , drop = FALSE]
  if (nrow(G) < min_individuals) {
    msg("only", nrow(G), "individuals pass QC"); return(NULL)
  }

  snps_kept <- snps[SNP %in% colnames(G)]
  snps_kept <- snps_kept[order(Genomic_Order)]
  final_idx <- match(snps_kept$SNP, colnames(G))
  if (anyNA(final_idx) || !length(final_idx)) {
    msg("internal SNP/column mismatch after filtering")
    return(NULL)
  }
  G <- G[, final_idx, drop = FALSE]
  colnames(G) <- snps_kept$SNP

  list(
    geno            = G,
    functional_snps = snps_kept,
    n_snps          = ncol(G),
    n_individuals   = nrow(G),
    n_rm_miss_snp   = n_rm_miss,
    n_rm_maf        = n_rm_maf,
    n_rm_mono       = n_rm_mono,
    n_rm_het_ind    = sum(!keep_het),
    n_rm_miss_ind   = sum(!keep_miss))
}


# ----------------------------------------------------------------------------
# PLINK writer
# ----------------------------------------------------------------------------

#' Write a hbr_genotype object out as PLINK files
#'
#' Writes PED and MAP files from the in-memory or streaming genotype,
#' then converts to PLINK binary format (`--make-bed`) and to a `.raw`
#' dosage file (`--recode A`). In streaming mode only one chunk of
#' columns is materialised at a time. The intermediate PED and MAP files
#' are removed after conversion.
#'
#' @param genotype A `hbr_genotype` object.
#' @param out_prefix Output PLINK prefix.
#' @param plink_path Path to the PLINK executable.
#' @param verbose Print progress.
#'
#' @return `out_prefix`, invisibly.
#'
#' @noRd
hbr_write_plink <- function(genotype, out_prefix, plink_path = "plink",
                            verbose = TRUE) {

  vcat <- function(...) if (isTRUE(verbose)) cat(...)

  if (!inherits(genotype, "hbr_genotype"))
    stop("hbr_write_plink(): `genotype` must be a hbr_genotype object.",
         call. = FALSE)

  ids <- genotype$row_names[genotype$row_keep]
  v   <- genotype$variants

  # If drop_ambiguous was used in streaming mode, exclude dropped columns
  keep_cols <- setdiff(seq_len(nrow(v)),
                        genotype$drop_cols %||% integer(0))
  v_kept <- v[keep_cols]

  ped_file <- paste0(out_prefix, ".ped")
  map_file <- paste0(out_prefix, ".map")

  vcat("Writing PLINK PED/MAP for", length(ids), "individuals x",
       nrow(v_kept), "variants...\n")

  map_dt <- data.table::data.table(
    CHR = v_kept$CHR, SNP = v_kept$SNP, CM = 0,
    POS = as.integer(round(v_kept$POS)))
  data.table::fwrite(map_dt, map_file, sep = "\t", col.names = FALSE)

  # Open output file once
  con <- file(ped_file, "w")
  on.exit(try(close(con), silent = TRUE), add = TRUE)

  chunk_cols <- 20000L
  n_chunks   <- ceiling(length(keep_cols) / chunk_cols)
  k <- 0L
  for (s in seq.int(1L, length(keep_cols), by = chunk_cols)) {
    e <- min(length(keep_cols), s + chunk_cols - 1L)
    idx <- keep_cols[s:e]
    k <- k + 1L
    if (verbose && n_chunks > 1L)
      cat(sprintf("  chunk %d / %d\n", k, n_chunks))

    g  <- .hbr_pull_genotypes(genotype, idx)
    vc <- v[idx]

    ref <- vc$REF; alt <- vc$ALT; nr <- nrow(g); nc <- ncol(g)

    ref_mat <- matrix(rep(ref, each = nr), nrow = nr, ncol = nc)
    alt_mat <- matrix(rep(alt, each = nr), nrow = nr, ncol = nc)

    a1 <- ref_mat; a1[g == 2]   <- alt_mat[g == 2]; a1[is.na(g)] <- "0"
    a2 <- ref_mat; a2[g == 1]   <- alt_mat[g == 1]
    a2[g == 2]   <- alt_mat[g == 2]; a2[is.na(g)] <- "0"

    ap <- matrix(NA_character_, nrow = nr, ncol = 2L * nc)
    ap[, seq(1L, 2L * nc, by = 2L)] <- a1
    ap[, seq(2L, 2L * nc, by = 2L)] <- a2

    lines <- paste(ids, ids, 0, 0, 0, -9,
                   apply(ap, 1, paste, collapse = " "))
    writeLines(lines, con)
  }
  close(con)

  vcat("Converting to PLINK binary format...\n")
  .hbr_run_plink(plink_path,
                 c("--file", shQuote(out_prefix),
                   "--make-bed", "--out", shQuote(out_prefix)))
  if (!file.exists(paste0(out_prefix, ".bed")))
    stop("PLINK reported success but .bed was not created: ",
         out_prefix, call. = FALSE)

  .hbr_run_plink(plink_path,
                 c("--bfile", shQuote(out_prefix),
                   "--recode", "A", "--out", shQuote(out_prefix)))
  if (!file.exists(paste0(out_prefix, ".raw")))
    stop("PLINK reported success but .raw was not created: ",
         out_prefix, call. = FALSE)

  unlink(c(ped_file, map_file))

  invisible(out_prefix)
}


# ----------------------------------------------------------------------------
# Subset in memory
# ----------------------------------------------------------------------------

#' Subset a hbr_genotype to phenotyped individuals
#'
#' Computes the row index vector of phenotyped individuals and stores it
#' in `genotype$row_keep`. In streaming mode no matrix is touched; in
#' in-memory mode the matrix is physically subset (fast, one-off).
#'
#' @param genotype A `hbr_genotype` object.
#' @param phenotype_data Data frame with an ID column.
#' @param pheno_id_col ID column in `phenotype_data`.
#' @param min_overlap_frac Minimum overlap; stops below.
#' @param verbose Print progress.
#'
#' @return A list with `genotype` and `id_report`.
#'
#' @noRd
.hbr_subset_genotype <- function(genotype, phenotype_data,
                                 pheno_id_col = "ID",
                                 min_overlap_frac = 0.5,
                                 verbose = TRUE) {

  id_report <- hbr_validate_individual_ids(phenotype_data, genotype,
                                           pheno_id_col, min_overlap_frac)
  if (verbose) print(id_report)
  if (id_report$status == "FAIL")
    stop(".hbr_subset_genotype(): individual-ID overlap below threshold (",
         id_report$summary$Overlap_pct_of_phenotyped, "% matched).",
         call. = FALSE)

  genotype$row_names   <- hbr_clean_id(genotype$row_names)
  geno_ids_clean       <- genotype$row_names
  
  keep <- which(geno_ids_clean %in% id_report$matched_ids)
  if (!length(keep))
    stop(".hbr_subset_genotype(): no individuals remain after matching.")
  
  genotype$row_keep      <- keep
  genotype$individuals   <- geno_ids_clean[keep]
  genotype$n_individuals <- length(keep)
  
  if (!isTRUE(genotype$streaming)) {
    genotype$genotypes <- genotype$genotypes[keep, , drop = FALSE]
    
    rownames(genotype$genotypes) <- genotype$row_names
  }

  if (verbose)
    cat("Subset to", genotype$n_individuals,
        "phenotyped individuals.\n")

  list(genotype = genotype, id_report = id_report)
}
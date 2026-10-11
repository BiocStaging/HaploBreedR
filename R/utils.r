# ============================================================================
# R/utils.r
# Cross-cutting helpers: ID cleaning, safe coercion, PLINK invocation,
# logging, genetics utilities, and HGVS formatting.
#
# ----------------------------------------------------------------------------
# Core helpers
# ----------------------------------------------------------------------------

#' Null-coalescing operator
#'
#' Returns `b` when `a` is `NULL` or a scalar `NA`, otherwise `a`.
#'
#' @param a,b Values.
#'
#' @return `a` or `b`.
#'
#' @noRd
`%||%` <- function(a, b) {
  if (is.null(a) || (length(a) == 1 && is.na(a))) b else a
}

#' Log elapsed time
#'
#' Prints elapsed time since `start_time` in seconds, minutes, or hours.
#'
#' @param start_time Start time.
#' @param step Label to print.
#' @param verbose Print.
#'
#' @return Elapsed time in seconds, invisibly.
#'
#' @noRd
.hbr_log_time <- function(start_time, step, verbose = TRUE) {
  if (!isTRUE(verbose)) return(invisible(NULL))
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  if (elapsed < 60) {
    cat(sprintf("    [time] %-42s %6.2f s\n", step, elapsed))
  } else if (elapsed < 3600) {
    cat(sprintf("    [time] %-42s %6.2f s  (%.1f min)\n",
                step, elapsed, elapsed / 60))
  } else {
    cat(sprintf("    [time] %-42s %6.2f s  (%.2f h)\n",
                step, elapsed, elapsed / 3600))
  }
  invisible(elapsed)
}

#' Run PLINK
#'
#' Wraps [system2()] and errors on non-zero exit status.
#'
#' @param plink_path Path to the PLINK executable.
#' @param args Argument vector.
#'
#' @return `TRUE`, invisibly.
#'
#' @noRd
.hbr_run_plink <- function(plink_path, args) {
  status <- tryCatch(
    system2(plink_path, args = args),
    error = function(e) -1L)
  if (!identical(status, 0L) && !identical(status, 0)) {
    stop("PLINK call failed (exit status ", status,
         "). Args: ", paste(args, collapse = " "),
         ". Confirm `plink_path` is correct and PLINK is installed.",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Clean a GFF3 gene ID
#'
#' Strips `gene:`, `ID=`, and trailing `;...` from attribute strings.
#'
#' @param x Character vector.
#'
#' @return Character vector.
#'
#' @noRd
.hbr_clean_gene_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^gene:", "", x, ignore.case = TRUE)

  has_gene_id <- !is.na(x) & grepl("(^|;)gene_id=", x, ignore.case = TRUE)
  if (any(has_gene_id)) {
    x[has_gene_id] <- sub(".*(?:^|;)gene_id=([^;]+).*", "\\1",
                          x[has_gene_id], perl = TRUE)
  }

  has_id <- !is.na(x) & !has_gene_id &
    grepl("(^|;)ID=", x, ignore.case = TRUE)
  if (any(has_id)) {
    x[has_id] <- sub(".*(?:^|;)ID=([^;]+).*", "\\1",
                     x[has_id], perl = TRUE)
  }

  x <- sub(";.*$", "", x)
  x <- sub("^gene:", "", x, ignore.case = TRUE)
  x <- gsub('^"|"$', "", x)
  x <- trimws(x)
  x[is.na(x) | x == "" | x %in% c("NA", "NaN", "NULL", ".")] <- NA_character_
  x
}

#' Drop all-blank columns
#'
#' Removes columns that are entirely `NA`, empty, or `"."`.
#'
#' @param dt A data.table.
#' @param keep_always Columns to keep regardless.
#'
#' @return A data.table.
#'
#' @noRd
.hbr_drop_blank_cols <- function(dt, keep_always = character(0)) {
  dt <- data.table::as.data.table(dt)
  if (!nrow(dt) || !ncol(dt)) return(dt)

  is_blank <- vapply(dt, function(x) {
    if (is.list(x)) {
      return(all(vapply(x, function(e) {
        is.null(e) ||
          (length(e) == 1L && (is.na(e) || !nzchar(as.character(e))))
      }, logical(1))))
    }
    if (is.character(x)) {
      return(all(is.na(x) | !nzchar(x) |
                   toupper(trimws(x)) %in% c("NA", "NAN", ".", "NULL")))
    }
    if (is.numeric(x) || is.logical(x)) {
      return(all(is.na(x) | !is.finite(x)))
    }
    FALSE
  }, logical(1))

  drop <- is_blank & !names(dt) %in% keep_always
  if (any(drop)) dt <- dt[, !drop, with = FALSE]
  dt
}

#' Clean an individual ID to a canonical form
#'
#' Trims whitespace, replaces hyphens with underscores, and collapses
#' runs of whitespace to single underscores.
#'
#' @param x Character vector of IDs.
#'
#' @return Character vector, same length as `x`.
#'
#' @noRd
hbr_clean_id <- function(x) {
  x <- trimws(as.character(x))
  gsub("\\s+", "_", gsub("-", "_", x))
}

#' Safe scalar extraction
#'
#' Returns `NA_real_` if the input is not a finite scalar.
#'
#' @param x Any object.
#'
#' @return Numeric scalar.
#'
#' @noRd
hbr_safe_num <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  if (is.list(x)) {
    for (v in x) if (is.numeric(v) && length(v) == 1) return(v)
    return(NA_real_)
  }
  if (is.numeric(x) && length(x) == 1) return(x)
  NA_real_
}

#' First element of each list entry
#'
#' Vectorised `[[1]]` for list columns.
#'
#' @param x List or vector.
#'
#' @return Character vector.
#'
#' @noRd
hbr_first_chr <- function(x) {
  xl <- as.list(x)
  vapply(xl, function(v) if (length(v)) as.character(v[1]) else NA_character_,
         character(1))
}

#' Strip an ID prefix
#'
#' Removes `gene:` / `transcript:` / `mRNA:` / `CDS:`.
#'
#' @param x Character vector.
#'
#' @return Character vector.
#'
#' @noRd
hbr_strip_prefix <- function(x) gsub("^(gene|transcript|mRNA|CDS):", "", x)


# ----------------------------------------------------------------------------
# Genetics utilities
# ----------------------------------------------------------------------------

#' Minor-allele frequency from a 0/1/2 genotype vector
#'
#' @param x Numeric vector.
#'
#' @return Numeric scalar in [0, 0.5].
#'
#' @noRd
hbr_calc_maf <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(0)
  p <- mean(x) / 2
  min(p, 1 - p)
}

#' Heterozygosity rate from a 0/1/2 genotype vector
#'
#' @param x Numeric vector.
#'
#' @return Numeric scalar in [0, 1].
#'
#' @noRd
hbr_calc_het <- function(x) {
  nm <- !is.na(x)
  if (!sum(nm)) return(1)
  sum(x[nm] == 1) / sum(nm)
}

#' Benjamini-Hochberg FDR correction with NA-safe handling
#'
#' @param p_values Numeric vector.
#' @param method Adjustment method.
#'
#' @return Adjusted p-values.
#'
#' @noRd
hbr_apply_fdr_correction <- function(p_values, method = "BH") {
  if (is.null(p_values) || length(p_values) == 0) return(numeric(0))
  na_pos <- is.na(p_values)
  p_clean <- p_values[!na_pos]
  if (length(p_clean) == 0) return(p_values)
  adj_clean <- stats::p.adjust(p_clean, method = method)
  out <- p_values
  out[!na_pos] <- adj_clean
  out
}

#' Standardize a chromosome name
#'
#' Strips `chr`/`chromosome` prefixes, collapses whitespace, normalizes
#' leading zeros.
#'
#' @param x Character vector.
#'
#' @return Character vector.
#'
#' @noRd
hbr_standardize_chr <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^chromosome[_\\.:-]*", "", x, ignore.case = TRUE)
  x <- sub("^chr[_\\.:-]*",        "", x, ignore.case = TRUE)
  x <- gsub("\\s+", "", x)
  is_numeric <- grepl("^[0-9]+$", x)
  x[is_numeric] <- as.character(as.integer(x[is_numeric]))
  x
}

#' Normalize an OSA1-style rice gene ID
#'
#' @param x Character vector.
#'
#' @examples
#' hbr_normalize_osa1_gene_id("LOC_Os01g0100100")
#' hbr_normalize_osa1_gene_id("Os01g0100100.1")
#' @return Character vector.
#'
#' @noRd
hbr_normalize_osa1_gene_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^LOC_", "", x)
  x <- sub("\\.[0-9]+$", "", x)
  x <- sub("-[0-9]+$",   "", x)
  x
}

#' Map one-letter amino-acid codes to three-letter names
#'
#' @param x Character vector of one-letter codes.
#'
#' @return Character vector of three-letter names.
#'
#' @noRd
hbr_aa_full_name <- function(x) {
  map <- c(
    A = "Ala", R = "Arg", N = "Asn", D = "Asp", C = "Cys",
    Q = "Gln", E = "Glu", G = "Gly", H = "His", I = "Ile",
    L = "Leu", K = "Lys", M = "Met", F = "Phe", P = "Pro",
    S = "Ser", T = "Thr", W = "Trp", Y = "Tyr", V = "Val",
    `*` = "Stop")
  out <- map[x]
  out[is.na(out)] <- x[is.na(out)]
  unname(out)
}

#' Build an amino-acid change label
#'
#' @param aa_ref Reference amino acid.
#' @param aa_alt Alternative amino acid.
#' @param codon_num Codon number.
#' @param full_name Use three-letter names.
#'
#' @return Character scalar.
#'
#' @noRd
hbr_aa_change_label <- function(aa_ref, aa_alt, codon_num, full_name = FALSE) {
  if (isTRUE(full_name)) {
    paste0(hbr_aa_full_name(aa_ref), codon_num, hbr_aa_full_name(aa_alt))
  } else {
    paste0(aa_ref, codon_num, aa_alt)
  }
}

#' Display label for a variant
#'
#' @param snps_dt A data.table with `Region`, `REF`, `ALT`, and possibly
#'   `AA_REF`, `AA_ALT`, `Codon_Num`.
#' @param full_name Use three-letter amino-acid codes.
#'
#' @return Character vector.
#'
#' @noRd
hbr_variant_display_label <- function(snps_dt, full_name = TRUE) {
  vapply(seq_len(nrow(snps_dt)), function(j) {
    region <- snps_dt$Region[j]
    if (!is.na(region) && region == "CDS_nonsynonymous") {
      hbr_aa_change_label(snps_dt$AA_REF[j], snps_dt$AA_ALT[j],
                          snps_dt$Codon_Num[j], full_name = full_name)
    } else if (!is.na(region)) {
      paste0(snps_dt$Region[j], ":", snps_dt$REF[j], ">", snps_dt$ALT[j])
    } else {
      paste0(snps_dt$REF[j], ">", snps_dt$ALT[j])
    }
  }, character(1))
}

#' Pack a genotype matrix to one byte per call
#'
#' Converts a numeric or integer dosage matrix (values 0, 1, 2; `NA`
#' allowed) into a raw matrix with 1 byte per cell (`255` encodes `NA`).
#' Returns the input unchanged if the matrix contains non-integer values
#' outside [0, 254]. Packing is done in column blocks to bound peak
#' memory.
#'
#' @param m Numeric or integer matrix of genotype dosages.
#' @param chunk_cols Number of columns per block. `NULL` (default) picks
#'   a block size targeting ~5 million cells per block.
#'
#' @return A raw matrix of class `hbr_packed_geno`, or `m` unchanged.
#'
#' @noRd
hbr_pack_genotypes <- function(m, chunk_cols = NULL) {
  if (is.null(m) || !is.matrix(m) || is.raw(m) || !is.numeric(m)) return(m)
  nr <- nrow(m); nc <- ncol(m)
  if (nr == 0L || nc == 0L) return(m)
  if (is.null(chunk_cols)) chunk_cols <- max(1L, as.integer(5e6 %/% nr))
  chunk_cols <- max(1L, as.integer(chunk_cols))

  r <- raw(as.double(nr) * as.double(nc))
  dim(r) <- c(nr, nc)
  for (st in seq.int(1L, nc, by = chunk_cols)) {
    en  <- min(nc, st + chunk_cols - 1L)
    blk <- m[, seq.int(st, en), drop = FALSE]
    if (!is.integer(blk) && any(blk != round(blk), na.rm = TRUE)) return(m)
    if (any(blk < 0 | blk > 254, na.rm = TRUE)) return(m)
    x <- as.integer(blk)
    x[is.na(x)] <- 255L
    r[, seq.int(st, en)] <- as.raw(x)
  }
  dimnames(r) <- dimnames(m)
  structure(r, class = "hbr_packed_geno")
}

#' Restore a genotype matrix packed by [hbr_pack_genotypes()]
#'
#' Unpacks a `hbr_packed_geno` matrix back to integer or numeric dosage.
#' Unpacking is done in column blocks to bound peak memory. Matrices
#' that are not packed are returned unchanged, so the function is safe
#' to call on any genotype matrix.
#'
#' @param x A packed matrix (or an ordinary matrix, returned unchanged).
#' @param as_integer Return an integer matrix (default) or numeric.
#' @param chunk_cols Columns per block. `NULL` (default) targets
#'   ~5 million cells per block.
#'
#' @return A matrix with `NA` for missing calls.
#'
#' @examples
#' # A plain matrix is returned unchanged
#' m <- matrix(c(0, 1, 2, NA), nrow = 2)
#' hbr_unpack_genotypes(m)
#'
#' @noRd
hbr_unpack_genotypes <- function(x, as_integer = TRUE, chunk_cols = NULL) {
  if (!inherits(x, "hbr_packed_geno")) return(x)
  d  <- dim(x); dn <- dimnames(x)
  xr <- unclass(x)
  if (is.null(chunk_cols))
    chunk_cols <- max(1L, as.integer(5e6 %/% max(1L, d[1])))
  out <- if (isTRUE(as_integer)) matrix(NA_integer_, d[1], d[2])
         else matrix(NA_real_, d[1], d[2])
  for (st in seq.int(1L, d[2], by = chunk_cols)) {
    en <- min(d[2], st + chunk_cols - 1L)
    v  <- as.integer(xr[, seq.int(st, en), drop = FALSE])
    v[v == 255L] <- NA_integer_
    out[, seq.int(st, en)] <- if (isTRUE(as_integer)) v else as.numeric(v)
  }
  dimnames(out) <- dn
  out
}

#' Write `pipeline_result_for_breeding.rds`
#'
#' Accepts either an in-memory numeric matrix or an `hbr_genotype`
#' streaming source in the `geno_mat` argument. When the argument is an
#' `hbr_genotype` the packing is streamed from the underlying BEDMatrix
#' so peak memory stays at one chunk. The function never loses the run:
#' if packing fails, the RDS is written without the genotype matrix and
#' a warning says so; the matrix can always be re-read from
#' `plink_prefix`.
#'
#' @param geno_mat Numeric matrix OR `hbr_genotype`.
#' @inheritParams .hbr_write_breeding_rds (see code)
#'
#' @return A list with `rds`, `reference_rds`, `genotype_saved`.
#'
#' @noRd
.hbr_write_breeding_rds <- function(superior_results, geno_mat, reference,
                                    ref_paths, plink_prefix, pipeline_summary,
                                    output_dir,
                                    rds_compress  = c("fast", "none", "xz"),
                                    rds_reference = c("embedded", "separate", "none"),
                                    verbose = TRUE) {
  
  rds_compress  <- match.arg(rds_compress)
  rds_reference <- match.arg(rds_reference)
  t0       <- Sys.time()
  rds_path <- file.path(output_dir, "pipeline_result_for_breeding.rds")
  ref_path <- file.path(output_dir, "pipeline_reference.rds")
  
  save_rds_fast <- function(obj, path) {
    if (rds_compress == "none") {
      saveRDS(obj, path, compress = FALSE)
    } else if (rds_compress == "xz") {
      saveRDS(obj, path, compress = "xz")
    } else {
      con <- gzfile(path, "wb", compression = 1L)
      on.exit(try(close(con), silent = TRUE), add = TRUE)
      serialize(obj, con)
      close(con)
    }
    invisible(path)
  }
  
  is_seq <- if (is.list(reference))
    vapply(reference, function(x) inherits(x, "XStringSet"), logical(1))
  else logical(0)
  ref_main <- reference; dropped_ref <- character(0)
  ref_file_out <- NA_character_
  if (rds_reference != "embedded" && is.list(reference)) {
    ref_main    <- reference[!is_seq]
    dropped_ref <- names(reference)[is_seq]
    if (rds_reference == "separate") {
      save_rds_fast(reference, ref_path)
      ref_file_out <- ref_path
    }
  }
  
  invisible(gc(verbose = FALSE))
  
  geno_ids <- if (inherits(geno_mat, "hbr_genotype"))
    geno_mat$individuals
  else if (!is.null(geno_mat)) rownames(geno_mat)
  else NULL
  
  geno_dim <- if (inherits(geno_mat, "hbr_genotype"))
    c(geno_mat$n_individuals, length(geno_mat$canonical_snp_names))
  else if (!is.null(geno_mat)) dim(geno_mat)
  else NULL
  
  geno_saved <- !is.null(geno_mat)
  packed <- tryCatch(
    hbr_pack_genotypes(geno_mat),
    error = function(e) {
      warning("Genotype packing failed (", conditionMessage(e),
              "): the RDS will not contain the genotype matrix.",
              call. = FALSE)
      geno_saved <<- FALSE
      NULL
    })
  
  build <- function(g) list(
    superior_results  = superior_results,
    geno_mat          = g,
    geno_packed       = inherits(g, "hbr_packed_geno"),
    geno_saved        = geno_saved && !is.null(g),
    geno_ids          = geno_ids,
    geno_dim          = geno_dim,
    plink_prefix      = plink_prefix,
    reference         = ref_main,
    reference_mode    = rds_reference,
    reference_dropped = dropped_ref,
    reference_file    = ref_file_out,
    reference_paths   = ref_paths,
    pipeline_summary  = pipeline_summary)
  
  ok <- tryCatch({ save_rds_fast(build(packed), rds_path); TRUE },
                 error = function(e) {
                   warning("Writing the RDS with the genotype matrix failed (",
                           conditionMessage(e), "); retrying without it.",
                           call. = FALSE)
                   FALSE })
  if (!ok) {
    geno_saved <- FALSE
    save_rds_fast(build(NULL), rds_path)
    warning("RDS saved WITHOUT the genotype matrix; re-read it from ",
            "`plink_prefix` ('", plink_prefix, "').", call. = FALSE)
  }
  
  if (verbose) {
    cat(sprintf("\nBreeding RDS saved: %s (%.1f MB)%s\n", rds_path,
                file.size(rds_path) / 1e6,
                if (geno_saved) "" else "  [no genotype matrix]"))
    if (!is.na(ref_file_out))
      cat(sprintf("Reference RDS saved: %s (%.1f MB)\n", ref_file_out,
                  file.size(ref_file_out) / 1e6))
    cat(sprintf("RDS time: %.1fs | compression: %s | reference: %s\n",
                as.numeric(difftime(Sys.time(), t0, units = "secs")),
                rds_compress, rds_reference))
  }
  invisible(list(rds = rds_path, reference_rds = ref_file_out,
                 genotype_saved = geno_saved))
}



#' Load the LD-region to gene mapping from a pipeline output directory
#'
#' @param output_dir Trait directory.
#'
#' @return A list with `gene_region`, `region_meta`, `genes_per_region`,
#'   or `NULL`.
#'
#' @noRd
.hbr_load_ld_regions <- function(output_dir) {

  sm_path <- file.path(output_dir, "publication_table_functional.csv")
  if (!file.exists(sm_path))
    sm_path <- file.path(output_dir, "superior_genes",
                         "publication_table_functional.csv")
  if (!file.exists(sm_path)) return(NULL)

  tbl <- data.table::fread(sm_path, na.strings = c("", "NA"))
  if (!"LD_Region_ID" %in% names(tbl)) return(NULL)

  gc <- intersect(c("RAPDB_Gene_ID", "MSU7_Gene_ID", "MSU_Gene_ID",
                    "Gene_ID", "Gene"), names(tbl))[1]
  if (is.na(gc)) return(NULL)

  cols <- unique(c(gc, "LD_Region_ID", "LD_Region_Start", "LD_Region_End",
                   "Chr", "Is_Superior", "Superior_Haplotype",
                   "P_value", "Cohens_d"))
  cols <- intersect(cols, names(tbl))

  out <- tbl[, ..cols]
  data.table::setnames(out, gc, "Gene", skip_absent = TRUE)
  for (cn in c("LD_Region_Start", "LD_Region_End", "Chr", "Superior_Haplotype",
               "P_value", "Cohens_d"))
    if (!cn %in% names(out)) out[, (cn) := NA]
  if (!"Is_Superior" %in% names(out)) out[, Is_Superior := NA]
  out <- out[!is.na(LD_Region_ID) & nzchar(as.character(LD_Region_ID))]
  if (!nrow(out)) return(NULL)

  exploded <- out[, {
    ids <- trimws(strsplit(as.character(LD_Region_ID), ";")[[1]])
    ids <- ids[nzchar(ids)]
    if (!length(ids)) NULL
    else data.table::data.table(
      Gene = Gene, LD_Region_ID = ids,
      LD_Region_Start = LD_Region_Start,
      LD_Region_End   = LD_Region_End,
      Chr             = if ("Chr" %in% names(.SD)) Chr else NA_character_)
  }, by = seq_len(nrow(out))][, seq_len := NULL]

  exploded <- unique(exploded, by = c("Gene", "LD_Region_ID"))

  meta <- out[, .(Gene, Is_Superior, Superior_Haplotype, P_value, Cohens_d)]
  meta <- unique(meta, by = "Gene")
  exploded <- merge(exploded, meta, by = "Gene", all.x = TRUE)

  region_meta <- unique(exploded[!is.na(LD_Region_ID) & nzchar(LD_Region_ID),
                                 .(LD_Region_ID, Chr,
                                   LD_Region_Start, LD_Region_End)])
  region_meta[, LD_Region_Length := LD_Region_End - LD_Region_Start + 1]

  genes_per_region <- exploded[
    Is_Superior == TRUE,
    .(N_Superior_Genes_In_Region = data.table::uniqueN(Gene),
      Superior_Genes = paste(sort(unique(Gene)), collapse = ";")),
    by = LD_Region_ID]

  all_genes_per_region <- exploded[
    , .(N_Total_Genes_In_Region = data.table::uniqueN(Gene)),
    by = LD_Region_ID]

  genes_per_region <- merge(genes_per_region, all_genes_per_region,
                            by = "LD_Region_ID", all.x = TRUE)

  list(gene_region      = exploded,
       region_meta      = region_meta,
       genes_per_region = genes_per_region)
}

#' Per-donor LD-region summary
#'
#' @param donor_details Per-donor gene table.
#' @param ld_regions Output of [`.hbr_load_ld_regions()`].
#'
#' @return A list with `per_region` and `per_donor`.
#'
#' @noRd
.hbr_donor_region_summary <- function(donor_details, ld_regions) {

  if (is.null(donor_details) || !nrow(donor_details) ||
      is.null(ld_regions)) {
    return(list(per_region = data.table::data.table(),
                per_donor  = data.table::data.table()))
  }

  merged <- merge(donor_details,
                  ld_regions$gene_region[, .(Gene, LD_Region_ID)],
                  by = "Gene", all.x = TRUE)
  merged <- merged[!is.na(LD_Region_ID) & nzchar(LD_Region_ID)]
  if (!nrow(merged)) {
    return(list(per_region = data.table::data.table(),
                per_donor  = data.table::data.table()))
  }

  merged <- merge(merged, ld_regions$genes_per_region,
                  by = "LD_Region_ID", all.x = TRUE)

  per_region <- merged[, .(
    N_Superior_In_Region       = data.table::uniqueN(Gene),
    N_Superior_Genes_In_Region = N_Superior_Genes_In_Region[1],
    N_Total_Genes_In_Region    = N_Total_Genes_In_Region[1],
    Complete_Region = data.table::uniqueN(Gene) ==
      N_Superior_Genes_In_Region[1],
    Genes_Carried = paste(sort(unique(Gene)), collapse = "; ")
  ), by = .(Donor, LD_Region_ID)]

  per_donor <- per_region[, .(
    N_LD_Regions                = .N,
    N_Complete_LD_Regions       = sum(Complete_Region, na.rm = TRUE),
    N_Total_Superior_In_Regions = sum(N_Superior_In_Region, na.rm = TRUE),
    N_Total_Region_Genes        = sum(N_Superior_Genes_In_Region,
                                       na.rm = TRUE),
    Complete_Regions            = if (any(Complete_Region))
      paste(LD_Region_ID[Complete_Region], collapse = "; ")
      else NA_character_,
    Regions_Carried             = paste(LD_Region_ID, collapse = "; ")
  ), by = Donor]

  list(per_region = per_region, per_donor = per_donor)
}


# ----------------------------------------------------------------------------
# HGVS formatting
# ----------------------------------------------------------------------------

#' One-to-three-letter amino-acid code map
#'
#' @noRd
.hbr_aa_1to3 <- c(
  A = "Ala", R = "Arg", N = "Asn", D = "Asp", C = "Cys", Q = "Gln",
  E = "Glu", G = "Gly", H = "His", I = "Ile", L = "Leu", K = "Lys",
  M = "Met", F = "Phe", P = "Pro", S = "Ser", T = "Thr", W = "Trp",
  Y = "Tyr", V = "Val", `*` = "Ter")

#' Format a single amino-acid change token as HGVS
#'
#' @param token Character scalar.
#'
#' @return Character scalar or `NA_character_`.
#'
#' @noRd
.hbr_format_hgvs_token <- function(token) {
  token <- trimws(token)
  if (!nzchar(token) || grepl(":", token)) return(NA_character_)

  m3 <- regmatches(token,
                    regexec("^([A-Za-z]{3})([0-9]+)([A-Za-z]{3})$",
                            token))[[1]]
  if (length(m3) == 4) return(paste0("p.", m3[2], m3[3], m3[4]))

  m1 <- regmatches(token,
                    regexec("^([A-Za-z])([0-9]+)([A-Za-z])$",
                            token))[[1]]
  if (length(m1) == 4) {
    ref3 <- .hbr_aa_1to3[toupper(m1[2])]
    alt3 <- .hbr_aa_1to3[toupper(m1[4])]
    if (!is.na(ref3) && !is.na(alt3))
      return(paste0("p.", ref3, m1[3], alt3))
  }
  NA_character_
}

#' Format a semicolon-joined Functional_Changes string as HGVS
#'
#' @param s Character scalar.
#'
#' @return A single string, or `NA_character_`.
#'
#' @noRd
.hbr_format_hgvs_string <- function(s) {
  if (is.na(s) || !nzchar(s)) return(NA_character_)
  toks <- strsplit(s, ";")[[1]]
  hgvs <- vapply(toks, .hbr_format_hgvs_token, character(1))
  hgvs <- hgvs[!is.na(hgvs)]
  if (!length(hgvs)) return(NA_character_)
  paste(hgvs, collapse = "; ")
}

#' Run an expression with a fixed RNG seed, restoring state afterwards
#'
#' BiocCheck discourages direct \code{set.seed()} calls because they
#' overwrite the user's RNG state. This helper snapshots
#' \code{.Random.seed}, sets the requested seed for the duration of
#' \code{expr}, and restores the original state on exit.
#'
#' @param seed Integer seed.
#' @param expr Expression to evaluate.
#'
#' @return The value of \code{expr}.
#'
#' @noRd
.with_seed <- function(seed, expr) {
  withr::with_seed(seed, force(expr))
}

#' Summary of the bundled rice genetic map
#'
#' @return A one-row data.frame with the number of anchors, chromosomes,
#'   and total genetic length.
#' @examples
#' hbr_rice_genetic_map_summary()
#' @noRd
hbr_rice_genetic_map_summary <- function() {
  m <- hbr_rice_genetic_map()
  per <- stats::aggregate(cm ~ chr, data = m,
                          FUN = function(v) diff(range(v)))
  names(per)[2] <- "length_cM"
  data.frame(
    n_markers  = nrow(m),
    n_chrom    = length(unique(m$chr)),
    total_cM   = round(sum(per$length_cM), 1),
    source     = attr(m, "source"),
    stringsAsFactors = FALSE)
}

#' Format a scalar for printing
#'
#' Returns the first element as a character string, or `"NA"` when the
#' value is `NULL`, zero-length, or `NA`. Internal helper used by the
#' print methods of the classed pipeline objects.
#'
#' @param v Any scalar value.
#' @return A character scalar.
#' @noRd
.hbr_fmt <- function(v) {
  if (is.null(v) || !length(v) || is.na(v[1])) "NA" else as.character(v[1])
}
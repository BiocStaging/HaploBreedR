# ============================================================================
# R/ld-regions.r
# LD-decay estimation, LD-window construction, and gene annotation of
# LD regions.
#
# ----------------------------------------------------------------------------
# LD decay
# ----------------------------------------------------------------------------

#' Compute LD-decay distance from PLINK data
#'
#' Runs PLINK `--r2` on the given dataset, bins SNP pairs by distance,
#' and returns the distance at which mean or median r-squared drops
#' below the requested threshold (interpolated between the two flanking
#' bins). Inputs are validated against missing or malformed PLINK files.
#'
#' All console output is gated by `verbose`. `TRUE` (default) prints a
#' single line reporting the estimated decay distance, the method used,
#' and the number of SNP pairs and bins; `FALSE` silences everything
#' except errors and warnings.
#'
#' @param plink_prefix PLINK file prefix or `.bed` path.
#' @param plink Path to the PLINK executable.
#' @param max_distance Maximum SNP-pair distance in bp.
#' @param bin_size Distance bin size in bp.
#' @param decay_threshold r-squared threshold for decay.
#' @param use_median Use median r-squared per bin instead of mean.
#' @param extra_plink_args Optional extra PLINK arguments.
#' @param verbose Print a one-line summary. `TRUE` (default) reports the
#'   decay distance, method, pair and bin counts; `FALSE` is silent.
#'
#' @return A structured list with `decay_distance_bp`,
#'   `decay_distance_kb`, `bin_summary`, and diagnostic fields.
#'
#' @noRd
hbr_calc_ld_decay <- function(plink_prefix,
                              plink = "plink",
                              max_distance = 500000,
                              bin_size = 10000,
                              decay_threshold = 0.20,
                              use_median = FALSE,
                              extra_plink_args = NULL,
                              verbose = TRUE) {

  if (!is.character(plink_prefix) || length(plink_prefix) != 1L ||
      is.na(plink_prefix) || !nzchar(plink_prefix))
    stop("plink_prefix must be a single PLINK prefix.", call. = FALSE)

  if (!is.character(plink) || length(plink) != 1L ||
      is.na(plink) || !nzchar(plink))
    stop("plink must be a valid PLINK executable path.", call. = FALSE)

  if (!is.numeric(max_distance) || length(max_distance) != 1L ||
      is.na(max_distance) || max_distance <= 0)
    stop("max_distance must be > 0.", call. = FALSE)

  if (!is.numeric(bin_size) || length(bin_size) != 1L ||
      is.na(bin_size) || bin_size <= 0)
    stop("bin_size must be > 0.", call. = FALSE)

  if (bin_size > max_distance)
    stop("bin_size cannot exceed max_distance.", call. = FALSE)

  if (!is.numeric(decay_threshold) || length(decay_threshold) != 1L ||
      is.na(decay_threshold) || decay_threshold <= 0 || decay_threshold >= 1)
    stop("decay_threshold must be between 0 and 1.", call. = FALSE)

  if (!is.logical(use_median) || length(use_median) != 1L || is.na(use_median))
    stop("use_median must be TRUE or FALSE.", call. = FALSE)

  if (!is.null(extra_plink_args) && !is.character(extra_plink_args))
    stop("extra_plink_args must be NULL or a character vector.",
         call. = FALSE)

  prefix <- if (grepl("\\.bed$", plink_prefix, ignore.case = TRUE))
    sub("\\.bed$", "", plink_prefix, ignore.case = TRUE)
  else plink_prefix

  bed_file <- paste0(prefix, ".bed")
  bim_file <- paste0(prefix, ".bim")
  fam_file <- paste0(prefix, ".fam")

  missing_files <- c(bed_file, bim_file, fam_file)[
    !file.exists(c(bed_file, bim_file, fam_file))]
  if (length(missing_files))
    stop("Incomplete PLINK dataset.\nMissing file(s):\n",
         paste(missing_files, collapse = "\n"), call. = FALSE)

  plink_executable <- Sys.which(plink)
  if (!nzchar(plink_executable)) {
    if (file.exists(plink))
      plink_executable <- normalizePath(plink, mustWork = TRUE)
    else
      stop("PLINK executable not found: ", plink, call. = FALSE)
  }

  temp_dir <- tempfile(pattern = "HaploBreedR_LD_")
  dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
  on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

  ld_prefix <- file.path(temp_dir, "ld_decay")

  max_distance_kb <- ceiling(max_distance / 1000)

  plink_args <- c("--bfile", prefix,
                  "--r2",
                  "--ld-window-kb", as.character(max_distance_kb),
                  "--ld-window", "999999",
                  "--ld-window-r2", "0",
                  "--out", ld_prefix)
  if (!is.null(extra_plink_args))
    plink_args <- c(plink_args, extra_plink_args)

  plink_output <- tryCatch(
    system2(command = plink_executable, args = plink_args,
            stdout = TRUE, stderr = TRUE),
    error = function(e)
      stop("Failed to execute PLINK: ", conditionMessage(e), call. = FALSE))

  ld_file <- paste0(ld_prefix, ".ld")

  if (!file.exists(ld_file))
    stop("PLINK LD calculation failed; .ld not created.\n\n",
         "PLINK output:\n", paste(plink_output, collapse = "\n"),
         call. = FALSE)
  if (file.info(ld_file)$size == 0)
    stop("PLINK generated an empty LD file.", call. = FALSE)

  ld <- data.table::fread(ld_file, header = TRUE, showProgress = FALSE)
  if (!nrow(ld))
    stop("PLINK LD file contains no SNP-pair records.", call. = FALSE)

  required_columns <- c("CHR_A", "BP_A", "SNP_A", "CHR_B", "BP_B",
                        "SNP_B", "R2")
  if (!all(required_columns %in% names(ld))) {
    if (ncol(ld) >= 7L) {
      ld <- ld[, seq_len(7L)]
      data.table::setnames(ld, required_columns)
    } else {
      stop("Unable to identify PLINK LD columns.", call. = FALSE)
    }
  }

  ld <- ld[, .(CHR_A = as.character(CHR_A),
               BP_A  = suppressWarnings(as.numeric(BP_A)),
               SNP_A = as.character(SNP_A),
               CHR_B = as.character(CHR_B),
               BP_B  = suppressWarnings(as.numeric(BP_B)),
               SNP_B = as.character(SNP_B),
               R2    = suppressWarnings(as.numeric(R2)))]

  ld <- ld[!is.na(CHR_A) & !is.na(CHR_B) &
             !is.na(BP_A) & !is.na(BP_B) & !is.na(R2)]
  if (!nrow(ld))
    stop("No valid SNP-pair records after filtering.", call. = FALSE)

  invalid_r2 <- ld[R2 < 0 | R2 > 1]
  if (nrow(invalid_r2) > 0) {
    warning(nrow(invalid_r2), " SNP pairs have R2 outside [0,1] and will ",
            "be removed.", call. = FALSE)
    ld <- ld[R2 >= 0 & R2 <= 1]
  }
  if (!nrow(ld))
    stop("No valid R2 values remain.", call. = FALSE)

  ld <- ld[CHR_A == CHR_B & BP_A != BP_B]
  if (!nrow(ld))
    stop("No valid intra-chromosomal SNP pairs found.", call. = FALSE)

  ld[, Distance := abs(BP_B - BP_A)]
  ld <- ld[Distance > 0 & Distance <= max_distance]
  if (!nrow(ld))
    stop("No SNP pairs remain within max_distance.", call. = FALSE)

  ld[, Distance_Bin := floor(Distance / bin_size) * bin_size]
  ld[Distance_Bin < bin_size, Distance_Bin := bin_size]

  bin_summary <- ld[, .(
    Mean_R2       = mean(R2, na.rm = TRUE),
    Median_R2     = stats::median(R2, na.rm = TRUE),
    N_Pairs       = .N,
    Min_Distance  = min(Distance, na.rm = TRUE),
    Max_Distance  = max(Distance, na.rm = TRUE)),
    by = Distance_Bin]
  data.table::setorder(bin_summary, Distance_Bin)

  metric_col <- if (use_median) "Median_R2" else "Mean_R2"
  decay_bin <- bin_summary[get(metric_col) < decay_threshold]

  if (nrow(decay_bin) > 0) {
    first_idx <- which(bin_summary[[metric_col]] < decay_threshold)[1]

    if (first_idx == 1L) {
      decay_distance <- bin_summary$Distance_Bin[1]
      method_used <- if (use_median) "median_R2_threshold"
                     else "mean_R2_threshold"
    } else {
      x1 <- bin_summary$Distance_Bin[first_idx - 1]
      x2 <- bin_summary$Distance_Bin[first_idx]
      y1 <- bin_summary[[metric_col]][first_idx - 1]
      y2 <- bin_summary[[metric_col]][first_idx]

      if (is.finite(y1) && is.finite(y2) && y1 != y2) {
        decay_distance <- x1 + (decay_threshold - y1) * (x2 - x1) / (y2 - y1)
        method_used <- if (use_median) "median_R2_interpolated"
                       else "mean_R2_interpolated"
      } else {
        decay_distance <- x2
        method_used <- if (use_median) "median_R2_first_bin"
                       else "mean_R2_first_bin"
      }
    }
  } else {
    decay_distance <- NA_real_
    method_used <- "threshold_not_reached"
    warning("LD did not fall below r2 = ", decay_threshold,
            " within max_distance = ",
            format(max_distance, big.mark = ",", scientific = FALSE),
            " bp. LD decay cannot be estimated.", call. = FALSE)
  }

  if (!is.na(decay_distance))
    decay_distance <- max(bin_size, min(decay_distance, max_distance))
  decay_distance <- as.numeric(decay_distance)

  total_pairs           <- nrow(ld)
  n_bins                <- nrow(bin_summary)
  max_observed_distance <- max(ld$Distance, na.rm = TRUE)
  mean_r2_overall       <- mean(ld$R2, na.rm = TRUE)
  median_r2_overall     <- stats::median(ld$R2, na.rm = TRUE)

  if (isTRUE(verbose))
    message(sprintf(
      "LD decay: %.1f kb at r2 < %.2f (%s) | %s SNP pairs, %d bins",
      decay_distance / 1000, decay_threshold, method_used,
      format(total_pairs, big.mark = ","), n_bins))

  list(decay_distance_bp     = decay_distance,
       decay_distance_kb     = decay_distance / 1000,
       bin_summary           = bin_summary,
       threshold             = decay_threshold,
       max_distance          = max_distance,
       bin_size              = bin_size,
       use_median            = use_median,
       method                = method_used,
       n_pairs               = total_pairs,
       n_bins                = n_bins,
       max_observed_distance = max_observed_distance,
       mean_r2               = mean_r2_overall,
       median_r2             = median_r2_overall,
       plink_prefix          = normalizePath(prefix, winslash = "/",
                                             mustWork = FALSE))
}


# ----------------------------------------------------------------------------
# LD-window construction
# ----------------------------------------------------------------------------

#' Create LD-defined genomic regions from GWAS SNPs
#'
#' Expands each GWAS SNP into a symmetric window of ±`ld_decay_distance`
#' base pairs, merges overlapping (and optionally adjacent) windows into
#' LD regions, and preserves full provenance from the original GWAS rows.
#'
#' @param gwas_snps Data frame with SNP, chromosome, and position.
#' @param ld_decay_distance Window half-width in bp.
#' @param snp_col,chr_col,pos_col Column names in `gwas_snps`.
#' @param merge_overlapping Merge overlapping windows.
#' @param merge_adjacent Also merge touching windows.
#' @param verbose Print progress.
#'
#' @return A `data.table` of merged regions with columns
#'   `LD_Region_ID`, `Chromosome`, `Region_Start`, `Region_End`,
#'   `Region_Length`, `GWAS_SNP_Count`, `GWAS_SNPs`, `GWAS_Positions`,
#'   `LD_Decay_bp`, `LD_Decay_kb`, `Original_Window_IDs`,
#'   `Original_GWAS_Rows`.
#'
#' @noRd
hbr_create_ld_windows <- function(gwas_snps,
                                  ld_decay_distance,
                                  snp_col = "SNP",
                                  chr_col = "Chromosome",
                                  pos_col = "Position",
                                  merge_overlapping = TRUE,
                                  merge_adjacent = TRUE,
                                  verbose = TRUE) {

  if (is.null(gwas_snps) || nrow(gwas_snps) == 0)
    stop("gwas_snps is empty.", call. = FALSE)

  gwas_snps <- data.table::as.data.table(data.table::copy(gwas_snps))

  required_cols <- c(snp_col, chr_col, pos_col)
  missing_cols  <- setdiff(required_cols, names(gwas_snps))
  if (length(missing_cols))
    stop("gwas_snps is missing required columns: ",
         paste(missing_cols, collapse = ", "), call. = FALSE)

  if (length(ld_decay_distance) != 1L ||
      !is.numeric(ld_decay_distance) ||
      is.na(ld_decay_distance) ||
      !is.finite(ld_decay_distance) ||
      ld_decay_distance <= 0)
    stop("`ld_decay_distance` must be a single positive numeric value.",
         call. = FALSE)

  ld_decay_distance <- as.numeric(ld_decay_distance)
  n_input <- nrow(gwas_snps)

  gwas_snps[, Original_GWAS_Row := seq_len(.N)]
  gwas_snps[, .SNP_ID_INTERNAL :=
              trimws(as.character(get(snp_col)))]
  gwas_snps[, .CHR_INTERNAL :=
              trimws(as.character(get(chr_col)))]
  gwas_snps[, .POS_INTERNAL :=
              suppressWarnings(as.numeric(get(pos_col)))]

  valid_rows <- (!is.na(gwas_snps$.SNP_ID_INTERNAL) &
                   nzchar(gwas_snps$.SNP_ID_INTERNAL) &
                   !is.na(gwas_snps$.CHR_INTERNAL) &
                   nzchar(gwas_snps$.CHR_INTERNAL) &
                   !is.na(gwas_snps$.POS_INTERNAL) &
                   is.finite(gwas_snps$.POS_INTERNAL) &
                   gwas_snps$.POS_INTERNAL > 0)
  n_invalid <- sum(!valid_rows)
  if (n_invalid > 0) {
    if (verbose)
      message("Invalid GWAS records removed: ", n_invalid)
    gwas_snps <- gwas_snps[valid_rows]
  }
  n_valid <- nrow(gwas_snps)
  if (n_valid == 0)
    stop("No valid GWAS SNP records remain after validation.",
         call. = FALSE)

  normalize_chr <- function(x) {
    x <- trimws(as.character(x))
    x <- sub("^[Cc][Hh][Rr]", "", x)
    numeric_x <- suppressWarnings(as.integer(x))
    idx <- !is.na(numeric_x)
    x[idx] <- as.character(numeric_x[idx])
    x
  }
  gwas_snps[, .CHR_INTERNAL := normalize_chr(.CHR_INTERNAL)]

  inconsistent_snps <- gwas_snps[
    , .(n_positions = data.table::uniqueN(.POS_INTERNAL),
        n_chr       = data.table::uniqueN(.CHR_INTERNAL)),
    by = .SNP_ID_INTERNAL][n_positions > 1 | n_chr > 1]
  if (nrow(inconsistent_snps) > 0)
    stop("QC failure: ", nrow(inconsistent_snps),
         " SNP IDs have multiple chromosome/position assignments.",
         call. = FALSE)

  unique_snp_count      <- data.table::uniqueN(gwas_snps$.SNP_ID_INTERNAL)
  duplicate_snp_records <- n_valid - unique_snp_count

  if (verbose) {
    message("GWAS SNP records supplied: ", n_input)
    message("Valid GWAS SNP records: ", n_valid)
    message("Unique GWAS SNPs: ", unique_snp_count)
    message("Repeated SNP records: ", duplicate_snp_records)
    if (n_invalid > 0)
      message("Invalid records removed: ", n_invalid)
  }

  unique_gwas_snps <- gwas_snps[
    , .(GWAS_SNP   = .SNP_ID_INTERNAL,
        Chromosome = .CHR_INTERNAL,
        Position   = .POS_INTERNAL,
        Original_GWAS_Rows = list(unique(Original_GWAS_Row))),
    by = .(.SNP_ID_INTERNAL, .CHR_INTERNAL, .POS_INTERNAL)]
  unique_gwas_snps[, c(".SNP_ID_INTERNAL", ".CHR_INTERNAL",
                       ".POS_INTERNAL") := NULL]

  if (nrow(unique_gwas_snps) != unique_snp_count)
    stop("QC failure: expected ", unique_snp_count,
         " unique SNPs but obtained ", nrow(unique_gwas_snps), ".",
         call. = FALSE)

  windows <- data.table::copy(unique_gwas_snps)
  windows[, Region_Start := pmax(1, Position - ld_decay_distance)]
  windows[, Region_End   := Position + ld_decay_distance]
  windows[, LD_Decay_bp  := ld_decay_distance]
  windows[, LD_Decay_kb  := ld_decay_distance / 1000]
  windows[, Region_Length := Region_End - Region_Start + 1]
  windows[, Original_Window_ID := paste0("GWAS_Window_", seq_len(.N))]

  if (nrow(windows) != unique_snp_count)
    stop("QC failure: ", nrow(windows), " windows for ",
         unique_snp_count, " unique SNPs.", call. = FALSE)

  input_unique_snps <- sort(unique(as.character(gwas_snps$.SNP_ID_INTERNAL)))
  window_snps       <- sort(unique(as.character(windows$GWAS_SNP)))
  if (!identical(input_unique_snps, window_snps))
    stop("QC failure: SNP identity was not preserved.", call. = FALSE)

  if (!merge_overlapping) {
    windows[, LD_Region_ID := paste0("GWAS_LD_Region_", seq_len(.N))]
    windows[, GWAS_SNP_Count := 1L]
    windows[, GWAS_SNPs := lapply(GWAS_SNP, function(x) x)]
    windows[, GWAS_Positions := lapply(Position, function(x) x)]
    windows[, Original_Window_IDs := lapply(Original_Window_ID,
                                            function(x) x)]
    data.table::setcolorder(windows,
                            c("LD_Region_ID", "Chromosome", "Region_Start",
                              "Region_End", "Region_Length", "GWAS_SNP_Count",
                              "GWAS_SNPs", "GWAS_Positions",
                              "LD_Decay_bp", "LD_Decay_kb",
                              "Original_Window_IDs", "Original_GWAS_Rows"))
    return(windows[])
  }

  chr_values  <- unique(windows$Chromosome)
  chr_numeric <- suppressWarnings(as.numeric(chr_values))
  chr_order   <- if (all(!is.na(chr_numeric)))
    chr_values[order(chr_numeric)] else sort(chr_values)
  windows[, Chromosome := factor(Chromosome, levels = chr_order)]
  data.table::setorder(windows, Chromosome, Region_Start, Region_End,
                       Position, GWAS_SNP)
  windows[, Chromosome := as.character(Chromosome)]

  windows[, Cum_End := cummax(Region_End), by = Chromosome]
  windows[, Previous_Cum_End := data.table::shift(Cum_End, type = "lag"),
          by = Chromosome]

  if (merge_adjacent) {
    windows[, New_Region := is.na(Previous_Cum_End) |
               Region_Start > Previous_Cum_End + 1]
  } else {
    windows[, New_Region := is.na(Previous_Cum_End) |
               Region_Start > Previous_Cum_End]
  }

  windows[, Region_Group := cumsum(New_Region), by = Chromosome]
  if (any(is.na(windows$Region_Group)))
    stop("QC failure: some windows were not assigned to an LD region.",
         call. = FALSE)

  regions <- windows[
    , .(Region_Start = min(Region_Start),
        Region_End   = max(Region_End),
        GWAS_SNPs    = list(unique(as.character(GWAS_SNP))),
        GWAS_Positions = list(unique(as.numeric(Position))),
        Original_Window_IDs = list(unique(as.character(Original_Window_ID))),
        Original_GWAS_Rows = list(sort(unique(unlist(
          Original_GWAS_Rows, use.names = FALSE))))),
    by = .(Chromosome, Region_Group)]

  regions[, Region_Length := Region_End - Region_Start + 1]
  regions[, GWAS_SNP_Count := lengths(GWAS_SNPs)]
  regions[, LD_Decay_bp := ld_decay_distance]
  regions[, LD_Decay_kb := ld_decay_distance / 1000]
  data.table::setorder(regions, Chromosome, Region_Start, Region_End)
  regions[, LD_Region_ID := paste0("GWAS_LD_Region_", seq_len(.N))]
  regions[, Region_Group := NULL]
  data.table::setcolorder(regions,
                          c("LD_Region_ID", "Chromosome", "Region_Start",
                            "Region_End", "Region_Length", "GWAS_SNP_Count",
                            "GWAS_SNPs", "GWAS_Positions",
                            "LD_Decay_bp", "LD_Decay_kb",
                            "Original_Window_IDs", "Original_GWAS_Rows"))

  output_unique_snps <- sort(unique(unlist(regions$GWAS_SNPs,
                                            use.names = FALSE)))
  if (!identical(input_unique_snps, output_unique_snps))
    stop("QC failure: unique GWAS SNP preservation failed.", call. = FALSE)

  if (sum(regions$GWAS_SNP_Count) != unique_snp_count)
    stop("QC failure: region SNP counts do not sum to unique SNP count.",
         call. = FALSE)

  output_original_rows <- sort(unique(unlist(regions$Original_GWAS_Rows,
                                               use.names = FALSE)))
  input_original_rows  <- sort(gwas_snps$Original_GWAS_Row)
  if (!identical(input_original_rows, output_original_rows))
    stop("QC failure: original GWAS row provenance was not preserved.",
         call. = FALSE)

  if (verbose) {
    message("")
    message("LD-region construction completed.")
    message("----------------------------------")
    message("GWAS SNP records supplied: ", n_input)
    message("Valid GWAS SNP records: ", n_valid)
    message("Unique GWAS SNPs: ", unique_snp_count)
    message("SNP-centered windows: ", nrow(windows))
    message("Merged LD regions: ", nrow(regions))
    message("LD-decay distance: ",
            format(ld_decay_distance, big.mark = ",", scientific = FALSE),
            " bp (", round(ld_decay_distance / 1000, 2), " kb)")
    message("SNP preservation QC: PASSED")
    message("Original-row provenance QC: PASSED")
    message("")
  }

  regions[]
}


# ----------------------------------------------------------------------------
# Region-gene annotation
# ----------------------------------------------------------------------------

#' Annotate LD regions with genes
#'
#' Overlaps a set of LD regions with a gene table and labels each
#' region-gene pair with a `Gene_Relationship` (`GENE_CONTAINED`,
#' `OVERLAPS_LEFT`, `OVERLAPS_RIGHT`, `CONTAINS_REGION`, or
#' `OVERLAPS`). Region list columns (e.g. `GWAS_SNPs`) are preserved
#' through the overlap. When `keep_intergenic = TRUE`, regions with no
#' overlapping gene are retained as `INTERGENIC`.
#'
#' @param regions A `data.table` of LD regions.
#' @param genes Gene table with ID, chromosome, start, end, and strand
#'   columns (column names auto-detected).
#' @param chr_col,start_col,end_col Region coordinate column names.
#' @param keep_intergenic Keep regions with no overlapping gene.
#' @param verbose Print progress.
#' @examples
#' regions <- data.table::data.table(
#'   LD_Region_ID = "R1",
#'   Chromosome   = "1",
#'   Region_Start = 1000L,
#'   Region_End   = 3000L)
#' genes <- data.table::data.table(
#'   Gene_ID         = "gene1",
#'   Gene_Chromosome = "1",
#'   Gene_Start      = 1000L,
#'   Gene_End        = 2500L,
#'   Gene_Strand     = "+")
#' hbr_annotate_regions_with_genes(regions, genes, verbose = FALSE)
#' @return A `data.table` of region-gene overlaps, or `NULL`.
#'
#' @export
hbr_annotate_regions_with_genes <- function(
    regions, genes,
    chr_col = "Chromosome",
    start_col = "Region_Start",
    end_col = "Region_End",
    keep_intergenic = TRUE,
    verbose = TRUE) {

  if (is.null(regions) || !nrow(regions)) return(NULL)
  if (is.null(genes) || !nrow(genes)) return(NULL)

  regions <- data.table::copy(data.table::as.data.table(regions))
  genes   <- data.table::copy(data.table::as.data.table(genes))

  region_id_col <- if ("LD_Region_ID" %in% names(regions))
    "LD_Region_ID" else NULL

  original_region_metadata <- data.table::copy(regions)
  list_cols_regions <- names(regions)[vapply(regions, is.list, logical(1))]

  if (verbose && length(list_cols_regions))
    cat("Preserving list columns separately:",
        paste(list_cols_regions, collapse = ", "), "\n")

  regions_atomic <- data.table::copy(regions)
  if (length(list_cols_regions))
    regions_atomic <- regions_atomic[, !list_cols_regions, with = FALSE]

  id_candidates <- intersect(
    c("Gene_ID", "gene_id", "GeneID", "geneID", "ID",
      "Gene", "RAPDB_Gene_ID", "MSU_Gene_ID"), names(genes))
  if (!length(id_candidates))
    stop("No gene ID column found. Available: ",
         paste(names(genes), collapse = ", "), call. = FALSE)

  chr_candidates <- intersect(
    c("Gene_Chromosome", "Chromosomes", "Chr", "chr",
      "Chromosome", "CHR"), names(genes))
  if (!length(chr_candidates))
    stop("No chromosome column found.", call. = FALSE)

  start_candidates <- intersect(
    c("Gene_Start", "Gene_Start2", "Start", "start"), names(genes))
  if (!length(start_candidates))
    stop("No start column found.", call. = FALSE)

  end_candidates <- intersect(
    c("Gene_End", "Gene_End2", "End", "end"), names(genes))
  if (!length(end_candidates))
    stop("No end column found.", call. = FALSE)

  data.table::setnames(genes, id_candidates[1],    "Gene_ID",
                       skip_absent = TRUE)
  data.table::setnames(genes, chr_candidates[1],   "Gene_Chromosome",
                       skip_absent = TRUE)
  data.table::setnames(genes, start_candidates[1], "Gene_Start",
                       skip_absent = TRUE)
  data.table::setnames(genes, end_candidates[1],   "Gene_End",
                       skip_absent = TRUE)

  if (!"Gene_Name" %in% names(genes)) {
    if ("gene_name" %in% names(genes))
      data.table::setnames(genes, "gene_name", "Gene_Name")
    else genes[, Gene_Name := Gene_ID]
  }

  if (!"Gene_Function" %in% names(genes)) {
    if ("gene_function" %in% names(genes))
      data.table::setnames(genes, "gene_function", "Gene_Function")
    else if ("Description" %in% names(genes))
      data.table::setnames(genes, "Description", "Gene_Function")
    else genes[, Gene_Function := NA_character_]
  }

  if (!"Gene_Strand" %in% names(genes)) {
    if ("Strand" %in% names(genes))
      data.table::setnames(genes, "Strand", "Gene_Strand")
    else if ("strand" %in% names(genes))
      data.table::setnames(genes, "strand", "Gene_Strand")
    else genes[, Gene_Strand := NA_character_]
  }

  list_cols_genes <- names(genes)[vapply(genes, is.list, logical(1))]
  if (length(list_cols_genes)) {
    if (verbose)
      cat("Removing list columns from genes:",
          paste(list_cols_genes, collapse = ", "), "\n")
    genes <- genes[, !list_cols_genes, with = FALSE]
  }

  required_region_cols <- c(chr_col, start_col, end_col)
  missing_region <- setdiff(required_region_cols, names(regions_atomic))
  if (length(missing_region))
    stop("regions missing: ",
         paste(missing_region, collapse = ", "), call. = FALSE)

  regions_atomic[, .Region_Chr   := as.character(get(chr_col))]
  regions_atomic[, .Region_Start :=
                   suppressWarnings(as.integer(get(start_col)))]
  regions_atomic[, .Region_End   :=
                   suppressWarnings(as.integer(get(end_col)))]

  genes[, .Gene_Chr   := as.character(Gene_Chromosome)]
  genes[, .Gene_Start := suppressWarnings(as.integer(Gene_Start))]
  genes[, .Gene_End   := suppressWarnings(as.integer(Gene_End))]

  normalize_chr <- function(x) {
    x <- trimws(as.character(x))
    x <- sub("^[Cc][Hh][Rr]", "", x)
    numeric_x <- suppressWarnings(as.integer(x))
    idx <- !is.na(numeric_x)
    x[idx] <- as.character(numeric_x[idx])
    x
  }

  regions_atomic[, .Region_Chr := normalize_chr(.Region_Chr)]
  genes[,          .Gene_Chr   := normalize_chr(.Gene_Chr)]

  regions_atomic <- regions_atomic[
    !is.na(.Region_Chr) & !is.na(.Region_Start) &
      !is.na(.Region_End) & .Region_Start <= .Region_End]
  genes <- genes[
    !is.na(.Gene_Chr) & !is.na(.Gene_Start) &
      !is.na(.Gene_End) & .Gene_Start <= .Gene_End]

  if (!nrow(regions_atomic)) {
    warning("No valid regions"); return(NULL)
  }
  if (!nrow(genes)) {
    warning("No valid genes"); return(NULL)
  }

  common_chr <- intersect(unique(regions_atomic$.Region_Chr),
                          unique(genes$.Gene_Chr))
  if (!length(common_chr)) {
    warning("No common chromosomes"); return(NULL)
  }

  if (verbose)
    cat("Common chromosomes:",
        paste(sort(common_chr), collapse = ", "), "\n")

  regions_atomic <- regions_atomic[.Region_Chr %in% common_chr]
  genes          <- genes         [.Gene_Chr   %in% common_chr]

  regions_atomic[, `:=`(.join_start = .Region_Start,
                         .join_end   = .Region_End)]
  genes[,          `:=`(.join_start = .Gene_Start,
                         .join_end   = .Gene_End)]

  atomic_check <- function(x, col) {
    if (!col %in% names(x)) return(x)
    if (is.list(x[[col]]))
      stop("Internal error: overlap column '", col,
           "' is still a list column.", call. = FALSE)
    x
  }
  for (col in c(".Region_Chr", ".join_start", ".join_end"))
    regions_atomic <- atomic_check(regions_atomic, col)
  for (col in c(".Gene_Chr", ".join_start", ".join_end"))
    genes <- atomic_check(genes, col)

  data.table::setkey(genes, .Gene_Chr, .join_start, .join_end)

  if (verbose) cat("Finding genomic region x gene overlaps...\n")

  overlaps <- data.table::foverlaps(
    regions_atomic, genes,
    by.x = c(".Region_Chr", ".join_start", ".join_end"),
    by.y = c(".Gene_Chr",   ".join_start", ".join_end"),
    type = "any", nomatch = 0L)

  if (is.null(overlaps) || !nrow(overlaps)) {
    if (!keep_intergenic) return(NULL)

    result <- data.table::copy(regions_atomic)
    result[, `:=`(
      Gene_ID           = NA_character_,
      Gene_Name         = NA_character_,
      Gene_Function     = NA_character_,
      Gene_Chromosome   = NA_character_,
      Gene_Start        = NA_integer_,
      Gene_End          = NA_integer_,
      Gene_Strand       = NA_character_,
      Gene_Relationship = "INTERGENIC",
      Annotation_Status = "Intergenic")]

    helper_cols <- intersect(
      c(".Region_Chr", ".Region_Start", ".Region_End",
        ".Gene_Chr", ".Gene_Start", ".Gene_End",
        ".join_start", ".join_end"), names(result))
    if (length(helper_cols)) result[, (helper_cols) := NULL]

    if (!is.null(region_id_col) &&
        region_id_col %in% names(original_region_metadata)) {
      metadata_cols <- intersect(names(original_region_metadata),
                                 list_cols_regions)
      if (length(metadata_cols)) {
        metadata <- original_region_metadata[
          , c(region_id_col, metadata_cols), with = FALSE]
        result <- merge(result, metadata, by = region_id_col,
                        all.x = TRUE, sort = FALSE)
      }
    }
    return(result[])
  }

  if (!"Gene_ID" %in% names(overlaps) && "i.Gene_ID" %in% names(overlaps))
    data.table::setnames(overlaps, "i.Gene_ID", "Gene_ID")

  if (!"Region_Start" %in% names(overlaps))
    overlaps[, Region_Start := .Region_Start]
  if (!"Region_End" %in% names(overlaps))
    overlaps[, Region_End := .Region_End]
  if (!"Chromosome" %in% names(overlaps))
    overlaps[, Chromosome := .Region_Chr]
  if (!is.null(region_id_col) && !"LD_Region_ID" %in% names(overlaps))
    overlaps[, LD_Region_ID := get(region_id_col)]

  overlaps[, Gene_Relationship := data.table::fcase(
    Gene_Start >= Region_Start & Gene_End <= Region_End,
      "GENE_CONTAINED",
    Gene_Start < Region_Start & Gene_End >= Region_Start &
      Gene_End <= Region_End,
      "OVERLAPS_LEFT",
    Gene_Start >= Region_Start & Gene_Start <= Region_End &
      Gene_End > Region_End,
      "OVERLAPS_RIGHT",
    Gene_Start < Region_Start & Gene_End > Region_End,
      "CONTAINS_REGION",
    default = "OVERLAPS")]

  overlaps[, Annotation_Status := "Annotated"]

  helper_cols <- intersect(
    c(".Region_Chr", ".Region_Start", ".Region_End",
      ".Gene_Chr", ".Gene_Start", ".Gene_End",
      ".join_start", ".join_end",
      "i..join_start", "i..join_end"), names(overlaps))
  if (length(helper_cols)) {
    if (verbose)
      cat("Removing helper columns:",
          paste(helper_cols, collapse = ", "), "\n")
    overlaps[, (helper_cols) := NULL]
  }

  dup_names <- duplicated(names(overlaps))
  if (any(dup_names))
    overlaps <- overlaps[, !dup_names, with = FALSE]

  if (!is.null(region_id_col) &&
      region_id_col %in% names(original_region_metadata)) {
    metadata_cols <- intersect(list_cols_regions,
                               names(original_region_metadata))
    if (length(metadata_cols)) {
      metadata <- original_region_metadata[
        , c(region_id_col, metadata_cols), with = FALSE]
      metadata <- metadata[, .SD[1], by = region_id_col]
      idx <- match(overlaps[[region_id_col]], metadata[[region_id_col]])
      for (col in metadata_cols)
        overlaps[, (col) := metadata[[col]][idx]]
    }
  }

  standard_cols <- c(
    "LD_Region_ID", "Chromosome", "Region_Start", "Region_End",
    "Region_Length", "GWAS_SNP_Count", "GWAS_SNPs", "GWAS_Positions",
    "LD_Decay_bp", "LD_Decay_kb",
    "Original_Window_IDs", "Original_GWAS_Rows",
    "Gene_ID", "Gene_Name", "Gene_Function",
    "Gene_Chromosome", "Gene_Start", "Gene_End", "Gene_Strand",
    "Gene_Relationship", "Annotation_Status")

  existing_standard <- intersect(standard_cols, names(overlaps))
  additional_cols   <- setdiff(names(overlaps), existing_standard)
  data.table::setcolorder(overlaps, c(existing_standard, additional_cols))

  dedup_cols <- intersect(
    c("LD_Region_ID", "Gene_ID", "Chromosome",
      "Region_Start", "Region_End"), names(overlaps))
  if (length(dedup_cols))
    overlaps <- overlaps[!duplicated(overlaps[, ..dedup_cols])]

  if (verbose) {
    cat("Annotated",
        data.table::uniqueN(overlaps$LD_Region_ID),
        "regions with",
        data.table::uniqueN(overlaps$Gene_ID),
        "genes.\n")
  }

  overlaps[]
}

#' Parse gene features from a GFF3 into a data.table
#'
#' Reads a GFF3 file, extracts `gene`-type features, decodes the
#' `ID`/`gene_id`/`locus_tag`, `Name`, and `Note`/`description`/
#' `product` attributes, standardizes chromosome names via
#' [hbr_standardize_chr()], and returns a deduplicated gene table.
#'
#' @param gff3 Path to a GFF3 file.
#' @param verbose Print progress.
#' @examples
#' regions <- data.table::data.table(
#'   LD_Region_ID = "R1",
#'   Chromosome   = "1",
#'   Region_Start = 1000L,
#'   Region_End   = 3000L)
#' genes <- data.table::data.table(
#'   Gene_ID         = "gene1",
#'   Gene_Chromosome = "1",
#'   Gene_Start      = 1000L,
#'   Gene_End        = 2500L,
#'   Gene_Strand     = "+")
#' hbr_annotate_regions_with_genes(regions, genes, verbose = FALSE)
#' @return A `data.table` with `Gene_ID`, `Gene_Name`,
#'   `Gene_Function`, `Gene_Chromosome`, `Gene_Start`, `Gene_End`,
#'   `Gene_Strand`. Empty if no gene features are present.
#'
#' @export
hbr_read_gff3_genes <- function(gff3, verbose = TRUE) {

  if (missing(gff3) || is.null(gff3) || !is.character(gff3) ||
      length(gff3) != 1L || !nzchar(gff3))
    stop("A valid GFF3 file path must be provided.", call. = FALSE)
  if (!file.exists(gff3))
    stop("GFF3 file not found: ", gff3, call. = FALSE)

  if (verbose) message("Reading GFF3 gene annotation...")

  gff <- data.table::fread(
    file = gff3, sep = "\t", header = FALSE, comment.char = "#",
    quote = "", fill = TRUE, showProgress = verbose)

  if (ncol(gff) < 9L)
    stop("Invalid GFF3 file: expected at least 9 columns, found ",
         ncol(gff), ".", call. = FALSE)

  data.table::setnames(
    gff, old = names(gff)[seq_len(9L)],
    new = c("Chromosome", "Source", "Type", "Start", "End",
            "Score", "Strand", "Phase", "Attributes"))

  genes <- gff[
    !is.na(Type) & tolower(trimws(Type)) == "gene",
    .(Gene_Chromosome = as.character(Chromosome),
      Gene_Start      = suppressWarnings(as.integer(Start)),
      Gene_End        = suppressWarnings(as.integer(End)),
      Gene_Strand     = as.character(Strand),
      Attributes      = as.character(Attributes))]
  rm(gff)

  if (!nrow(genes)) {
    warning("No gene features found in GFF3.", call. = FALSE)
    return(data.table::data.table(
      Gene_ID = character(), Gene_Name = character(),
      Gene_Function = character(), Gene_Chromosome = character(),
      Gene_Start = integer(), Gene_End = integer(),
      Gene_Strand = character()))
  }

  genes <- genes[!is.na(Gene_Start) & !is.na(Gene_End) &
                   Gene_Start > 0L & Gene_End > 0L]
  if (!nrow(genes)) {
    warning("No valid gene coordinates found in GFF3.", call. = FALSE)
    return(data.table::data.table(
      Gene_ID = character(), Gene_Name = character(),
      Gene_Function = character(), Gene_Chromosome = character(),
      Gene_Start = integer(), Gene_End = integer(),
      Gene_Strand = character()))
  }

  attrs <- genes$Attributes

  extract_attr <- function(x, keys) {
    pattern <- paste0("(?:^|;)", "(?:", paste(keys, collapse = "|"),
                      ")=", "([^;]*)")
    value <- sub(pattern, "\\1", x, perl = TRUE)
    value[!grepl(pattern, x, perl = TRUE)] <- NA_character_
    value <- tryCatch(utils::URLdecode(value), error = function(e) value)
    value
  }

  gene_id       <- extract_attr(attrs, c("ID", "gene_id", "locus_tag"))
  gene_name     <- extract_attr(attrs, c("Name", "gene_name"))
  gene_name[is.na(gene_name) | !nzchar(gene_name)] <-
    gene_id[is.na(gene_name) | !nzchar(gene_name)]
  gene_function <- extract_attr(attrs, c("Note", "description", "product"))

  genes_dt <- data.table::data.table(
    Gene_ID         = as.character(gene_id),
    Gene_Name       = as.character(gene_name),
    Gene_Function   = as.character(gene_function),
    Gene_Chromosome = as.character(genes$Gene_Chromosome),
    Gene_Start      = genes$Gene_Start,
    Gene_End        = genes$Gene_End,
    Gene_Strand     = as.character(genes$Gene_Strand))

  genes_dt[, Gene_Chromosome := hbr_standardize_chr(Gene_Chromosome)]

  data.table::setkey(genes_dt, Gene_ID, Gene_Chromosome, Gene_Start, Gene_End)
  genes_dt <- unique(genes_dt)
  data.table::setkey(genes_dt, NULL)

  if (verbose) {
    message("Read ", format(nrow(genes_dt), big.mark = ","),
            " unique gene features from GFF3.")
    message("Chromosomes detected: ",
            paste(sort(unique(genes_dt$Gene_Chromosome)), collapse = ", "))
  }

  genes_dt[]
}

#' Annotate LD regions with MSU7 genes
#'
#' Parses an MSU7 GFF3 with [hbr_read_gff3_genes()] and overlaps the
#' resulting gene table with the supplied LD regions using
#' [hbr_annotate_regions_with_genes()]. The output is tagged with
#' `Annotation_Source = "MSU"` and `Annotation_Build = "MSU7"`.
#'
#' @inheritParams hbr_annotate_regions_with_genes
#' @param gff3 Path to an MSU7 GFF3 file.
#' @return A `data.table` of annotated regions, or `NULL`.
#' @examples
#' toy <- system.file("extdata", "toy", package = "HaploBreedR")
#' regions <- data.table::data.table(
#'   LD_Region_ID = "R1", Chromosome = "1",
#'   Region_Start = 1000L, Region_End = 3000L)
#' hbr_annotate_regions_with_genes_msu(
#'   regions, file.path(toy, "toy.gff3"), verbose = FALSE)
#'
#' @noRd
hbr_annotate_regions_with_genes_msu <- function(
    regions, gff3,
    chr_col = "Chromosome",
    start_col = "Region_Start",
    end_col = "Region_End",
    keep_intergenic = TRUE,
    verbose = TRUE) {

  if (verbose) cat("Parsing MSU GFF3...\n")
  genes <- hbr_read_gff3_genes(gff3 = gff3, verbose = verbose)
  if (is.null(genes) || !nrow(genes)) return(NULL)

  annotated <- hbr_annotate_regions_with_genes(
    regions = regions, genes = genes,
    chr_col = chr_col, start_col = start_col, end_col = end_col,
    keep_intergenic = keep_intergenic, verbose = verbose)

  if (!is.null(annotated) && nrow(annotated) > 0)
    annotated[, `:=`(Annotation_Source = "MSU", Annotation_Build = "MSU7")]
  annotated
}

#' Annotate LD regions with RAP-DB genes
#'
#' Parses an IRGSP-1.0 GFF3 with [hbr_read_gff3_genes()] and overlaps
#' the resulting gene table with the supplied LD regions using
#' [hbr_annotate_regions_with_genes()]. The output is tagged with
#' `Annotation_Source = "RAPDB"` and `Annotation_Build = "IRGSP-1.0"`.
#'
#' @inheritParams hbr_annotate_regions_with_genes
#' @param gff3 Path to an IRGSP-1.0 GFF3 file.
#' @return A `data.table` of annotated regions, or `NULL`.
#' @examples
#' toy <- system.file("extdata", "toy", package = "HaploBreedR")
#' regions <- data.table::data.table(
#'   LD_Region_ID = "R1", Chromosome = "1",
#'   Region_Start = 1000L, Region_End = 3000L)
#' hbr_annotate_regions_with_genes_rapdb(
#'   regions, file.path(toy, "toy.gff3"), verbose = FALSE)
#'
#' @noRd
hbr_annotate_regions_with_genes_rapdb <- function(
    regions, gff3,
    chr_col = "Chromosome",
    start_col = "Region_Start",
    end_col = "Region_End",
    keep_intergenic = TRUE,
    verbose = TRUE) {

  if (verbose) cat("Parsing RAP-DB GFF3...\n")
  genes <- hbr_read_gff3_genes(gff3 = gff3, verbose = verbose)
  if (is.null(genes) || !nrow(genes)) return(NULL)

  annotated <- hbr_annotate_regions_with_genes(
    regions = regions, genes = genes,
    chr_col = chr_col, start_col = start_col, end_col = end_col,
    keep_intergenic = keep_intergenic, verbose = verbose)

  if (!is.null(annotated) && nrow(annotated) > 0)
    annotated[, `:=`(Annotation_Source = "RAPDB",
                      Annotation_Build = "IRGSP-1.0")]
  annotated
}
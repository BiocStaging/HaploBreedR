# ============================================================================
# R/breeding-decision.r
# End-to-end breeding decision support: multi-trait parent ranking, cross
# recommendations, and a written report. Includes the honest CV helper,
# report builder, Excel workbook writer, and the shared RDS writer.
#
# All internal helpers use the .hbr_ prefix.
# ============================================================================


# ----------------------------------------------------------------------------
# Global helpers
# ----------------------------------------------------------------------------

#' @noRd
`%||%` <- function(a, b)
  if (is.null(a) || (length(a) == 1L && is.na(a))) b else a

#' Assert a condition
#'
#' Stops with `msg` if `cond` is not `TRUE`.
#'
#' @param cond Condition to test.
#' @param msg Error message.
#'
#' @return `TRUE`, invisibly.
#'
#' @noRd
.hbr_assert <- function(cond, msg) {
  if (!isTRUE(cond)) stop(msg, call. = FALSE)
  invisible(TRUE)
}

#' Log elapsed time
#'
#' Prints elapsed time since `t0` in seconds, minutes, or hours.
#'
#' @param t0 A `POSIXct` start time.
#' @param label Label to print.
#' @param verbose Print the message.
#'
#' @return The elapsed time in seconds, invisibly.
#'
#' @noRd
.hbr_log <- function(t0, label, verbose = TRUE) {
  if (!isTRUE(verbose)) return(invisible(NULL))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (dt < 60)        message(sprintf("  [%-30s] %.1f s", label, dt))
  else if (dt < 3600) message(sprintf("  [%-30s] %.1f min", label, dt / 60))
  else                message(sprintf("  [%-30s] %.2f h", label, dt / 3600))
  invisible(dt)
}

#' Normalize a chromosome name
#'
#' Strips `chr`/`chromosome` prefixes and leading zeros.
#'
#' @param x Character vector of chromosome names.
#'
#' @return Normalized character vector.
#'
#' @noRd
.hbr_norm_chr <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^(chr|Chr|CHR)", "", x)
  x <- sub("^chromosome[_-]?", "", x, ignore.case = TRUE)
  x <- sub("^0+([0-9]+)$", "\\1", x)
  toupper(x)
}

#' Total size of a directory
#'
#' @param dir Directory path.
#'
#' @return Recursive sum of file sizes in bytes.
#'
#' @noRd
.hbr_dir_size <- function(dir) {
  f <- list.files(dir, recursive = TRUE, full.names = TRUE)
  if (!length(f)) 0 else sum(file.size(f), na.rm = TRUE)
}


# ----------------------------------------------------------------------------
# Trait probing
# ----------------------------------------------------------------------------

#' Probe a HaploBreedR trait directory
#'
#' Reads `pipeline_summary` from RDS or CSV and returns trait metadata.
#'
#' @param source Trait directory or RDS path.
#' @param verbose Print progress.
#'
#' @return A list with `source`, `rds_path`, `trait_col`, `direction`,
#'   `mode`, `genes_tested`, `superior_count`, `n_perm`, `ld_decay_bp`.
#'
#' @noRd
.hbr_probe_trait <- function(source, verbose = TRUE) {
  rds_path <- NULL
  if (dir.exists(source)) {
    cand <- file.path(source, "pipeline_result_for_breeding.rds")
    if (file.exists(cand)) rds_path <- cand
  } else if (file.exists(source) &&
             grepl("\\.rds$", source, ignore.case = TRUE)) {
    rds_path <- source
    source   <- dirname(source)
  }
  .hbr_assert(dir.exists(source),
              sprintf("Trait source missing: %s", source))

  info <- list(
    source         = normalizePath(source, winslash = "/",
                                    mustWork = FALSE),
    rds_path       = rds_path,
    trait_col      = NA_character_,
    direction      = NA_character_,
    mode           = NA_character_,
    genes_tested   = NA_integer_,
    superior_count = NA_integer_,
    n_perm         = NA_integer_,
    ld_decay_bp    = NA_real_)

  if (!is.null(rds_path)) {
    x  <- readRDS(rds_path)
    ps <- as.list(x$pipeline_summary[1, ])
    info$trait_col      <- as.character(ps$Trait     %||% NA_character_)
    info$direction      <- as.character(ps$Direction %||% NA_character_)
    info$mode           <- as.character(ps$Mode      %||% NA_character_)
    info$genes_tested   <- suppressWarnings(as.integer(ps$Genes_Tested))
    info$superior_count <- suppressWarnings(
      as.integer(ps$Superior_Haplotypes))
    info$ld_decay_bp    <- suppressWarnings(as.numeric(ps$LD_Decay_bp))
    sr <- x$superior_results
    if (!is.null(sr) && length(sr)) {
      nv <- vapply(sr, function(g) {
        v <- g$test_res$permutation_n_valid
        if (is.null(v) || !length(v)) NA_integer_ else as.integer(v[1])
      }, integer(1))
      nv <- nv[!is.na(nv) & nv > 0L]
      if (length(nv)) info$n_perm <- as.integer(max(nv))
    }
    rm(x); invisible(gc(verbose = FALSE))
  } else {
    ps_csv <- file.path(source, "pipeline_summary.csv")
    if (file.exists(ps_csv)) {
      ps <- data.table::fread(ps_csv)
      info$trait_col      <- as.character(ps$Trait[1])
      info$direction      <- as.character(ps$Direction[1])
      info$mode           <- as.character(ps$Mode[1])
      info$genes_tested   <- suppressWarnings(
        as.integer(ps$Genes_Tested[1]))
      info$superior_count <- suppressWarnings(
        as.integer(ps$Superior_Haplotypes[1]))
      info$ld_decay_bp    <- suppressWarnings(
        as.numeric(ps$LD_Decay_Bp[1]))
    }
  }

  if (verbose)
    message(sprintf("  %s: col=%s dir=%s mode=%s n_perm=%s tested=%s sup=%s",
                    basename(source), info$trait_col, info$direction,
                    info$mode, info$n_perm,
                    info$genes_tested, info$superior_count))
  info
}


# ----------------------------------------------------------------------------
# Reading pipeline outputs
# ----------------------------------------------------------------------------

#' Read a pipeline summary table
#'
#' Loads `haplotype_summary_functional.csv` (or
#' `publication_table_functional.csv`) from a trait directory.
#'
#' @param source Trait directory.
#'
#' @return A `data.table`, or `NULL`.
#'
#' @noRd
.hbr_summary_table <- function(source) {
  hsf <- file.path(source, "haplotype_summary_functional.csv")
  pub <- file.path(source, "publication_table_functional.csv")
  src <- if (file.exists(hsf)) hsf else pub
  if (!file.exists(src)) return(NULL)
  dt <- tryCatch(data.table::fread(src, fill = TRUE),
                 error = function(e) NULL)
  if (is.null(dt) || nrow(dt) < 10L) {
    dt <- tryCatch(read.table(src, header = TRUE, sep = ",", quote = "\"",
                              fill = TRUE, comment.char = "",
                              stringsAsFactors = FALSE),
                   error = function(e) NULL)
    if (!is.null(dt)) dt <- data.table::as.data.table(dt)
  }
  dt
}

#' Superior-haplotype table
#'
#' Reads per-gene superior info from the RDS or the summary CSV.
#'
#' @param info Output of [`.hbr_probe_trait()`].
#'
#' @return A `data.table`.
#'
#' @noRd
.hbr_superior_table <- function(info) {
  if (!is.null(info$rds_path) && file.exists(info$rds_path)) {
    x  <- readRDS(info$rds_path)
    sr <- x$superior_results
    rm(x); invisible(gc(verbose = FALSE))
    if (!is.null(sr) && length(sr)) {
      return(data.table::rbindlist(lapply(names(sr), function(g) {
        s <- sr[[g]]$superior
        if (is.null(s)) return(NULL)
        data.table::data.table(
          gene        = g,
          haplotype   = as.character(s$haplotype),
          nt_sequence = as.character(s$nt_sequence),
          aa_sequence = as.character(s$aa_sequence),
          cohens_d    = as.numeric(s$cohens_d),
          is_superior = as.logical(s$is_superior))
      }), fill = TRUE))
    }
  }
  sm <- .hbr_summary_table(info$source)
  if (is.null(sm)) return(data.table::data.table())
  gcol <- intersect(c("MSU7_Gene_ID","GFF3_Gene_ID","MSU_Gene_ID","Gene"),
                    names(sm))[1]
  if (is.na(gcol) || !"Is_Superior" %in% names(sm))
    return(data.table::data.table())
  keep <- sm$Is_Superior %in% TRUE
  if (!any(keep)) return(data.table::data.table())
  data.table::data.table(
    gene        = as.character(sm[[gcol]][keep]),
    haplotype   = as.character(sm$Superior_Haplotype[keep]),
    nt_sequence = if ("Superior_NT_Seq" %in% names(sm))
                    as.character(sm$Superior_NT_Seq[keep])
                  else NA_character_,
    aa_sequence = if ("Superior_AA_Seq" %in% names(sm))
                    as.character(sm$Superior_AA_Seq[keep])
                  else NA_character_,
    cohens_d    = if ("Cohens_d" %in% names(sm))
                    as.numeric(sm$Cohens_d[keep]) else NA_real_,
    is_superior = TRUE)
}

#' Gene to LD-region map
#'
#' @param info Output of [`.hbr_probe_trait()`].
#'
#' @return A `data.table` with `gene` and `region`.
#'
#' @noRd
.hbr_ld_region_map <- function(info) {
  sm <- .hbr_summary_table(info$source)
  if (is.null(sm) || !"LD_Region_ID" %in% names(sm))
    return(data.table::data.table(gene = character(0),
                                   region = character(0)))
  gcol <- intersect(c("MSU7_Gene_ID","GFF3_Gene_ID","MSU_Gene_ID","Gene"),
                    names(sm))[1]
  if (is.na(gcol))
    return(data.table::data.table(gene = character(0),
                                   region = character(0)))
  dt <- data.table::data.table(
    gene   = as.character(sm[[gcol]]),
    region = as.character(sm$LD_Region_ID))
  dt <- dt[!is.na(gene) & nzchar(gene) &
            !is.na(region) & nzchar(region)]
  unique(dt, by = "gene")
}

#' Extract and cache the reference
#'
#' Pulls `reference` out of a pipeline RDS and stores it once.
#'
#' @param info Output of [`.hbr_probe_trait()`].
#' @param cache_path Cache file path.
#' @param verbose Print progress.
#'
#' @return The reference object.
#'
#' @noRd
.hbr_extract_reference <- function(info, cache_path, verbose = TRUE) {
  if (file.exists(cache_path)) {
    ref <- tryCatch(readRDS(cache_path), error = function(e) NULL)
    if (!is.null(ref)) return(ref)
    unlink(cache_path)
  }
  if (is.null(info$rds_path) || !file.exists(info$rds_path))
    stop("No RDS for reference extraction. ",
         "At least one trait must supply a ",
         "pipeline_result_for_breeding.rds.", call. = FALSE)
  if (verbose) message("Extracting reference (one-time, ~30 s)")
  x <- readRDS(info$rds_path); ref <- x$reference
  rm(x); invisible(gc(verbose = FALSE))
  .hbr_assert(!is.null(ref), "RDS has no reference.")
  tmp <- paste0(cache_path, ".partial_", Sys.getpid())
  saveRDS(ref, tmp, compress = FALSE)
  file.rename(tmp, cache_path)
  ref
}

#' Read per-individual phenotypes
#'
#' @param info Output of [`.hbr_probe_trait()`].
#' @param trait_col Trait column name.
#' @param verbose Print progress.
#'
#' @return A `data.table` with `ID` and the trait column.
#'
#' @noRd
.hbr_phenotypes <- function(info, trait_col, verbose = TRUE) {
  ph <- NULL
  if (!is.null(info$rds_path) && file.exists(info$rds_path)) {
    x  <- readRDS(info$rds_path)
    sr <- x$superior_results
    if (!is.null(sr) && length(sr)) {
      rows <- lapply(sr, function(g) {
        d <- g$test_res$data
        if (is.null(d) || !trait_col %in% names(d) || !"ID" %in% names(d))
          return(NULL)
        data.table::data.table(
          ID    = as.character(d$ID),
          value = suppressWarnings(as.numeric(d[[trait_col]])))
      })
      rows <- rows[!vapply(rows, is.null, logical(1))]
      if (length(rows)) {
        ph <- data.table::rbindlist(rows)
        ph <- ph[!is.na(value)]
        ph <- unique(ph, by = "ID")
      }
    }
    rm(x); invisible(gc(verbose = FALSE))
  }
  if (is.null(ph)) {
    files <- c(
      list.files(file.path(info$source, "superior_genes"),
                 pattern = "_individuals\\.csv$", full.names = TRUE),
      list.files(info$source, pattern = "_individuals\\.csv$",
                 full.names = TRUE))
    files <- unique(files)
    for (f in files) {
      d <- tryCatch(data.table::fread(f), error = function(e) NULL)
      if (is.null(d) || !"ID" %in% names(d) ||
          !trait_col %in% names(d)) next
      p <- data.table::data.table(
        ID    = as.character(d$ID),
        value = suppressWarnings(as.numeric(d[[trait_col]])))
      p <- p[!is.na(value)]
      ph <- if (is.null(ph)) p else
        unique(data.table::rbindlist(list(ph, p)), by = "ID")
    }
  }
  .hbr_assert(!is.null(ph),
              sprintf("No phenotype values for column '%s'.", trait_col))
  data.table::setnames(ph, "value", trait_col)
  if (verbose)
    message(sprintf("  %s: %d individuals", trait_col, nrow(ph)))
  ph
}


# ----------------------------------------------------------------------------
# Genotype + GRM
# ----------------------------------------------------------------------------

#' Read PLINK files into the breeding-decision data shape
#'
#' Thin adapter around [`.hbr_read_genotype_plink()`] so the
#' breeding-decision module benefits from the chunked reader
#' (memory-bounded) and the automatic `.raw` generation when it is
#' missing. Returns the trio `G`, `ids`, `snps` expected by the caller.
#'
#' @param plink_prefix PLINK file prefix.
#' @param plink_path Path to the PLINK executable.
#' @param verbose Print progress.
#'
#' @return A list with `G` (numeric matrix), `ids` (character), `snps`
#'   (data.table).
#'
#' @noRd
.hbr_read_plink <- function(plink_prefix, plink_path = "plink",
                            verbose = TRUE) {
  if (verbose) message("Reading PLINK .raw (chunked)")
  res <- .hbr_read_genotype_plink(plink_prefix, plink_path)
  list(G    = res$genotypes,
       ids  = res$individuals,
       snps = res$variants)
}

#' Build a genome-wide GRM
#'
#' VanRaden-style inbred-line genomic relationship matrix.
#'
#' @param G Numeric dosage matrix.
#' @param maf_min MAF filter.
#' @param verbose Print progress.
#'
#' @return A list with `G`, `A`, `d`, `keep`, `p`, `n_used`, `n_total`.
#'
#' @noRd
.hbr_build_grm <- function(G, maf_min = 0.05, verbose = TRUE) {
  if (verbose) message("Building genome-wide GRM")
  p <- colMeans(G, na.rm = TRUE) / 2
  maf <- pmin(p, 1 - p)
  keep <- !is.na(maf) & maf >= maf_min
  Z <- sweep(G[, keep, drop = FALSE], 2L, 2 * p[keep], "-")
  Z[is.na(Z)] <- 0
  A <- tcrossprod(Z)
  d <- 2 * (2 * sum(p[keep] * (1 - p[keep])))
  Gr <- A / d
  if (verbose)
    message(sprintf("  %d / %d SNPs pass MAF >= %.3f; diag mean %.3f",
                    sum(keep), ncol(G), maf_min, mean(diag(Gr))))
  list(G = Gr, A = A, d = d, keep = keep, p = p[keep],
       n_used = sum(keep), n_total = ncol(G))
}

#' Build a feature GRM
#'
#' Same as [`.hbr_build_grm()`] but restricted to a SNP subset.
#'
#' @param G Numeric dosage matrix.
#' @param feature_snps SNP IDs.
#' @param maf_min MAF filter.
#'
#' @return A list with `A_S`, `d_S`, `keep_idx`.
#'
#' @noRd
.hbr_feature_grm <- function(G, feature_snps, maf_min = 0.05) {
  n <- nrow(G)
  idx <- match(feature_snps, colnames(G))
  idx <- idx[!is.na(idx)]
  if (!length(idx))
    return(list(A_S = matrix(0, n, n), d_S = 0,
                keep_idx = integer(0)))
  Z <- G[, idx, drop = FALSE]
  p <- colMeans(Z, na.rm = TRUE) / 2
  maf <- pmin(p, 1 - p)
  keep <- !is.na(maf) & maf >= maf_min
  if (!any(keep))
    return(list(A_S = matrix(0, n, n), d_S = 0,
                keep_idx = integer(0)))
  Z <- sweep(Z[, keep, drop = FALSE], 2L, 2 * p[keep], "-")
  Z[is.na(Z)] <- 0
  A_S <- tcrossprod(Z)
  d_S <- 2 * (2 * sum(p[keep] * (1 - p[keep])))
  list(A_S = A_S, d_S = d_S, keep_idx = idx[keep])
}


# ----------------------------------------------------------------------------
# Reconstruction of tested-gene haplotype codes
# ----------------------------------------------------------------------------

#' Reconstruct haplotype codes for tested genes
#'
#' Fast path: reads `functional_snps_all_tested_genes.csv`. Fallback:
#' per-gene variant scan via [`.hbr_variant_scan()`].
#'
#' @param info Trait info.
#' @param trait_name Trait label.
#' @param G Numeric dosage matrix.
#' @param ref Reference object.
#' @param params Gate parameters.
#' @param verbose Print progress.
#'
#' @return A list with `nt_seqs`, `hap_map`, `func_map`, `tested`,
#'   `n_genes`.
#'
#' @noRd
.hbr_reconstruct <- function(info, trait_name, G, ref, params,
                              verbose = TRUE) {
  sm <- .hbr_summary_table(info$source)
  if (is.null(sm)) return(NULL)

  gcols <- intersect(c("MSU7_Gene_ID","GFF3_Gene_ID","MSU_Gene_ID","Gene"),
                     names(sm))
  if (!length(gcols)) return(NULL)
  counts <- vapply(gcols, function(cn) {
    v <- sm[[cn]]
    length(unique(v[!is.na(v) & nzchar(as.character(v))]))
  }, integer(1))
  gcol <- gcols[which.max(counts)]
  tested <- unique(as.character(sm[[gcol]]))
  tested <- tested[!is.na(tested) & nzchar(tested)]
  if (verbose)
    message(sprintf("  %s: %d tested genes (col %s)",
                    trait_name, length(tested), gcol))

  src_csv <- file.path(info$source,
                        "functional_snps_all_tested_genes.csv")
  if (file.exists(src_csv)) {
    if (verbose) message("  Loading pipeline functional SNP CSV")
    func_map <- data.table::as.data.table(
      data.table::fread(src_csv, fill = TRUE))
    if (!"gene" %in% names(func_map)) {
      g1 <- intersect(c("Gene_ID","MSU7_Gene_ID","GFF3_Gene_ID","Gene"),
                      names(func_map))[1]
      if (!is.na(g1)) data.table::setnames(func_map, g1, "gene")
    }
    if ("gene" %in% names(func_map))
      func_map <- func_map[gene %in% tested]
  } else {
    if (verbose)
      message("  Variant identification per gene (fallback)")
    if (is.null(ref))
      stop("No functional_snps_all_tested_genes.csv in ", info$source,
           " and no reference available for fallback scan.", call. = FALSE)
    func_map <- .hbr_variant_scan(tested, info, ref, verbose = verbose)
    if (is.null(func_map) || !nrow(func_map)) {
      if (verbose) message("  No functional SNPs identified.")
      return(NULL)
    }
  }

  if (!"Genomic_Order" %in% names(func_map))
    func_map[, Genomic_Order := ifelse(Strand %in% "-", -POS, POS)]
  func_map[is.na(Genomic_Order), Genomic_Order := POS]

  idx <- match(as.character(func_map$SNP), colnames(G))
  if (anyNA(idx)) {
    if (verbose)
      message(sprintf("  %d SNPs not in G; dropping", sum(is.na(idx))))
    ok <- !is.na(idx); func_map <- func_map[ok]; idx <- idx[ok]
  }
  if (!length(idx)) return(NULL)

  Gm <- G[, idx, drop = FALSE]
  colnames(Gm) <- as.character(func_map$SNP)
  rownames(Gm) <- rownames(G)

  func_map <- func_map[order(gene, Genomic_Order)]
  sub_list <- split(func_map, by = "gene", keep.by = FALSE,
                    sorted = FALSE)
  genes_u  <- names(sub_list)
  nt_seqs  <- vector("list", length(genes_u)); names(nt_seqs) <- genes_u
  hap_map  <- vector("list", length(genes_u))

  for (gi in seq_along(genes_u)) {
    g   <- genes_u[gi]
    sub <- sub_list[[gi]]
    sub <- sub[order(Genomic_Order)]

    gd <- tryCatch(
      suppressMessages(hbr_extract_genotypes(
        sub, Gm,
        max_miss_snp    = params$max_miss_snp,
        maf_min         = params$maf_min,
        max_het_ind     = params$max_het_ind,
        max_miss_ind    = params$max_miss_ind,
        min_individuals = params$min_individuals_qc,
        gene_label      = g)),
      error = function(e) NULL)
    if (is.null(gd)) next

    fs <- gd$functional_snps
    if (!"AA_REF" %in% names(fs))
      fs[, AA_REF := ifelse(
        !is.na(AA_Change) & grepl("^[A-Z*][0-9]+[A-Z*]$", AA_Change),
        substr(AA_Change, 1, 1), NA_character_)]
    if (!"AA_ALT" %in% names(fs))
      fs[, AA_ALT := ifelse(
        !is.na(AA_Change) & grepl("^[A-Z*][0-9]+[A-Z*]$", AA_Change),
        substr(AA_Change, nchar(AA_Change), nchar(AA_Change)),
        NA_character_)]
    gd$functional_snps <- fs

    hm <- tryCatch(
      suppressMessages(hbr_build_haplotype_matrix(gd)),
      error = function(e) NULL)
    if (is.null(hm)) next

    hd <- tryCatch(
      suppressMessages(hbr_define_haplotypes(
        hm, min_individuals = params$min_ind_per_hap,
        min_frequency       = params$min_hap_freq)),
      error = function(e) NULL)
    if (is.null(hd)) next

    v <- rep(NA_character_, nrow(G)); names(v) <- rownames(G)
    jj <- match(as.character(hd$haplotypes$ID), rownames(G))
    ok <- !is.na(jj)
    v[jj[ok]] <- as.character(hd$haplotypes$NT_Seq[ok])
    nt_seqs[[g]] <- v

    hm2 <- data.table::as.data.table(hd$haplotype_map)
    hm2[, gene := g]
    hap_map[[gi]] <- hm2
  }

  nt_seqs <- nt_seqs[!vapply(nt_seqs, is.null, logical(1))]
  hap_map <- data.table::rbindlist(hap_map, fill = TRUE)

  list(nt_seqs = nt_seqs, hap_map = hap_map,
       func_map = func_map, tested = tested,
       n_genes = length(nt_seqs))
}

#' Per-gene variant scan (fallback for .hbr_reconstruct)
#'
#' Reads `snps.rds` next to the trait directory and calls
#' [hbr_get_functional_variants()] gene by gene. QC thresholds are
#' applied downstream in [`.hbr_reconstruct()`] via
#' [hbr_extract_genotypes()].
#'
#' @param tested Character vector of gene IDs.
#' @param info Trait info.
#' @param ref Reference object.
#' @param verbose Print progress.
#'
#' @return A `data.table` of functional SNPs across genes.
#'
#' @noRd
.hbr_variant_scan <- function(tested, info, ref, verbose = TRUE) {

  snps_rds <- file.path(info$source, "snps.rds")
  if (!file.exists(snps_rds)) {
    stop("Fallback variant scan requires snps.rds in ", info$source,
         ", or a functional_snps_all_tested_genes.csv in the same ",
         "directory. Neither found.", call. = FALSE)
  }

  snps <- readRDS(snps_rds)
  ref_alt <- data.table::data.table(
    SNP = as.character(snps$SNP),
    CHR = as.character(snps$CHR),
    POS = as.integer(snps$POS),
    REF = as.character(snps$REF),
    ALT = as.character(snps$ALT))
  ref_alt <- ref_alt[!is.na(POS) & !is.na(REF) & !is.na(ALT)]
  ref_alt[, CHR_N := .hbr_norm_chr(CHR)]
  data.table::setkey(ref_alt, CHR_N, POS)

  cds_dt    <- data.table::as.data.table(ref$cds)
  utr_dt    <- data.table::as.data.table(ref$utr)
  cds_seqs  <- ref$cds_seqs
  pep_seqs  <- ref$pep_seqs
  all_genes <- data.table::as.data.table(ref$genes)
  all_genes[, CHR_N := .hbr_norm_chr(Chr)]

  gene_range <- all_genes[gene_id %in% tested &
                           !is.na(Start) & !is.na(End),
                          .(gene_id, CHR_N,
                            start = as.numeric(Start),
                            end   = as.numeric(End))]
  data.table::setkey(gene_range, gene_id)
  ref_alt_by_chr <- split(ref_alt, by = "CHR_N", keep.by = FALSE)

  t0 <- Sys.time()
  parts <- vector("list", length(tested))
  for (i in seq_along(tested)) {
    g  <- tested[i]
    gr <- gene_range[g]
    if (!nrow(gr)) next
    block <- ref_alt_by_chr[[gr$CHR_N[1]]]
    if (is.null(block) || !nrow(block)) next
    ref_alt_gene <- block[POS >= (gr$start[1] - 5000) &
                          POS <= (gr$end[1]   + 5000)]
    if (!nrow(ref_alt_gene)) next

    res <- tryCatch({
      obj <- NULL
      invisible(utils::capture.output(
        obj <- suppressMessages(suppressWarnings(
          hbr_get_functional_variants(
            rapdb_id            = g,
            cds_dt              = cds_dt,
            utr_dt              = utr_dt,
            cds_seqs            = cds_seqs,
            pep_seqs            = pep_seqs,
            all_genes           = all_genes,
            ref_alt             = ref_alt_gene,
            use_all_transcripts = TRUE,
            regions_to_include  = c("CDS_nonsynonymous",
                                     "5UTR", "3UTR")))),
        type = "output"))
      obj
    }, error = function(e) NULL)

    if (is.null(res) || is.null(res$functional_snps) ||
        !nrow(res$functional_snps)) next

    fs <- data.table::as.data.table(res$functional_snps)
    fs[, gene := g]
    keep <- intersect(c("gene","SNP","POS","REF","ALT","Chr","Strand",
                         "Region","AA_Change","Variant_Label",
                         "Codon_Num","Is_Nonsynonymous","Genomic_Order",
                         "AA_REF","AA_ALT"),
                      names(fs))
    parts[[i]] <- fs[, ..keep]

    if (verbose && i %% 100L == 0L)
      message(sprintf("    variant scan: %d / %d (%.1f s)",
                      i, length(tested),
                      as.numeric(difftime(Sys.time(), t0,
                                          units = "secs"))))
  }
  data.table::rbindlist(parts[!vapply(parts, is.null, logical(1))],
                         fill = TRUE)
}


# ----------------------------------------------------------------------------
# REML and mixed model fitting
# ----------------------------------------------------------------------------

#' Single-kernel REML
#'
#' One-variance-component REML via profile likelihood.
#'
#' @param y Numeric response.
#' @param X Fixed-effect design matrix.
#' @param K Kinship matrix.
#'
#' @return A list with fitted components.
#'
#' @noRd
.hbr_reml <- function(y, X, K) {
  n <- length(y); p <- ncol(X); I <- diag(n)
  neg_ll <- function(log_delta) {
    delta <- exp(log_delta)
    W <- K + delta * I
    Wi <- tryCatch(solve(W), error = function(e) NULL)
    if (is.null(Wi)) return(1e10)
    XVX <- crossprod(X, Wi %*% X)
    XVXi <- tryCatch(solve(XVX), error = function(e) NULL)
    if (is.null(XVXi)) return(1e10)
    beta <- XVXi %*% crossprod(X, Wi %*% y)
    r <- y - X %*% beta
    rWr <- as.numeric(crossprod(r, Wi %*% r))
    if (!is.finite(rWr) || rWr <= 0) return(1e10)
    ldW <- as.numeric(determinant(W,   logarithm = TRUE)$modulus)
    ldX <- as.numeric(determinant(XVX, logarithm = TRUE)$modulus)
    -(-0.5 * (ldW + ldX + (n - p) * log(rWr / (n - p))))
  }
  opt <- optimize(neg_ll, c(-15, 15))
  delta_hat <- exp(opt$minimum)
  W <- K + delta_hat * I; Wi <- solve(W)
  XVX <- crossprod(X, Wi %*% X)
  beta <- solve(XVX, crossprod(X, Wi %*% y))
  r <- y - X %*% beta
  rWr <- as.numeric(crossprod(r, Wi %*% r))
  s2g <- rWr / (n - p); s2e <- delta_hat * s2g
  g_hat <- as.vector(K %*% (Wi %*% r))
  list(sigma2_g = s2g, sigma2_e = s2e, beta = as.vector(beta),
       g_hat = g_hat, logLik = -opt$objective, W_inv = Wi,
       y = y, X = X, K = K, n = n, p = p)
}

#' Two-kernel REML fit
#'
#' @param y Numeric response.
#' @param X Fixed-effect design matrix.
#' @param G_R Random-effect kinship matrix.
#' @param K_H Optional second (feature) kinship matrix.
#' @param w_grid Weight grid for the mixture.
#'
#' @return A fitted-model list.
#'
#' @noRd
.hbr_fit <- function(y, X = NULL, G_R, K_H = NULL,
                     w_grid = seq(0, 1, by = 0.1)) {
  n <- length(y)
  if (is.null(X)) X <- matrix(1, n, 1)
  keep <- !is.na(y)
  if (!all(keep)) {
    y <- y[keep]; X <- X[keep, , drop = FALSE]
    G_R <- G_R[keep, keep, drop = FALSE]
    if (!is.null(K_H)) K_H <- K_H[keep, keep, drop = FALSE]
    n <- length(y)
  }
  scG <- mean(diag(G_R))
  scH <- if (!is.null(K_H)) mean(diag(K_H)) else NA_real_
  G_R <- G_R / scG
  if (!is.null(K_H)) K_H <- K_H / scH

  if (is.null(K_H)) {
    fit <- .hbr_reml(y, X, G_R)
    fit$w_hat     <- NA_real_
    fit$sigma2_R  <- fit$sigma2_g
    fit$sigma2_H  <- 0
    fit$K_w       <- G_R
    fit$K_H       <- NULL
    fit$G_R       <- G_R
    fit$scale_G   <- scG
    fit$scale_H   <- scH
    return(fit)
  }

  ll <- vapply(w_grid, function(w)
    .hbr_reml(y, X, (1 - w) * G_R + w * K_H)$logLik, numeric(1))
  bw <- w_grid[which.max(ll)]
  if (length(w_grid) >= 2L) {
    opt <- tryCatch(
      optimize(function(w)
        .hbr_reml(y, X, (1 - w) * G_R + w * K_H)$logLik,
        c(max(0, bw - 0.1), min(1, bw + 0.1)), maximum = TRUE),
      error = function(e) NULL)
    if (!is.null(opt) && opt$objective > max(ll)) bw <- opt$maximum
  }
  Kw <- (1 - bw) * G_R + bw * K_H
  fit <- .hbr_reml(y, X, Kw)
  fit$w_hat     <- bw
  fit$sigma2_R  <- (1 - bw) * fit$sigma2_g
  fit$sigma2_H  <- bw * fit$sigma2_g
  fit$K_w       <- Kw
  fit$K_H       <- K_H
  fit$G_R       <- G_R
  fit$scale_G   <- scG
  fit$scale_H   <- scH
  fit$profile_w <- data.frame(w = w_grid, logLik = ll)
  fit
}

#' Reliability and PEV
#'
#' @param fit Fitted model from [`.hbr_fit()`].
#'
#' @return A list with `reliability` and `PEV`.
#'
#' @noRd
.hbr_reliability <- function(fit) {
  K <- fit$K_w
  KWinvK <- K %*% (fit$W_inv %*% K)
  PEV <- fit$sigma2_g * (diag(K) - diag(KWinvK))
  r <- 1 - PEV / (fit$sigma2_g * diag(K))
  list(reliability = r, PEV = PEV)
}


# ----------------------------------------------------------------------------
# Per-trait full-panel fit
# ----------------------------------------------------------------------------

#' Per-trait full-panel fit
#'
#' @param trait Trait label.
#' @param y Numeric response.
#' @param ids Individual IDs.
#' @param grm Output of [`.hbr_build_grm()`].
#' @param G Numeric dosage matrix.
#' @param recon Output of [`.hbr_reconstruct()`].
#' @param sup_tbl Superior-haplotype table.
#' @param maf_min MAF filter for feature GRM.
#'
#' @return A list with `fit_M0`, `fit_M2`, GEBVs, reliabilities, etc.
#'
#' @noRd
.hbr_fit_panel <- function(trait, y, ids, grm, G, recon, sup_tbl,
                            maf_min = 0.05) {
  names(y) <- ids
  X <- matrix(1, length(y), 1)

  fit_m0 <- .hbr_fit(y, X, G_R = grm$G, K_H = NULL)

  fit_m2 <- NULL
  n_features <- 0L
  if (!is.null(sup_tbl) && nrow(sup_tbl) &&
      !is.null(recon) && !is.null(recon$func_map) &&
      nrow(recon$func_map)) {
    sup_genes <- unique(as.character(sup_tbl$gene))
    feat_snps <- unique(as.character(
      recon$func_map$SNP[recon$func_map$gene %in% sup_genes]))
    feat_snps <- feat_snps[!is.na(feat_snps)]
    if (length(feat_snps) >= 5L) {
      sp <- .hbr_feature_grm(G, feat_snps, maf_min = maf_min)
      if (sp$d_S > 0) {
        A_R <- grm$A - sp$A_S
        d_R <- grm$d - sp$d_S
        G_R <- A_R / d_R
        K_H <- sp$A_S / sp$d_S
        fit_m2 <- .hbr_fit(y, X, G_R = G_R, K_H = K_H)
        n_features <- length(sp$keep_idx)
      }
    }
  }

  rel_m0 <- .hbr_reliability(fit_m0)
  rel_m2 <- if (!is.null(fit_m2)) .hbr_reliability(fit_m2) else NULL

  list(trait = trait, ids = ids, y = y,
       fit_M0 = fit_m0, fit_M2 = fit_m2,
       g_M0   = fit_m0$g_hat,
       g_M2   = if (!is.null(fit_m2)) fit_m2$g_hat else fit_m0$g_hat,
       rel_M0 = rel_m0$reliability,
       rel_M2 = if (!is.null(rel_m2)) rel_m2$reliability
                else rel_m0$reliability,
       PEV_M2 = if (!is.null(rel_m2)) rel_m2$PEV else rel_m0$PEV,
       sigma2_g = if (!is.null(fit_m2)) fit_m2$sigma2_g
                  else fit_m0$sigma2_g,
       sigma2_e = if (!is.null(fit_m2)) fit_m2$sigma2_e
                  else fit_m0$sigma2_e,
       w_hat    = if (!is.null(fit_m2)) fit_m2$w_hat else NA_real_,
       n_features = n_features)
}


# ----------------------------------------------------------------------------
# Gate parameters
# ----------------------------------------------------------------------------

#' Superior-haplotype gate parameters
#'
#' Builds the parameter list used by reconstruction and CV.
#'
#' @param info Trait info.
#' @param overrides Named list of overrides.
#'
#' @return A parameter list.
#'
#' @noRd
.hbr_gate_params <- function(info, overrides = list()) {
  p <- list(
    alpha              = 0.05,
    min_cohens_d       = 0.3,
    max_miss_snp       = 0.30,
    maf_min            = 0.00,
    max_het_ind        = 0.10,
    max_miss_ind       = 0.30,
    min_individuals_qc = 5L,
    min_ind_per_hap    = 3L,
    min_hap_freq       = 0,
    n_perm             = info$n_perm %||% 1000L,
    posthoc            = "tukey",
    pc_cols            = character(0),
    trait_threshold    = NULL,
    direction          = info$direction,
    mode               = info$mode)
  for (nm in names(overrides)) {
    v <- overrides[[nm]]
    if (!is.null(v)) p[[nm]] <- v
  }
  p
}


# ----------------------------------------------------------------------------
# Multi-trait index
# ----------------------------------------------------------------------------

#' Multi-trait selection index
#'
#' Weighted sum of scaled, direction-corrected GEBVs.
#'
#' @param fits Named list of per-trait fit objects.
#' @param weights Named numeric weights.
#' @param directions Named direction vector.
#' @param scale_by_sd Scale GEBVs by trait genetic SD.
#'
#' @return A `data.table` with per-individual index, SE, rank,
#'   contributions.
#'
#' @noRd
.hbr_index <- function(fits, weights, directions, scale_by_sd = TRUE) {
  traits <- intersect(names(fits), names(weights))
  .hbr_assert(length(traits) > 0L,
             "No traits in both `fits` and `weights`.")
  fits <- fits[traits]
  common <- Reduce(intersect, lapply(fits, function(f) f$ids))
  .hbr_assert(length(common) > 0L,
             "No common individuals across traits.")

  idx <- rep(0, length(common)); names(idx) <- common
  se2 <- rep(0, length(common)); names(se2) <- common
  contrib <- data.table::data.table(id = common)

  for (tr in traits) {
    f  <- fits[[tr]]
    j  <- match(common, f$ids)
    g  <- f$g_M2[j]
    pev <- f$PEV_M2[j]
    dir <- directions[[tr]] %||% "higher"
    w   <- as.numeric(weights[[tr]] %||% 1)
    sc  <- if (isTRUE(scale_by_sd)) sqrt(f$sigma2_g) else 1
    gg  <- g / sc
    pv  <- pev / sc^2
    if (identical(dir, "lower")) gg <- -gg
    idx <- idx + w * gg
    se2 <- se2 + w^2 * pv
    contrib[[tr]] <- gg
  }
  se <- sqrt(se2)
  ord <- order(-idx)
  rank <- integer(length(idx)); rank[ord] <- seq_along(idx)
  data.table::data.table(id = common, contrib[, -1, with = FALSE],
                          index = idx, SE = se, rank = rank)
}

#' Rank stability
#'
#' Bootstrap probability that each individual is in the top `k`.
#'
#' @param index_est Numeric index.
#' @param index_se Numeric standard error.
#' @param top_k Number of top individuals per draw.
#' @param n_draw Number of bootstrap draws.
#' @param seed Random seed.
#'
#' @return A numeric vector of probabilities.
#'
#' @noRd
.hbr_rank_stability <- function(index_est, index_se, top_k = 15L,
                                n_draw = 5000L, seed = 1L) {
  n <- length(index_est); hits <- numeric(n)
  draw_matrix <- .with_seed(seed, {
    m <- matrix(stats::rnorm(n * n_draw), nrow = n, ncol = n_draw)
    m
  })
  for (b in seq_len(n_draw)) {
    draw <- index_est + index_se * draw_matrix[, b]
    top  <- order(-draw)[seq_len(min(top_k, n))]
    hits[top] <- hits[top] + 1
  }
  hits / n_draw
}


# ----------------------------------------------------------------------------
# Cross scoring (LD-region aware)
# ----------------------------------------------------------------------------

#' LD-region complementarity
#'
#' Proportion of LD regions where exactly one parent carries the superior
#' haplotype.
#'
#' @param g1,g2 Parent IDs.
#' @param nt_seqs_lookup Gene to NT sequences lookup.
#' @param sup_dt Superior table.
#' @param ld_map LD-region map.
#'
#' @return A list with region counts and coverage.
#'
#' @noRd
.hbr_region_complement <- function(g1, g2, nt_seqs_lookup, sup_dt, ld_map) {
  if (is.null(sup_dt) || !nrow(sup_dt))
    return(list(n_regions = 0L, n_either = 0L, n_both = 0L,
                coverage = 0, complement = 0))
  s <- data.table::as.data.table(sup_dt)
  s <- merge(s, ld_map, by = "gene", all.x = TRUE)
  s[is.na(region) | !nzchar(region),
    region := paste0("__gene_", gene)]

  region_hits <- s[, {
    p1_carry <- FALSE; p2_carry <- FALSE
    for (g in unique(gene)) {
      seqs <- nt_seqs_lookup[[g]]
      if (is.null(seqs)) next
      s1 <- seqs[g1]; s2 <- seqs[g2]
      nts <- nt_sequence[gene == g]; nts <- nts[!is.na(nts)]
      for (nt in nts) {
        if (!is.na(s1) && s1 == nt) p1_carry <- TRUE
        if (!is.na(s2) && s2 == nt) p2_carry <- TRUE
      }
    }
    list(p1_carries = p1_carry, p2_carries = p2_carry)
  }, by = region]

  n_regions <- nrow(region_hits)
  n_either  <- sum(region_hits$p1_carries | region_hits$p2_carries)
  n_both    <- sum(region_hits$p1_carries & region_hits$p2_carries)
  list(n_regions = n_regions, n_either = n_either, n_both = n_both,
       coverage   = if (n_regions) n_either / n_regions else 0,
       complement = if (n_regions) (n_either - n_both) / n_regions else 0)
}

#' Rank pairwise crosses
#'
#' Enumerates and scores all pairs of candidate parents.
#'
#' @param parents Character vector of parent IDs.
#' @param idx The multi-trait index table.
#' @param G GRM.
#' @param nt_seqs_lookup Gene to NT sequences lookup.
#' @param sup_sets Superior sets per trait.
#' @param ld_regions LD-region maps per trait.
#' @param relatedness_max Relatedness ceiling.
#'
#' @return A `data.table` of scored crosses.
#'
#' @noRd
.hbr_pair_crosses <- function(parents, idx, G, nt_seqs_lookup,
                              sup_sets, ld_regions,
                              relatedness_max = 0.30) {
  n <- length(parents)
  .hbr_assert(n >= 2L, "Need >= 2 parents.")
  idx_m <- idx[match(parents, idx$id)]
  g_idx <- idx_m$index; names(g_idx) <- parents

  rows <- list()
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      p1 <- parents[i]; p2 <- parents[j]
      rel <- G[p1, p2]
      tr_summ <- lapply(names(sup_sets), function(tr) {
        .hbr_region_complement(p1, p2, nt_seqs_lookup,
                                sup_sets[[tr]], ld_regions[[tr]])
      })
      n_reg    <- sum(vapply(tr_summ, `[[`, integer(1), "n_regions"))
      n_either <- sum(vapply(tr_summ, `[[`, integer(1), "n_either"))
      n_both   <- sum(vapply(tr_summ, `[[`, integer(1), "n_both"))
      rows[[length(rows) + 1L]] <- data.table::data.table(
        p1 = p1, p2 = p2,
        mid_index   = (g_idx[p1] + g_idx[p2]) / 2,
        relatedness = rel,
        coverage    = if (n_reg) n_either / n_reg else 0,
        complement  = if (n_reg) (n_either - n_both) / n_reg else 0,
        n_regions   = n_reg)
    }
  }
  out <- data.table::rbindlist(rows)
  out <- out[relatedness <= relatedness_max]
  out[, score := mid_index + 2 * complement]
  data.table::setorder(out, -score)
  out[, rank := .I]
  out[]
}


# ----------------------------------------------------------------------------
# Trait resolution
# ----------------------------------------------------------------------------

#' Resolve `trait_dirs` input
#'
#' Accepts a named character vector or a list of `list(name=, source=)`.
#'
#' @param trait_dirs Trait directory input.
#' @param verbose Print progress.
#'
#' @return A named list of trait specifications.
#'
#' @noRd
.hbr_resolve_traits <- function(trait_dirs, verbose = TRUE) {
  if (is.character(trait_dirs) && !is.null(names(trait_dirs))) {
    .hbr_assert(all(nzchar(names(trait_dirs))),
               "All trait names must be non-empty.")
    out <- vector("list", length(trait_dirs))
    for (i in seq_along(trait_dirs)) {
      info <- .hbr_probe_trait(trait_dirs[[i]], verbose = verbose)
      out[[i]] <- list(
        name      = names(trait_dirs)[i],
        source    = info$source,
        trait_col = info$trait_col,
        direction = info$direction,
        info      = info)
    }
    names(out) <- names(trait_dirs)
    return(out)
  }
  if (is.list(trait_dirs)) {
    .hbr_assert(all(vapply(trait_dirs, function(x)
      is.list(x) && !is.null(x$name) && !is.null(x$source), logical(1))),
      "Each element must have `name` and `source`.")
    out <- vector("list", length(trait_dirs))
    for (i in seq_along(trait_dirs)) {
      t <- trait_dirs[[i]]
      info <- .hbr_probe_trait(t$source, verbose = verbose)
      out[[i]] <- list(
        name      = t$name,
        source    = info$source,
        trait_col = t$trait_col %||% info$trait_col,
        direction = t$direction %||% info$direction,
        info      = info)
    }
    names(out) <- vapply(out, function(x) x$name, character(1))
    return(out)
  }
  stop("`trait_dirs` must be a named character vector or a list of ",
       "list(name=, source=) elements.", call. = FALSE)
}


# ----------------------------------------------------------------------------
# Per-individual GEBV table writer
# ----------------------------------------------------------------------------

#' Write the per-individual GEBV table
#'
#' Emits `gebv_per_individual.csv` with raw, scaled, and observed columns.
#'
#' @param result The hgsi_result list.
#' @param output_dir Output directory.
#'
#' @return The output `data.table`, invisibly.
#'
#' @noRd
.hbr_write_gebv_table <- function(result, output_dir) {
  ids <- result$fits[[1]]$ids
  trs <- names(result$fits)
  out <- data.table::data.table(id = ids)

  for (tr in trs) {
    f   <- result$fits[[tr]]
    j   <- match(ids, f$ids)
    dir <- result$directions[[tr]]
    sc  <- sqrt(f$sigma2_g)

    out[[paste0(tr, "_GEBV_M0")]]     <- f$g_M0[j]
    out[[paste0(tr, "_GEBV_M2")]]     <- f$g_M2[j]
    out[[paste0(tr, "_GEBV_scaled")]] <- if (identical(dir, "lower"))
      -f$g_M2[j] / sc else f$g_M2[j] / sc
    out[[paste0(tr, "_rel")]]         <- f$rel_M2[j]

    ph <- result$phenos[[tr]]
    if (!is.null(ph)) {
      jj <- match(ids, hbr_clean_id(ph$ID))
      out[[paste0(tr, "_observed")]] <- ph[[tr]][jj]
    }
  }

  idx_m <- result$index[match(ids, id)]
  out[, `:=`(index     = idx_m$index,
              index_SE  = idx_m$SE,
              rank      = idx_m$rank,
              P_top15   = idx_m$P_top15)]
  data.table::fwrite(out, file.path(output_dir, "gebv_per_individual.csv"))
  invisible(out)
}


# ----------------------------------------------------------------------------
# Breeding RDS writer (shared)
# ----------------------------------------------------------------------------

#' Write `pipeline_result_for_breeding.rds` (and optionally the reference)
#'
#' Shared by [haplobreed_analysis()] and [hbr_save_breeding_rds()]. Never
#' loses the run: if the genotype matrix cannot be packed or written (for
#' example for lack of memory) the RDS is written WITHOUT it and a
#' warning says so; the matrix can always be re-read from
#' `plink_prefix`.
#'
#' @param superior_results Superior-gene result list.
#' @param geno_mat Numeric dosage matrix (or packed `hbr_packed_geno`).
#' @param reference Reference object.
#' @param ref_paths Reference path list.
#' @param plink_prefix PLINK file prefix.
#' @param pipeline_summary Pipeline summary data.table.
#' @param output_dir Output directory.
#' @param rds_compress Compression setting.
#' @param rds_reference Reference storage mode.
#' @param verbose Print progress.
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
  geno_ids <- if (!is.null(geno_mat)) rownames(geno_mat) else NULL
  geno_dim <- if (!is.null(geno_mat)) dim(geno_mat) else NULL
  geno_saved <- !is.null(geno_mat)
  packed <- tryCatch(hbr_pack_genotypes(geno_mat), error = function(e) {
    warning("Genotype packing failed (", conditionMessage(e),
            "): the RDS will not contain the genotype matrix.", call. = FALSE)
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


# ----------------------------------------------------------------------------
# CV helper: fold-specific gate
# ----------------------------------------------------------------------------

#' Fold-specific superior-haplotype gate
#'
#' Recomputes the superior call on a training subset.
#'
#' @param train_ids Training IDs.
#' @param y Numeric response named by ID.
#' @param ids All IDs.
#' @param nt_seqs Gene to NT sequences.
#' @param sup_tbl Superior table.
#' @param params Gate parameters.
#' @param direction Trait direction.
#' @param seed Random seed.
#'
#' @return A `data.table` of gate outcomes.
#'
#' @noRd
.hbr_gate_fold <- function(train_ids, y, ids, nt_seqs, sup_tbl,
                           params, direction, seed) {
  genes <- names(nt_seqs)
  # Seed only where randomness is drawn, later in the function
  if (!length(genes)) return(data.table::data.table())

  idx_train <- match(train_ids, ids)
  ok <- !is.na(idx_train)
  idx_train <- idx_train[ok]
  y_tr <- y[idx_train]

  keep_obs <- !is.na(y_tr)
  idx_train <- idx_train[keep_obs]
  y_tr <- y_tr[keep_obs]
  if (length(y_tr) < 20L) return(data.table::data.table())

  out <- vector("list", length(genes))
  for (gi in seq_along(genes)) {
    g <- genes[gi]
    seqs <- nt_seqs[[g]]
    if (is.null(seqs)) next
    codes <- seqs[idx_train]
    if (length(unique(stats::na.omit(codes))) < 2L) next

    dt <- data.table::data.table(ID = ids[idx_train],
                                  hap = as.character(codes),
                                  y = y_tr)
    dt <- dt[!is.na(hap) & !is.na(y)]
    if (nrow(dt) < 20L) next

    grp <- dt[, .(n = .N, mean_y = mean(y),
                  se = stats::sd(y) / sqrt(.N)), by = hap]
    grp <- grp[n >= params$min_ind_per_hap]
    if (nrow(grp) < 2L) next

    dir <- direction
    grp <- grp[order(if (dir == "higher") -mean_y else mean_y)]
    best <- grp[1]
    second <- grp[2]

    g1 <- dt[hap == best$hap]$y
    g2 <- dt[hap == second$hap]$y
    tt <- tryCatch(stats::t.test(g1, g2, var.equal = TRUE),
                    error = function(e) NULL)
    p_val <- if (!is.null(tt)) tt$p.value else NA_real_

    psd <- sqrt(((length(g1) - 1) * stats::var(g1) +
                  (length(g2) - 1) * stats::var(g2)) /
                 (length(g1) + length(g2) - 2))
    d_val <- if (is.finite(psd) && psd > 0)
      abs(mean(g1) - mean(g2)) / psd else NA_real_

    gate_a <- is.finite(p_val) && p_val < params$alpha
    gate_c <- is.finite(d_val) && d_val >= params$min_cohens_d
    is_sup <- gate_a && gate_c

    out[[gi]] <- data.table::data.table(
      gene = g, is_superior = is_sup,
      superior_hap = as.character(best$hap),
      p_value = p_val, cohens_d = d_val,
      seed = seed)
  }
  data.table::rbindlist(out[!vapply(out, is.null, logical(1))], fill = TRUE)
}


# ----------------------------------------------------------------------------
# CV helper: honest nested CV
# ----------------------------------------------------------------------------

#' Honest nested CV for one trait
#'
#' Outer folds with fold-specific gate replay (M0, M2, M3).
#'
#' @param trait Trait label.
#' @param fit Full-data fit.
#' @param trait_spec Trait specification.
#' @param geno Genotype list.
#' @param grm GRM object.
#' @param recon Reconstruction object.
#' @param sup_tbl Superior table.
#' @param params Gate parameters.
#' @param k Number of folds.
#' @param repeats Number of repeats.
#' @param n_perm Permutation iterations inside the fold gate.
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return A list with `predictions` and metadata.
#'
#' @noRd
.hbr_cv_one_honest <- function(trait, fit, trait_spec, geno, grm,
                               recon, sup_tbl, params,
                               k = 3L, repeats = 1L, n_perm = 50L,
                               seed = 1L, verbose = TRUE) {

  ids <- fit$ids
  y   <- fit$y
  G_all <- grm$G
  X <- matrix(1, length(y), 1)

  feat_naive <- character(0)
  if (!is.null(recon) && !is.null(sup_tbl) && nrow(sup_tbl)) {
    sup_genes <- unique(as.character(sup_tbl$gene))
    if (!is.null(recon$func_map) && nrow(recon$func_map)) {
      feat_naive <- unique(as.character(
        recon$func_map$SNP[recon$func_map$gene %in% sup_genes]))
      feat_naive <- feat_naive[!is.na(feat_naive)]
    }
  }

  obs <- which(!is.na(y))
  fold_assignment <- .with_seed(seed, {
    sample(rep(seq_len(k), length.out = n_obs))
  })
  n_obs <- length(obs)
  if (n_obs < 30L) {
    if (verbose) message("  Too few observations for CV; skipping.")
    return(NULL)
  }

  preds_all <- list()
  for (r in seq_len(repeats)) {
    fold <- sample(rep(seq_len(k), length.out = n_obs))
    for (f in seq_len(k)) {
      va <- obs[fold == f]
      tr <- setdiff(obs, va)
      if (length(tr) < 20L || length(va) < 5L) next

      fit_m0 <- .hbr_fit(y[tr], X[tr, , drop = FALSE],
                         G_R = G_all[tr, tr, drop = FALSE])
      G_VT <- G_all[va, tr, drop = FALSE] / fit_m0$scale_G
      r_T  <- fit_m0$y - fit_m0$X %*% fit_m0$beta
      g_m0 <- as.vector(G_VT %*% (fit_m0$W_inv %*% r_T))

      preds_all[[length(preds_all) + 1L]] <- data.table::data.table(
        rep = r, fold = f, method = "M0",
        id = ids[va], y = y[va], pred = g_m0)

      if (length(feat_naive) >= 5L) {
        sp <- .hbr_feature_grm(geno$G, feat_naive, maf_min = 0.05)
        if (sp$d_S > 0) {
          A_R <- grm$A - sp$A_S
          d_R <- grm$d - sp$d_S
          G_R_tt <- A_R[tr, tr, drop = FALSE] / d_R
          K_H_tt <- sp$A_S[tr, tr, drop = FALSE] / sp$d_S
          G_R_vt <- A_R[va, tr, drop = FALSE] / d_R
          K_H_vt <- sp$A_S[va, tr, drop = FALSE] / sp$d_S

          fit_m2 <- .hbr_fit(y[tr], X[tr, , drop = FALSE],
                             G_R = G_R_tt, K_H = K_H_tt)
          r_T2 <- fit_m2$y - fit_m2$X %*% fit_m2$beta
          g_vt <- G_R_vt / fit_m2$scale_G
          k_vt <- K_H_vt / fit_m2$scale_H
          Kw_vt <- (1 - fit_m2$w_hat) * g_vt + fit_m2$w_hat * k_vt
          g_m2 <- as.vector(Kw_vt %*% (fit_m2$W_inv %*% r_T2))

          preds_all[[length(preds_all) + 1L]] <- data.table::data.table(
            rep = r, fold = f, method = "M2",
            id = ids[va], y = y[va], pred = g_m2)
        }
      }

      if (!is.null(recon$nt_seqs)) {
        gate <- .hbr_gate_fold(
          train_ids = ids[tr],
          y = y, ids = ids,
          nt_seqs = recon$nt_seqs,
          sup_tbl = sup_tbl,
          params = params, direction = trait_spec$direction,
          seed = seed + r * 1000L + f)

        feat_fold <- character(0)
        if (nrow(gate) && any(gate$is_superior)) {
          sup_fold <- gate[is_superior == TRUE]$gene
          feat_fold <- unique(as.character(
            recon$func_map$SNP[recon$func_map$gene %in% sup_fold]))
          feat_fold <- feat_fold[!is.na(feat_fold)]
        }

        if (length(feat_fold) >= 5L) {
          sp <- .hbr_feature_grm(geno$G, feat_fold, maf_min = 0.05)
          if (sp$d_S > 0) {
            A_R <- grm$A - sp$A_S
            d_R <- grm$d - sp$d_S
            G_R_tt <- A_R[tr, tr, drop = FALSE] / d_R
            K_H_tt <- sp$A_S[tr, tr, drop = FALSE] / sp$d_S
            G_R_vt <- A_R[va, tr, drop = FALSE] / d_R
            K_H_vt <- sp$A_S[va, tr, drop = FALSE] / sp$d_S

            fit_m3 <- .hbr_fit(y[tr], X[tr, , drop = FALSE],
                               G_R = G_R_tt, K_H = K_H_tt)
            r_T3 <- fit_m3$y - fit_m3$X %*% fit_m3$beta
            g_vt <- G_R_vt / fit_m3$scale_G
            k_vt <- K_H_vt / fit_m3$scale_H
            Kw_vt <- (1 - fit_m3$w_hat) * g_vt + fit_m3$w_hat * k_vt
            g_m3 <- as.vector(Kw_vt %*% (fit_m3$W_inv %*% r_T3))

            preds_all[[length(preds_all) + 1L]] <- data.table::data.table(
              rep = r, fold = f, method = "M3",
              id = ids[va], y = y[va], pred = g_m3)
          }
        }
      }
    }
  }

  preds <- data.table::rbindlist(preds_all, fill = TRUE)
  list(trait = trait, predictions = preds,
       n_folds = k, repeats = repeats)
}

#' Compare CV methods against a baseline
#'
#' Bootstrap CIs on delta-PA and delta-RMSE vs the M0 baseline.
#'
#' @param cv Output of [`.hbr_cv_one_honest()`].
#' @param baseline Baseline method name.
#' @param n_boot Number of bootstrap draws.
#' @param seed Random seed.
#'
#' @return A `data.table` of comparisons.
#'
#' @noRd
.hbr_cv_compare <- function(cv, baseline = "M0", n_boot = 2000L, seed = 1L) {
  preds <- cv$predictions
  if (!nrow(preds)) return(data.table::data.table())

  agg <- preds[, .(y = mean(y, na.rm = TRUE),
                    pred = mean(pred, na.rm = TRUE)),
               by = .(id, method)]
  base <- agg[method == baseline, .(id, y, pred_base = pred)]

  methods <- setdiff(unique(agg$method), baseline)
  
  out <- list()
  for (m in methods) {
    cmp <- merge(base, agg[method == m, .(id, pred_m = pred)], by = "id")
    cmp <- cmp[!is.na(y) & !is.na(pred_base) & !is.na(pred_m)]
    n <- nrow(cmp)
    if (n < 20L) next

    pa_base <- stats::cor(cmp$y, cmp$pred_base)
    pa_m    <- stats::cor(cmp$y, cmp$pred_m)
    rmse_base <- sqrt(mean((cmp$y - cmp$pred_base)^2))
    rmse_m    <- sqrt(mean((cmp$y - cmp$pred_m)^2))

    boot_pa <- numeric(n_boot); boot_rmse <- numeric(n_boot)
    for (b in seq_len(n_boot)) {
      i <- .with_seed(seed + b, sample(n, n, replace = TRUE)); s <- cmp[i]
      boot_pa[b]   <- stats::cor(s$y, s$pred_m) -
                      stats::cor(s$y, s$pred_base)
      boot_rmse[b] <- sqrt(mean((s$y - s$pred_m)^2)) -
                      sqrt(mean((s$y - s$pred_base)^2))
    }

    out[[m]] <- data.table::data.table(
      method = m, baseline = baseline, n = n,
      d_PA = pa_m - pa_base,
      d_PA_lo = unname(stats::quantile(boot_pa,   0.025, na.rm = TRUE)),
      d_PA_hi = unname(stats::quantile(boot_pa,   0.975, na.rm = TRUE)),
      d_RMSE  = rmse_m - rmse_base,
      d_RMSE_lo = unname(stats::quantile(boot_rmse, 0.025, na.rm = TRUE)),
      d_RMSE_hi = unname(stats::quantile(boot_rmse, 0.975, na.rm = TRUE)))
  }
  data.table::rbindlist(out, fill = TRUE)
}


# ----------------------------------------------------------------------------
# Report builder
# ----------------------------------------------------------------------------

#' Build the breeding decision Markdown report
#'
#' Renders the full decision report from an `hgsi_result`.
#'
#' @param result The hgsi_result list.
#' @param verbose Print progress.
#'
#' @return A single character string containing the Markdown report.
#'
#' @noRd
.hbr_report <- function(result, verbose = TRUE) {
  if (verbose) message("Writing decision report")
  ev   <- result$evidence
  par  <- result$parents
  cr   <- result$crosses
  cf   <- result$counterfactual
  traits <- names(result$traits)

  fmt <- function(x, n = 3) formatC(x, format = "f", digits = n)

  lines <- character(0)
  add <- function(...) lines <<- c(lines, ...)

  add("# Haplotype-Informed Breeding Decision Report")
  add("")
  add(sprintf("**Generated:** %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  add(sprintf("**Traits:** %s", paste(traits, collapse = ", ")))
  add(sprintf("**Panel:** %d individuals", nrow(result$G$G)))
  add(sprintf("**SNPs:** %d total, %d passing MAF",
              result$G$n_total, result$G$n_used))
  add(sprintf("**MAF filter:** %.3f", result$maf_min))
  add(sprintf("**Relatedness ceiling:** %.2f",
              result$parameters$relatedness_max))
  add("")

  add("---")
  add("")
  add("## 1. Trait summary")
  add("")
  add("| Trait | Direction | Mode | Tested | Superior | n | h2 | w_hat |")
  add("|---|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(ev))) {
    add(sprintf("| %s | %s | %s | %s | %s | %d | %.3f | %s |",
                ev$trait[i], ev$direction[i], ev$mode[i],
                ev$tested[i], ev$superior[i], ev$n_obs[i],
                ev$h2[i],
                if (is.na(ev$w_hat[i])) "-" else fmt(ev$w_hat[i])))
  }
  add("")

  if (!is.null(result$cv) && any(!vapply(result$cv, is.null, logical(1)))) {
    add("### 1.1 Cross-validation")
    add("")
    add("| Trait | PA | RMSE |")
    add("|---|---|---|")
    for (cv in result$cv) {
      if (is.null(cv)) next
      add(sprintf("| %s | %.3f | %.3f |", cv$trait, cv$PA, cv$RMSE))
    }
    add("")
  }

  add("---")
  add("")
  add("## 2. Multi-trait index")
  add("")
  add(sprintf("Weights: %s",
              paste(sprintf("%s=%g", names(result$weights), result$weights),
                    collapse = ", ")))
  add(sprintf("Directions: %s",
              paste(sprintf("%s=%s", names(result$directions),
                            result$directions), collapse = ", ")))
  add("")

  add("---")
  add("")
  add("## 3. Parent decision sheet")
  add("")
  add("| Rank | ID | Index | P_top | Haplotype coverage | Mean relatedness | Cross count | Action |")
  add("|---|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(par))) {
    add(sprintf("| %d | %s | %.3f | %.3f | %.3f | %.3f | %d | %s |",
                par$rank[i], par$id[i], par$index[i], par$P_top15[i],
                par$haplo_cov[i], par$mean_rel[i], par$n_crosses[i],
                par$action[i]))
  }
  add("")

  n_adv <- sum(par$action == "advance")
  if (n_adv > 0L) {
    add(sprintf("**Advance (%d):** %s",
                n_adv, paste(par[action == "advance"]$id, collapse = ", ")))
  }
  add("")

  add("---")
  add("")
  add("## 4. Recommended crosses")
  add("")
  add("Coverage and complement are computed at the LD-region level.")
  add("")
  add("| Rank | Parent 1 | Parent 2 | Mid-index | Relatedness | Coverage | Complement | Score |")
  add("|---|---|---|---|---|---|---|---|")
  for (i in seq_len(nrow(cr))) {
    add(sprintf("| %d | %s | %s | %.3f | %.3f | %.3f | %.3f | %.3f |",
                cr$rank[i], cr$p1[i], cr$p2[i], cr$mid_index[i],
                cr$relatedness[i], cr$coverage[i],
                cr$complement[i], cr$score[i]))
  }
  add("")

  add("---")
  add("")
  add("## 5. Counterfactual -- GBLUP-only ranking")
  add("")
  add(sprintf("- Parents kept in both rankings: **%d / %d**",
              length(cf$kept), length(cf$top_m3)))
  add(sprintf("- Entered under GBLUP-only: **%d**", length(cf$entered_m0)))
  add(sprintf("- Exited under GBLUP-only: **%d**", length(cf$exited_m0)))
  add("")
  if (length(cf$entered_m0))
    add(sprintf("**Entered:** %s", paste(cf$entered_m0, collapse = ", ")))
  if (length(cf$exited_m0))
    add(sprintf("**Exited:** %s", paste(cf$exited_m0, collapse = ", ")))
  add("")

  add("---")
  add("")
  add("## 6. Provenance")
  add("")
  for (tr in traits) {
    t <- result$traits[[tr]]
    add(sprintf("- **%s**: source `%s`, direction=%s, trait_col=`%s`",
                tr, t$source, t$direction, t$trait_col))
  }
  add("")
  add(sprintf("PLINK prefix: `%s`", result$plink_prefix))
  add(sprintf("Output directory: `%s`", result$output_dir))
  add("")
  add("_Every numeric claim in this report was produced by an executed script._")

  paste(lines, collapse = "\n")
}


# ----------------------------------------------------------------------------
# Excel workbook
# ----------------------------------------------------------------------------

#' Format a numeric or character column for the Excel workbook
#'
#' Fixed-digit formatting with thousands separators.
#'
#' @param x Vector.
#' @param digits Number of decimal places.
#'
#' @return A character vector.
#'
#' @noRd
.hbr_fmt <- function(x, digits = 3) {
  if (is.numeric(x))
    formatC(x, format = "f", digits = digits, big.mark = ",")
  else as.character(x)
}

#' Format a p-value for Excel cells
#'
#' Scientific notation below 1e-4, otherwise fixed digits. Distinct from
#' the figure-annotation formatter [hbr_format_p()].
#'
#' @param p Numeric vector of p-values.
#'
#' @return A character vector.
#'
#' @noRd
.hbr_fmt_p_excel <- function(p) {
  vapply(p, function(x) {
    if (is.na(x)) return("")
    if (x < 1e-4) formatC(x, format = "e", digits = 2)
    else          formatC(x, format = "f", digits = 4)
  }, character(1))
}

#' @noRd
`%+%` <- function(a, b) paste0(a, b)

#' Write the decision workbook
#'
#' @param result The hgsi_result list.
#' @param path Output file.
#'
#' @return The output path, invisibly, or `NULL`.
#'
#' @noRd
.hbr_write_workbook <- function(result,
                                path = file.path(
                                  result$output_dir,
                                  sprintf("breeding_decision_%s.xlsx",
                                          format(Sys.Date(), "%Y%m%d")))) {

  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    warning("Package 'openxlsx' not installed. ",
            "install.packages('openxlsx') to write the workbook.",
            call. = FALSE)
    return(invisible(NULL))
  }

  wb <- openxlsx::createWorkbook()

  hdr_style <- openxlsx::createStyle(
    fontSize = 11, textDecoration = "bold",
    fgFill = "#1F3864", fontColour = "white",
    halign = "center", valign = "center",
    border = "TopBottomLeftRight")
  sub_style <- openxlsx::createStyle(
    fontSize = 11, textDecoration = "bold",
    fgFill = "#D9E1F2", fontColour = "black",
    halign = "left", valign = "center")

  add_table <- function(sheet, title, subtitle, df,
                        col_widths = NULL, freeze_col = 0L) {
    openxlsx::addWorksheet(wb, sheet)
    r <- 1L

    openxlsx::writeData(wb, sheet, title, startRow = r, startCol = 1)
    openxlsx::addStyle(wb, sheet, sub_style, rows = r, cols = 1,
                        stack = FALSE)
    r <- r + 1L

    if (!is.null(subtitle) && nzchar(subtitle)) {
      openxlsx::writeData(wb, sheet, subtitle, startRow = r, startCol = 1)
      r <- r + 1L
    }
    r <- r + 1L

    openxlsx::writeData(wb, sheet, df, startRow = r, startCol = 1,
                         headerStyle = hdr_style, keepNA = TRUE)

    openxlsx::addFilter(wb, sheet,
                         rows = r:(r + nrow(df)),
                         cols = seq_len(ncol(df)))
    if (freeze_col > 0L)
      openxlsx::freezePane(wb, sheet,
                            firstActiveRow = r + 1L,
                            firstActiveCol = freeze_col + 1L)

    if (is.null(col_widths)) {
      col_widths <- pmax(11, pmin(26,
        vapply(seq_len(ncol(df)), function(j) {
          w <- max(nchar(as.character(df[[j]])), na.rm = TRUE)
          max(w, nchar(names(df)[j]))
        }, numeric(1)) + 2))
    }
    openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(df)),
                            widths = col_widths)

    openxlsx::setRowHeights(wb, sheet, rows = r, heights = 22)
    invisible(r)
  }

  ev <- data.table::as.data.table(result$evidence)
  t1 <- data.frame(
    Trait                = ev$trait,
    Direction            = ev$direction,
    Mode                 = ev$mode,
    `Tested genes`       = ev$tested,
    `Superior genes`     = ev$superior,
    `N individuals`      = ev$n_obs,
    `h2`                 = .hbr_fmt(ev$h2, 3),
    `sigma2_g`           = .hbr_fmt(ev$sigma2_g, 2),
    `sigma2_e`           = .hbr_fmt(ev$sigma2_e, 2),
    `Feature share (w)`  = ifelse(is.na(ev$w_hat), "",
                                   .hbr_fmt(ev$w_hat, 3)),
    `N feature SNPs`     = ev$n_features,
    check.names = FALSE)
  add_table("1_Trait_summary",
    title = "Trait summary",
    subtitle = sprintf("Panel: %d individuals  |  MAF filter: %.3f  |  %d SNPs tested",
                       nrow(result$G$G), result$maf_min, result$G$n_total),
    df = t1, freeze_col = 1L)

  if (!is.null(result$cv_compare) &&
      any(!vapply(result$cv_compare, is.null, logical(1)))) {
    rows <- list()
    for (tr in names(result$cv_compare)) {
      cc <- result$cv_compare[[tr]]
      if (is.null(cc) || !nrow(cc)) next
      for (i in seq_len(nrow(cc))) {
        rows[[length(rows) + 1L]] <- data.frame(
          Trait          = tr,
          Contrast       = sprintf("%s - %s", cc$method[i], cc$baseline[i]),
          N              = cc$n[i],
          `d PA`         = .hbr_fmt(cc$d_PA[i], 4),
          `d PA lo`      = .hbr_fmt(cc$d_PA_lo[i], 4),
          `d PA hi`      = .hbr_fmt(cc$d_PA_hi[i], 4),
          `d RMSE`       = .hbr_fmt(cc$d_RMSE[i], 4),
          `d RMSE lo`    = .hbr_fmt(cc$d_RMSE_lo[i], 4),
          `d RMSE hi`    = .hbr_fmt(cc$d_RMSE_hi[i], 4),
          check.names = FALSE)
      }
    }
    if (length(rows)) {
      t2 <- do.call(rbind, rows)
      add_table("2_CV_evidence",
        title = "Cross-validation evidence (M2, M3 vs M0)",
        subtitle = paste0("Naive (M2) uses the pipeline's full-data superior set. ",
                          "Honest (M3) rediscoveries inside each outer training fold. ",
                          "95% bootstrap CI."),
        df = t2, freeze_col = 2L)
    }
  }

  ids <- result$fits[[1]]$ids
  trs <- names(result$fits)
  gebv <- data.table::data.table(id = ids)
  for (tr in trs) {
    f <- result$fits[[tr]]
    j <- match(ids, f$ids)
    dir <- result$directions[[tr]]
    sc  <- sqrt(f$sigma2_g)
    gebv[[paste0(tr, "_GEBV_M0")]]     <- f$g_M0[j]
    gebv[[paste0(tr, "_GEBV_M2")]]     <- f$g_M2[j]
    gebv[[paste0(tr, "_GEBV_scaled")]] <- if (identical(dir, "lower"))
      -f$g_M2[j] / sc else f$g_M2[j] / sc
    gebv[[paste0(tr, "_rel")]]         <- f$rel_M2[j]
  }
  idx_m <- result$index[match(ids, id)]
  gebv[, `:=`(index = idx_m$index, index_SE = idx_m$SE,
              rank = idx_m$rank, P_top15 = idx_m$P_top15)]
  data.table::setnames(gebv, "index", "MultiTrait_index")
  data.table::setnames(gebv, "index_SE", "MultiTrait_SE")
  data.table::setnames(gebv, "rank", "MultiTrait_rank")
  data.table::setnames(gebv, "P_top15", "P_top_15")

  gebv_out <- as.data.frame(gebv)
  num_cols <- setdiff(names(gebv_out), "id")
  for (cc in num_cols) gebv_out[[cc]] <- .hbr_fmt(gebv_out[[cc]], 4)

  add_table("3_GEBV_per_individual",
    title = "Genomic estimated breeding values per individual",
    subtitle = paste0("Raw (GEBV_M0 = GBLUP, GEBV_M2 = two-kernel) and ",
                       "scaled/direction-corrected (GEBV_scaled). ",
                       "Reliability per trait (rel)."),
    df = gebv_out, freeze_col = 1L)

  par <- data.table::as.data.table(result$parents)
  contrib_cols <- setdiff(names(result$index),
                           c("id","index","SE","rank","P_top15"))
  contrib_cols <- intersect(contrib_cols, names(par))

  t4 <- data.frame(
    Rank          = par$rank,
    ID            = par$id,
    `Multi-trait index` = .hbr_fmt(par$index, 3),
    `SE`          = .hbr_fmt(par$SE, 3),
    `P(top 15)`   = .hbr_fmt(par$P_top15, 3),
    `Haplotype coverage` = .hbr_fmt(par$haplo_cov, 3),
    `Mean relatedness`   = .hbr_fmt(par$mean_rel, 3),
    `Cross count` = par$n_crosses,
    Action        = par$action,
    check.names = FALSE)
  for (cc in contrib_cols) {
    t4[[paste0(cc, "_scaled")]] <- .hbr_fmt(par[[cc]], 3)
  }
  add_table("4_Parent_sheet",
    title = "Parent decision sheet",
    subtitle = sprintf("Top %d parents by multi-trait index. Contributions ",
                       nrow(par)) %+%
      "are per-trait GEBVs scaled by trait genetic SD and direction-corrected.",
    df = t4, freeze_col = 2L)

  cr <- data.table::as.data.table(result$crosses)
  t5 <- data.frame(
    Rank            = cr$rank,
    `Parent 1`      = cr$p1,
    `Parent 2`      = cr$p2,
    `Mid-parent index` = .hbr_fmt(cr$mid_index, 3),
    Relatedness     = .hbr_fmt(cr$relatedness, 3),
    `Haplotype coverage` = .hbr_fmt(cr$coverage, 3),
    Complementarity = .hbr_fmt(cr$complement, 3),
    Score           = .hbr_fmt(cr$score, 3),
    `N regions`     = cr$n_regions,
    check.names = FALSE)
  add_table("5_Cross_recommendations",
    title = "Recommended crosses",
    subtitle = paste0("Coverage and complementarity computed at the LD-region level. ",
                       "Relatedness ceiling applied."),
    df = t5, freeze_col = 3L)

  cf <- result$counterfactual
  cf_df <- data.frame(
    Category = c("Kept in both rankings",
                  "Entered under GBLUP-only",
                  "Exited under GBLUP-only",
                  "Top (haplotype-informed)",
                  "Top (GBLUP-only)"),
    Count = c(length(cf$kept), length(cf$entered_m0),
              length(cf$exited_m0), length(cf$top_m3), length(cf$top_m0)),
    IDs   = c(paste(cf$kept,       collapse = ", "),
              paste(cf$entered_m0, collapse = ", "),
              paste(cf$exited_m0,  collapse = ", "),
              paste(cf$top_m3,     collapse = ", "),
              paste(cf$top_m0,     collapse = ", ")),
    check.names = FALSE)
  add_table("6_Counterfactual",
    title = "Counterfactual -- what changes if haplotype features are ignored",
    subtitle = "Comparison of parent rankings with and without the feature model.",
    df = cf_df, freeze_col = 1L)

  if (!is.null(result$sup_sets) && length(result$sup_sets) > 0L) {
    all_sup <- data.table::rbindlist(
      lapply(names(result$sup_sets), function(tr) {
        s <- result$sup_sets[[tr]]
        if (is.null(s) || !nrow(s)) return(NULL)
        data.table::data.table(
          trait = tr,
          gene  = as.character(s$gene),
          haplotype = as.character(s$haplotype),
          nt_sequence = as.character(s$nt_sequence))
      }), fill = TRUE)

    multi <- all_sup[, .(
      n_traits  = data.table::uniqueN(trait),
      traits    = paste(sort(unique(trait)), collapse = ";"),
      n_unique_hap = data.table::uniqueN(nt_sequence),
      haplotypes   = paste(unique(haplotype), collapse = ";")),
      by = gene]
    multi <- multi[n_traits > 1L]
    multi[, concordant := n_unique_hap == 1L]

    if (nrow(multi)) {
      t7 <- data.frame(
        Gene           = multi$gene,
        `N traits`     = multi$n_traits,
        Traits         = multi$traits,
        Concordant     = ifelse(multi$concordant, "yes", "no"),
        Haplotypes     = multi$haplotypes,
        check.names = FALSE)
      add_table("7_Cross_trait_genes",
        title = "Genes superior in more than one trait",
        subtitle = "Concordant = same superior haplotype across all listed traits.",
        df = t7, freeze_col = 1L)
    }
  }

  prov <- data.frame(
    Field = c("Package", "Generated", "PLINK prefix",
               "Panel size", "Total SNPs", "SNPs after MAF",
               "MAF filter", "Relatedness ceiling",
               "Traits",
               "Weights",
               "Directions",
               "CV folds", "CV repeats", "CV n_perm",
               "Output directory"),
    Value = c("HaploBreedR / hbr_breeding_decision",
               format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
               result$plink_prefix,
               nrow(result$G$G), result$G$n_total, result$G$n_used,
               sprintf("%.3f", result$maf_min),
               sprintf("%.2f", result$parameters$relatedness_max),
               paste(names(result$traits), collapse = ", "),
               paste(sprintf("%s=%g", names(result$weights),
                             result$weights), collapse = ", "),
               paste(sprintf("%s=%s", names(result$directions),
                             result$directions), collapse = ", "),
               result$parameters$cv$folds,
               result$parameters$cv$repeats,
               result$parameters$cv$n_perm,
               result$output_dir),
    check.names = FALSE)
  add_table("8_Provenance",
    title = "Provenance and parameters",
    subtitle = "All values recorded at run time.",
    df = prov, freeze_col = 1L)

  src_rows <- lapply(names(result$traits), function(tr) {
    t <- result$traits[[tr]]
    data.frame(
      Trait           = tr,
      Source          = t$source,
      `Trait column`  = t$trait_col,
      Direction       = t$direction,
      Mode            = t$info$mode %||% "",
      `N tested`      = t$info$genes_tested %||% NA,
      `N superior`    = t$info$superior_count %||% NA,
      `n_perm`        = t$info$n_perm %||% NA,
      check.names = FALSE)
  })
  t9 <- do.call(rbind, src_rows)
  add_table("9_Sources",
    title = "Trait directory sources",
    subtitle = "Parameters read from each HaploBreedR pipeline output.",
    df = t9, freeze_col = 1L)

  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  message("Excel workbook written: ", path)
  invisible(path)
}


# ============================================================================
# MAIN EXPORTED FUNCTION
# ============================================================================

#' Haplotype-informed breeding decisions
#'
#' Given one or more HaploBreedR trait directories plus a PLINK prefix,
#' produces a multi-trait parent ranking, cross recommendations, and a
#' written report.
#'
#' @param trait_dirs Named character vector (names = short labels) or
#'   list of `list(name =, source =, trait_col =, direction =)`. Each
#'   source is a HaploBreedR trait directory (or path to its RDS).
#' @param plink_prefix PLINK prefix (no extension). `.bed/.bim/.fam` and
#'   `.raw` required.
#' @param output_dir Directory for saved objects and reports. Defaults to
#'   `breeding_decision_<YYYYMMDD>` in the current working directory.
#' @param weights Named numeric weights. Defaults to equal weights.
#' @param directions Override detected direction per trait.
#' @param maf_min MAF filter for the genome-wide GRM. Default 0.05.
#' @param relatedness_max Cross-pair relatedness ceiling. Default 0.30.
#' @param top_parents Number of parents to report. Default 20.
#' @param top_crosses Number of crosses to report. Default 15.
#' @param gate_params List of gate overrides passed to reconstruction and CV.
#' @param run_cv If TRUE, run honest nested CV per trait.
#' @param cv_folds,cv_repeats CV parameters.
#' @param cv_n_perm Permutation iterations inside the fold gate.
#' @param plink_path Path to the PLINK executable.
#' @param verbose Print progress.
#' @return A list of class `hgsi_result`.
#' @examples
#' \donttest{
#' # Requires a completed HaploBreedR pipeline run in the trait directory.
#' toy <- system.file("extdata", "toy", package = "HaploBreedR")
#' if (nzchar(toy) &&
#'     file.exists(file.path(toy, "pipeline_result_for_breeding.rds"))) {
#'   res <- hbr_breeding_decision(
#'     trait_dirs   = c(GY = toy),
#'     plink_prefix = file.path(toy, "toy"),
#'     output_dir   = tempdir())
#' }
#' }
#' @export
hbr_breeding_decision <- function(
    trait_dirs,
    plink_prefix,
    output_dir      = NULL,
    weights         = NULL,
    directions      = NULL,
    maf_min         = 0.05,
    relatedness_max = 0.30,
    top_parents     = 20L,
    top_crosses     = 15L,
    gate_params     = list(),
    run_cv          = FALSE,
    cv_folds        = 3L,
    cv_repeats      = 1L,
    cv_n_perm       = 50L,
    plink_path      = "plink",
    verbose         = TRUE) {

  t0 <- Sys.time()

  if (is.null(output_dir))
    output_dir <- file.path(getwd(),
      paste0("breeding_decision_", format(Sys.Date(), "%Y%m%d")))
  if (!dir.exists(output_dir))
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # ------------------------------------------------------------------
  # 1. Trait resolution
  # ------------------------------------------------------------------
  if (verbose) message("Step 1/8  Probing trait directories")
  traits <- .hbr_resolve_traits(trait_dirs, verbose = verbose)
  trait_names <- names(traits)
  .hbr_assert(length(trait_names) >= 1L, "Need at least one trait.")

  if (!is.null(directions)) {
    for (tr in names(directions)) {
      if (tr %in% trait_names) traits[[tr]]$direction <- directions[[tr]]
    }
  }
  for (tr in trait_names) {
    .hbr_assert(traits[[tr]]$direction %in% c("higher", "lower"),
               sprintf("Direction for '%s' not resolvable.", tr))
  }
  if (is.null(weights)) {
    weights <- stats::setNames(rep(1, length(trait_names)), trait_names)
  } else {
    .hbr_assert(all(trait_names %in% names(weights)),
               "`weights` must include all traits.")
  }
  directions_vec <- vapply(traits, function(t) t$direction, character(1))
  names(directions_vec) <- trait_names

  # ------------------------------------------------------------------
  # 2. Genotype + GRM
  # ------------------------------------------------------------------
  if (verbose) message("\nStep 2/8  Reading genotypes")
  geno  <- .hbr_read_plink(plink_prefix, plink_path, verbose = verbose)
  G_all <- geno$G
  ids   <- geno$ids
  snps  <- geno$snps

  if (verbose) message("\nStep 3/8  Building GRM")
  grm <- .hbr_build_grm(G_all, maf_min = maf_min, verbose = verbose)

  # ------------------------------------------------------------------
  # 3. Reference extraction
  # ------------------------------------------------------------------
  ref <- NULL
  ref_trait <- NULL
  for (tr in trait_names) {
    if (!is.null(traits[[tr]]$info$rds_path) &&
        file.exists(traits[[tr]]$info$rds_path)) {
      ref_trait <- tr; break
    }
  }
  if (!is.null(ref_trait)) {
    cache_path <- file.path(output_dir, "reference_cache.rds")
    ref <- .hbr_extract_reference(traits[[ref_trait]]$info,
                                   cache_path, verbose = verbose)
  } else {
    if (verbose)
      message("  No RDS available; using CSV fast path only")
  }

  # ------------------------------------------------------------------
  # 4. Per-trait processing
  # ------------------------------------------------------------------
  if (verbose) message("\nStep 4/8  Processing each trait")

  nt_seqs_lookup      <- list()
  hap_codes_per_trait <- list()
  sup_sets            <- list()
  ld_regions          <- list()
  fits                <- list()
  phenos              <- list()
  evidence            <- list()

  for (tr in trait_names) {
    if (verbose) message(sprintf("\n--- %s ---", tr))
    info   <- traits[[tr]]$info
    params <- .hbr_gate_params(info, overrides = gate_params)
    params$direction <- traits[[tr]]$direction

    ph <- .hbr_phenotypes(info, traits[[tr]]$trait_col, verbose = verbose)
    phenos[[tr]] <- ph

    sup_tbl <- .hbr_superior_table(info)
    ld_map  <- .hbr_ld_region_map(info)
    sup_sets[[tr]]   <- sup_tbl
    ld_regions[[tr]] <- ld_map
    if (verbose)
      message(sprintf("  %d superior genes; %d LD regions",
                      nrow(sup_tbl),
                      data.table::uniqueN(ld_map$region)))

    recon <- NULL
    if (!is.null(info$rds_path) ||
        file.exists(file.path(info$source,
                               "functional_snps_all_tested_genes.csv"))) {
      recon <- .hbr_reconstruct(info, tr, G_all, ref, params,
                                 verbose = verbose)
    }
    if (!is.null(recon)) {
      hap_codes_per_trait[[tr]] <- recon
      for (g in names(recon$nt_seqs))
        if (is.null(nt_seqs_lookup[[g]]))
          nt_seqs_lookup[[g]] <- recon$nt_seqs[[g]]
    }

    y <- rep(NA_real_, length(ids)); names(y) <- ids
    ph_c <- hbr_clean_id(ph$ID)
    idx <- match(ids, ph_c)
    ok <- !is.na(idx)
    y[ok] <- suppressWarnings(
      as.numeric(ph[[traits[[tr]]$trait_col]][idx[ok]]))

    fit <- .hbr_fit_panel(tr, y, ids, grm, G_all, recon, sup_tbl,
                           maf_min = maf_min)
    fits[[tr]] <- fit

    evidence[[tr]] <- data.table::data.table(
      trait      = tr,
      direction  = traits[[tr]]$direction,
      mode       = info$mode,
      tested     = info$genes_tested,
      superior   = info$superior_count,
      n_obs      = sum(!is.na(y)),
      sigma2_g   = fit$sigma2_g,
      sigma2_e   = fit$sigma2_e,
      h2         = fit$sigma2_g / (fit$sigma2_g + fit$sigma2_e),
      w_hat      = fit$w_hat,
      n_features = fit$n_features)
  }
  evidence <- data.table::rbindlist(evidence)

  cv_result <- NULL
  cv_compare <- NULL
  if (isTRUE(run_cv)) {
    if (verbose)
      message("\nStep 4b/8  Honest CV with fold-specific gate replay")
    if (!exists(".hbr_cv_one_honest", mode = "function")) {
      warning("run_cv = TRUE requires the CV helpers to be available. ",
              "Skipping CV.", call. = FALSE)
    } else {
      cv_result <- lapply(trait_names, function(tr) {
        tryCatch(
          .hbr_cv_one_honest(
            trait      = tr,
            fit        = fits[[tr]],
            trait_spec = traits[[tr]],
            geno       = geno,
            grm        = grm,
            recon      = hap_codes_per_trait[[tr]],
            sup_tbl    = sup_sets[[tr]],
            params     = .hbr_gate_params(traits[[tr]]$info,
                                           overrides = gate_params),
            k          = cv_folds,
            repeats    = cv_repeats,
            n_perm     = cv_n_perm,
            seed       = 1L,
            verbose    = verbose),
          error = function(e) {
            if (verbose)
              message("  CV failed for ", tr, ": ", conditionMessage(e))
            NULL
          })
      })
      names(cv_result) <- trait_names
      cv_compare <- lapply(cv_result, function(cv) {
        if (is.null(cv)) return(NULL)
        tryCatch(
          .hbr_cv_compare(cv, baseline = "M0", n_boot = 2000L, seed = 1L),
          error = function(e) NULL)
      })
      names(cv_compare) <- trait_names
    }
  }

  # ------------------------------------------------------------------
  # 5. Multi-trait index
  # ------------------------------------------------------------------
  if (verbose) message("\nStep 5/8  Multi-trait index")
  idx <- .hbr_index(fits, weights, directions_vec, scale_by_sd = TRUE)
  idx[, P_top15 := .hbr_rank_stability(index, SE, top_k = 15L, seed = 1L)]
  data.table::setorder(idx, -index)

  # ------------------------------------------------------------------
  # 6. Parent sheet
  # ------------------------------------------------------------------
  if (verbose) message("\nStep 6/8  Parent sheet")
  top_n   <- min(top_parents, nrow(idx))
  parents <- idx[seq_len(top_n)]

  parents[, mean_rel := vapply(id, function(x)
    mean(grm$G[x, setdiff(parents$id, x)]), numeric(1))]

  parents[, haplo_cov := vapply(id, function(p) {
    tot <- 0L; hit <- 0L
    for (tr in trait_names) {
      s <- merge(sup_sets[[tr]], ld_regions[[tr]], by = "gene",
                  all.x = TRUE)
      if (!nrow(s)) next
      s[is.na(region) | !nzchar(region),
        region := paste0("__g_", gene)]
      for (r in unique(s$region)) {
        sub <- s[region == r]
        nts <- unique(sub$nt_sequence)
        tot <- tot + 1L
        carry <- FALSE
        for (g in sub$gene) {
          seqs <- nt_seqs_lookup[[g]]
          if (is.null(seqs)) next
          s1 <- seqs[p]
          if (!is.na(s1) && s1 %in% nts) { carry <- TRUE; break }
        }
        if (carry) hit <- hit + 1L
      }
    }
    if (tot) hit / tot else NA_real_
  }, numeric(1))]

  parents[, action := data.table::fcase(
    P_top15 >= 0.85, "advance",
    P_top15 >= 0.50, "consider",
    default         = "hold")]

  # ------------------------------------------------------------------
  # 7. Crosses
  # ------------------------------------------------------------------
  if (verbose) message("\nStep 7/8  Cross recommendations")
  top_ids <- idx$id[seq_len(min(top_crosses + 5L, nrow(idx)))]
  crosses <- .hbr_pair_crosses(
    parents         = top_ids,
    idx             = idx,
    G               = grm$G,
    nt_seqs_lookup  = nt_seqs_lookup,
    sup_sets        = sup_sets,
    ld_regions      = ld_regions,
    relatedness_max = relatedness_max)
  crosses <- crosses[seq_len(min(top_crosses, .N))]

  parents[, n_crosses := vapply(id, function(x)
    sum(crosses$p1 == x | crosses$p2 == x), integer(1))]

  # ------------------------------------------------------------------
  # 8. Counterfactual
  # ------------------------------------------------------------------
  if (verbose) message("\nStep 8/8  Counterfactual (GBLUP only)")
  fits_m0 <- lapply(fits, function(f) { f$g_M2 <- f$g_M0; f })
  idx_m0  <- .hbr_index(fits_m0, weights, directions_vec,
                         scale_by_sd = TRUE)
  data.table::setorder(idx_m0, -index)

  top_k_used <- min(top_parents, nrow(idx), nrow(idx_m0))
  top_m3 <- idx$id[seq_len(top_k_used)]
  top_m0 <- idx_m0$id[seq_len(top_k_used)]
  counterfactual <- list(
    top_m3     = top_m3,
    top_m0     = top_m0,
    kept       = intersect(top_m3, top_m0),
    entered_m0 = setdiff(top_m0, top_m3),
    exited_m0  = setdiff(top_m3, top_m0))

  # ------------------------------------------------------------------
  # Compose result
  # ------------------------------------------------------------------
  result <- structure(
    list(
      traits         = traits,
      G              = grm,
      nt_seqs        = nt_seqs_lookup,
      hap_codes      = hap_codes_per_trait,
      sup_sets       = sup_sets,
      ld_regions     = ld_regions,
      fits           = fits,
      phenos         = phenos,
      evidence       = evidence,
      cv             = cv_result,
      cv_compare     = cv_compare,
      index          = idx,
      parents        = parents,
      crosses        = crosses,
      counterfactual = counterfactual,
      weights        = weights,
      directions     = directions_vec,
      output_dir     = output_dir,
      plink_prefix   = plink_prefix,
      maf_min        = maf_min,
      parameters     = list(
        relatedness_max = relatedness_max,
        top_parents     = top_parents,
        top_crosses     = top_crosses,
        gate_params     = gate_params,
        run_cv          = run_cv,
        cv              = list(folds = cv_folds,
                                repeats = cv_repeats,
                                n_perm = cv_n_perm))),
    class = c("hgsi_result", "list"))

  # ------------------------------------------------------------------
  # Write outputs
  # ------------------------------------------------------------------
  result$report <- tryCatch(
    .hbr_report(result, verbose = verbose),
    error = function(e) {
      warning("Report generation failed: ",
              conditionMessage(e), call. = FALSE)
      NULL
    })
  if (!is.null(result$report))
    writeLines(result$report,
               file.path(output_dir, "decision_report.md"))

  xlsx_path <- tryCatch(
    .hbr_write_workbook(result,
      path = file.path(output_dir,
        sprintf("breeding_decision_%s.xlsx",
                format(Sys.Date(), "%Y%m%d")))),
    error = function(e) {
      warning("Workbook export failed: ",
              conditionMessage(e), call. = FALSE)
      NULL
    })
  if (!is.null(xlsx_path)) result$workbook <- xlsx_path

  tryCatch(.hbr_write_gebv_table(result, output_dir),
           error = function(e)
             warning("GEBV table failed: ",
                     conditionMessage(e), call. = FALSE))

  if (!is.null(cv_result)) {
    for (tr in names(cv_result)) {
      if (is.null(cv_result[[tr]])) next
      tryCatch(
        data.table::fwrite(
          cv_result[[tr]]$predictions,
          file.path(output_dir,
                    paste0("cv_oof_gebv_", tr, ".csv"))),
        error = function(e) NULL)
    }
  }

  saveRDS(result, file.path(output_dir, "decision_result.rds"))

  if (verbose) {
    .hbr_log(t0, "hbr_breeding_decision", verbose)
    message(sprintf("\nOutputs written to %s", output_dir))
    message(sprintf("  Output size: %.1f MB",
                    .hbr_dir_size(output_dir) / 1e6))
    message("  decision_report.md")
    message("  breeding_decision_<date>.xlsx")
    message("  gebv_per_individual.csv")
    message("  decision_result.rds")
    if (!is.null(cv_result)) message("  cv_oof_gebv_<trait>.csv")
  }
  invisible(result)
}


#' Save the breeding RDS from an existing haplobreed_analysis() result
#'
#' Recovers the `pipeline_result_for_breeding.rds` (for example after an
#' in-run save failed) without re-running the pipeline, as long as the
#' result object is still in the R session.
#'
#' @param result Object returned by [haplobreed_analysis()].
#' @param output_dir Directory to write to. Defaults to the run's output
#'   directory recorded in `result$pipeline_summary$Output_Dir`.
#' @param rds_compress One of `"fast"`, `"none"`, `"xz"`.
#' @param rds_reference One of `"embedded"`, `"separate"`, `"none"`.
#' @param verbose Print progress.
#' @return Invisibly, the output of [`.hbr_write_breeding_rds()`].
#' @examples
#' result <- list(
#'   superior_results = list(),
#'   geno_mat         = matrix(c(0L, 1L, 2L, 0L), nrow = 2, ncol = 2,
#'                             dimnames = list(c("i1", "i2"),
#'                                             c("s1", "s2"))),
#'   reference        = list(genome = "toy", annotation = "CUSTOM"),
#'   reference_paths  = list(),
#'   plink_prefix     = "toy",
#'   pipeline_summary = data.frame(Output_Dir = tempdir(),
#'                                 stringsAsFactors = FALSE))
#' out <- hbr_save_breeding_rds(result, output_dir = tempdir())
#' file.exists(out$rds)
#' @export
hbr_save_breeding_rds <- function(result, output_dir = NULL,
                                  rds_compress  = c("fast", "none", "xz"),
                                  rds_reference = c("embedded", "separate", "none"),
                                  verbose = TRUE) {
  if (is.null(output_dir))
    output_dir <- result$pipeline_summary$Output_Dir[1]
  if (is.null(output_dir) || is.na(output_dir) || !dir.exists(output_dir))
    stop("Supply an existing `output_dir`.", call. = FALSE)
  .hbr_write_breeding_rds(
    superior_results = result$superior_results,
    geno_mat         = result$geno_mat,
    reference        = result$reference,
    ref_paths        = result$reference_paths,
    plink_prefix     = result$plink_prefix,
    pipeline_summary = result$pipeline_summary,
    output_dir       = output_dir,
    rds_compress     = match.arg(rds_compress),
    rds_reference    = match.arg(rds_reference),
    verbose          = verbose)
}

#' Print an hgsi_result object
#'
#' Compact one-screen summary of a haplotype-informed breeding decision:
#' trait count, ranked parents, recommended crosses, and output path.
#'
#' @param x An object of class `hgsi_result`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @examples
#' res <- structure(
#'   list(traits   = list(GY = list(name = "GY")),
#'        parents  = data.table::data.table(
#'          rank = 1L, id = "IND01", index = 0.5, SE = 0.1,
#'          P_top15 = 0.9, haplo_cov = 0.8, mean_rel = 0.1,
#'          n_crosses = 2L, action = "advance"),
#'        crosses  = data.table::data.table(
#'          rank = 1L, p1 = "IND01", p2 = "IND02",
#'          mid_index = 0.4, relatedness = 0.1,
#'          complement = 0.5, score = 1.4),
#'        evidence = data.table::data.table(
#'          trait = "GY", direction = "higher",
#'          mode = "GWAS", tested = 10L, superior = 2L,
#'          n_obs = 100L, sigma2_g = 0.5, sigma2_e = 0.5,
#'          h2 = 0.5, w_hat = 0.3, n_features = 5L),
#'        output_dir = tempdir()),
#'   class = "hgsi_result")
#' print(res)
#'
#' @export
print.hgsi_result <- function(x, ...) {
  cat("<hgsi_result>\n")
  cat(sprintf("  traits      : %s\n",
              paste(names(x$traits), collapse = ", ")))
  cat(sprintf("  parents     : %d ranked, %d to advance\n",
              nrow(x$parents),
              sum(x$parents$action == "advance")))
  cat(sprintf("  crosses     : %d recommended\n", nrow(x$crosses)))
  cat(sprintf("  output      : %s\n", x$output_dir))
  invisible(x)
}

#' Summarise an hgsi_result object
#'
#' Prints the `hgsi_result` header, the per-trait evidence table, the top
#' ten parents, and the top ten recommended crosses.
#'
#' @param object An object of class `hgsi_result`.
#' @param ... Ignored.
#' @return `object`, invisibly.
#' @examples
#' res <- structure(
#'   list(traits   = list(GY = list(name = "GY")),
#'        parents  = data.table::data.table(
#'          rank = 1L, id = "IND01", index = 0.5, SE = 0.1,
#'          P_top15 = 0.9, haplo_cov = 0.8, mean_rel = 0.1,
#'          n_crosses = 2L, action = "advance"),
#'        crosses  = data.table::data.table(
#'          rank = 1L, p1 = "IND01", p2 = "IND02",
#'          mid_index = 0.4, relatedness = 0.1,
#'          complement = 0.5, score = 1.4),
#'        evidence = data.table::data.table(
#'          trait = "GY", direction = "higher",
#'          mode = "GWAS", tested = 10L, superior = 2L,
#'          n_obs = 100L, sigma2_g = 0.5, sigma2_e = 0.5,
#'          h2 = 0.5, w_hat = 0.3, n_features = 5L),
#'        output_dir = tempdir()),
#'   class = "hgsi_result")
#' summary(res)
#'
#' @export
summary.hgsi_result <- function(object, ...) {
  print(object)
  cat("\nEvidence:\n")
  print(object$evidence)
  cat("\nTop parents:\n")
  print(object$parents[seq_len(min(10, nrow(object$parents))),
                       .(rank, id,
                         index     = round(index, 3),
                         P_top15   = round(P_top15, 3),
                         haplo_cov = round(haplo_cov, 3),
                         action)])
  cat("\nTop crosses:\n")
  print(object$crosses[seq_len(min(10, nrow(object$crosses))),
                       .(rank, p1, p2,
                         mid_index   = round(mid_index, 3),
                         relatedness = round(relatedness, 3),
                         complement  = round(complement, 3),
                         score       = round(score, 3))])
  invisible(object)
}

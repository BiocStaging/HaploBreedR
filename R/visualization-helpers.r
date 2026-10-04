# ============================================================================
# R/visualization-helpers.r
# Shared helpers for the plotting layer.
# ============================================================================


# ----------------------------------------------------------------------------
# 1. Column aliasing
# ----------------------------------------------------------------------------

#' Column-name alias registry
#'
#' Alternative names recognised for each canonical column by
#' [`.hbr_resolve_col()`].
#'
#' @noRd
.hbr_col_aliases <- list(
  gene_id            = c("MSU_Gene_ID", "RAPDB_Gene_ID", "MSU7_Gene_ID",
                         "Gene_ID", "Gene", "gene_id", "gene",
                         "Candidate_Gene_ID"),
  gene_function      = c("MSU_Gene_Function", "RAPDB_Gene_Function",
                         "Gene_Function", "gene_function"),
  chromosome         = c("Chr", "chr", "CHR", "Chromosome", "chromosome"),
  gene_start         = c("Gene_Start", "gene_start", "Start", "start"),
  gene_end           = c("Gene_End", "gene_end", "End", "end"),
  p_value            = c("P_value", "p_value", "P", "p", "PVAL", "Pvalue"),
  permutation_p      = c("Permutation_P", "permutation_p", "Perm_P"),
  adjusted_p         = c("Adjusted_P", "adjusted_p", "FDR_P", "q_value"),
  superior_haplotype = c("Superior_Haplotype", "superior_haplotype",
                         "Superior_Hap"),
  is_superior        = c("Is_Superior", "is_superior", "Superior"),
  cohens_d           = c("Cohens_d", "cohens_d", "Cohen_D"),
  improvement_pct    = c("Improvement_Pct", "improvement_pct"),
  superior_mean      = c("Superior_Mean", "superior_mean", "Superior_Hap_Mean"),
  descriptive_letter = c("Descriptive_Letter", "Letter", "letter",
                         "Significance_Group", "Duncan_Group",
                         "Tukey_Group"),
  gwas_snp           = c("GWAS_SNP", "gwas_snp", "Lead_SNP", "SNP"),
  gwas_position      = c("GWAS_Position", "GWAS_POS", "SNP_Position",
                         "Position", "Candidate_GWAS_SNP_Position"))

#' Resolve a canonical column name to an actual column
#'
#' @param dt A data.table.
#' @param canonical Canonical role name.
#' @param required Stop if the column is missing.
#'
#' @return The resolved column name or `NA_character_`.
#'
#' @noRd
.hbr_resolve_col <- function(dt, canonical, required = FALSE) {
  aliases <- if (canonical %in% names(.hbr_col_aliases))
    .hbr_col_aliases[[canonical]] else canonical

  avail_lower <- tolower(trimws(names(dt)))
  for (cand in c(canonical, aliases)) {
    hit <- which(avail_lower == tolower(cand))[1]
    if (!is.na(hit)) return(names(dt)[hit])
  }
  if (required)
    stop("Missing column (any of): ",
         paste(c(canonical, aliases), collapse = ", "),
         "\n  Available: ", paste(names(dt), collapse = ", "),
         call. = FALSE)
  NA_character_
}


# ----------------------------------------------------------------------------
# 2. Gene-list resolver
# ----------------------------------------------------------------------------

#' Resolve a `gene_list` argument
#'
#' Accepts a character vector of gene IDs, a summary CSV path, a
#' directory containing the pipeline outputs, or an in-memory
#' `data.frame`. Applies superior-gene filtering and optional sorting
#' and capping by `n_gene`.
#'
#' @param gene_list Gene IDs, CSV, directory, or data.frame.
#' @param summary_csv Summary CSV path.
#' @param output_dir Directory containing the pipeline outputs.
#' @param n_gene Maximum number of genes to keep.
#' @param sort_by Columns to sort by.
#' @param verbose Print progress.
#'
#' @return A character vector of gene IDs.
#'
#' @noRd
.hbr_resolve_genes <- function(gene_list   = NULL,
                               summary_csv = NULL,
                               output_dir  = NULL,
                               n_gene      = NULL,
                               sort_by     = c("Cohens_d", "Superior_Mean"),
                               verbose     = TRUE) {

  if (!is.null(n_gene)) {
    n_gene <- suppressWarnings(as.integer(n_gene))
    if (is.na(n_gene) || n_gene < 1L)
      stop("`n_gene` must be NULL or a positive integer.", call. = FALSE)
  }

  parse_sort_spec <- function(spec) {
    if (is.null(spec) || !length(spec)) return(NULL)
    if (!is.null(names(spec)) && any(nzchar(names(spec)))) {
      cols <- names(spec)
      dirs <- tolower(as.character(spec))
      dirs[!dirs %in% c("asc", "desc")] <- "desc"
      ords <- ifelse(dirs == "asc", 1L, -1L)
    } else {
      cols <- as.character(spec)
      ords <- ifelse(
        grepl("^P_|_P$|P_value|p_value|Permutation|Adjusted_P|FDR|q_value",
              cols, ignore.case = TRUE), 1L, -1L)
    }
    list(cols = cols, ords = ords)
  }
  sort_spec <- parse_sort_spec(sort_by)

  apply_filters <- function(s) {
    n_before <- nrow(s)
    sup_col <- .hbr_resolve_col(s, "is_superior", required = FALSE)
    if (!is.na(sup_col)) {
      keep <- as.logical(s[[sup_col]]) %in% TRUE
      s <- s[keep]
    } else if (verbose) {
      message("  No Is_Superior column found; using all rows.")
    }
    if (verbose && nrow(s) < n_before)
      message("  Filtered ", n_before, " -> ", nrow(s),
              " superior gene(s).")
    s
  }

  rank_and_cap <- function(dt, gene_col) {
    ids <- unique(as.character(dt[[gene_col]]))
    ids <- ids[!is.na(ids) & nzchar(ids)]

    if (!is.null(sort_spec)) {
      actual_cols <- vapply(sort_spec$cols, function(cc) {
        .hbr_resolve_col(dt, cc, required = FALSE)
      }, character(1))
      keep <- !is.na(actual_cols)
      actual_cols <- actual_cols[keep]
      actual_ords <- sort_spec$ords[keep]

      if (length(actual_cols)) {
        sub <- data.table::data.table(.gene = as.character(dt[[gene_col]]))
        for (k in seq_along(actual_cols))
          sub[[actual_cols[k]]] <- suppressWarnings(
            as.numeric(dt[[actual_cols[k]]]))
        sub <- sub[!is.na(.gene) & nzchar(.gene)]
        if (nrow(sub)) {
          notall_na <- Reduce(`|`, lapply(actual_cols, function(cn)
            !is.na(sub[[cn]])))
          sub <- sub[notall_na]
        }
        sub <- unique(sub, by = ".gene")

        if (nrow(sub)) {
          data.table::setorderv(sub, actual_cols, order = actual_ords,
                                na.last = TRUE)
          ids <- as.character(sub$.gene)
        }
      } else if (verbose) {
        message("  None of sort_by (",
                paste(sort_spec$cols, collapse = ", "),
                ") found in summary; using source order.")
      }
    }

    if (!is.null(n_gene)) ids <- utils::head(ids, n_gene)
    ids
  }

  if (is.data.frame(gene_list)) {
    gc <- .hbr_resolve_col(gene_list, "gene_id", required = FALSE)
    if (is.na(gc)) return(character(0))
    gene_list <- apply_filters(data.table::as.data.table(gene_list))
    return(rank_and_cap(gene_list, gc))
  }

  if (is.character(gene_list) && length(gene_list) >= 1L) {
    if (length(gene_list) == 1L) {
      p <- gene_list
      if (dir.exists(p)) {
        f <- file.path(p, "publication_table_functional.csv")
        if (file.exists(f)) {
          if (verbose) message("Reading summary CSV: ", f)
          s  <- data.table::fread(f)
          gc <- .hbr_resolve_col(s, "gene_id", required = TRUE)
          s  <- apply_filters(s)
          return(rank_and_cap(s, gc))
        }
      }
      if (file.exists(p) && grepl("\\.csv$", p, ignore.case = TRUE)) {
        if (verbose) message("Reading summary CSV: ", p)
        s  <- data.table::fread(p)
        gc <- .hbr_resolve_col(s, "gene_id", required = TRUE)
        s  <- apply_filters(s)
        return(rank_and_cap(s, gc))
      }
    }
    ids <- unique(as.character(gene_list))
    if (!is.null(n_gene)) ids <- utils::head(ids, n_gene)
    return(ids)
  }

  target <- summary_csv %||% output_dir
  if (!is.null(target)) {
    if (dir.exists(target)) {
      f <- file.path(target, "publication_table_functional.csv")
      if (!file.exists(f))
        stop("No publication_table_functional.csv in ", target,
             call. = FALSE)
      target <- f
    }
    if (file.exists(target) && grepl("\\.csv$", target, ignore.case = TRUE)) {
      s  <- data.table::fread(target)
      gc <- .hbr_resolve_col(s, "gene_id", required = TRUE)
      s  <- apply_filters(s)
      return(rank_and_cap(s, gc))
    }
  }
  character(0)
}


# ----------------------------------------------------------------------------
# 3. Genotype format detection and subset readers
# ----------------------------------------------------------------------------

#' Detect genotype file format
#'
#' @param geno_file Path or prefix.
#'
#' @return One of `"vcf"`, `"raw"`, `"csv"`, `"hapmap"`, `"plink"`.
#'
#' @noRd
.hbr_detect_geno_format <- function(geno_file) {
  if (missing(geno_file) || is.null(geno_file)) return(NA_character_)
  low <- tolower(geno_file)
  if (grepl("\\.vcf(\\.gz)?$", low))    return("vcf")
  if (grepl("\\.raw$", low))            return("raw")
  if (grepl("\\.csv$", low))            return("csv")
  if (grepl("\\.hmp(\\.txt)?$", low))   return("hapmap")
  if (grepl("\\.(bed|bim|fam)$", low))  return("plink")
  if (file.exists(paste0(geno_file, ".raw")))     return("raw")
  if (file.exists(paste0(geno_file, ".csv")))     return("csv")
  if (file.exists(paste0(geno_file, ".hmp.txt"))) return("hapmap")
  if (file.exists(paste0(geno_file, ".bed")))     return("plink")
  if (file.exists(paste0(geno_file, ".vcf")))     return("vcf")
  if (file.exists(paste0(geno_file, ".vcf.gz")))  return("vcf")
  stop("Cannot auto-detect genotype format for: ", geno_file, call. = FALSE)
}

#' Normalise a genotype file path for the given format
#'
#' @param geno_file Path or prefix.
#' @param geno_format Format label.
#'
#' @return A normalised path.
#'
#' @noRd
.hbr_normalize_geno_path <- function(geno_file, geno_format) {
  if (is.null(geno_file)) return(geno_file)
  low <- tolower(geno_file)
  switch(geno_format,
         raw = if (grepl("\\.raw$", low)) geno_file
               else paste0(geno_file, ".raw"),
         csv = if (grepl("\\.csv$", low)) geno_file
               else paste0(geno_file, ".csv"),
         vcf = if (grepl("\\.vcf(\\.gz)?$", low)) geno_file
               else if (file.exists(paste0(geno_file, ".vcf.gz")))
                 paste0(geno_file, ".vcf.gz")
               else paste0(geno_file, ".vcf"),
         hapmap = if (grepl("\\.(hmp|hmp\\.txt|txt)$", low)) geno_file
                  else if (file.exists(paste0(geno_file, ".hmp.txt")))
                    paste0(geno_file, ".hmp.txt")
                  else paste0(geno_file, ".txt"),
         plink = sub("\\.(bed|bim|fam)$", "", geno_file, ignore.case = TRUE),
         geno_file)
}

#' Read a subset of a PLINK .raw file
#'
#' @param path PLINK .raw file path.
#' @param req_snps Required SNPs.
#' @param req_ids Required individuals.
#' @param bim_path Optional .bim path for column renaming.
#'
#' @return A numeric matrix or `NULL`.
#'
#' @noRd
.hbr_read_geno_raw_subset <- function(path, req_snps, req_ids = NULL,
                                      bim_path = NULL) {
  if (!file.exists(path))
    stop("Genotype file not found: ", path, call. = FALSE)
  g <- data.table::fread(path)
  ids <- as.character(g$IID)
  geno_cols <- setdiff(names(g), c("FID", "IID", "PAT", "MAT", "SEX",
                                    "PHENOTYPE"))
  raw_snps <- sub("_[ACGTDI]+$", "", geno_cols, ignore.case = TRUE)
  if (!is.null(bim_path) && file.exists(bim_path) &&
      !any(raw_snps %in% req_snps)) {
    bim <- data.table::fread(bim_path, header = FALSE)
    if (nrow(bim) == length(geno_cols))
      raw_snps <- as.character(bim[[2]])
  }
  keep <- which(raw_snps %in% req_snps)
  if (!length(keep)) return(NULL)
  m <- as.matrix(g[, ..geno_cols][, keep, with = FALSE])
  storage.mode(m) <- "numeric"
  rownames(m) <- ids
  colnames(m) <- raw_snps[keep]
  if (!is.null(req_ids)) {
    idx <- match(req_ids, rownames(m))
    k <- !is.na(idx)
    if (any(k)) m <- m[idx[k], , drop = FALSE]
  }
  m
}

#' Read a subset of a CSV genotype matrix
#'
#' @param path CSV path with an `id` column.
#' @param req_snps Required SNPs.
#' @param req_ids Required individuals.
#'
#' @return A numeric matrix or `NULL`.
#'
#' @noRd
.hbr_read_geno_csv_subset <- function(path, req_snps, req_ids = NULL) {
  if (!file.exists(path))
    stop("Genotype file not found: ", path, call. = FALSE)
  g <- data.table::fread(path, na.strings = c("", "NA"))
  ids <- as.character(g$id)
  snp_cols <- setdiff(names(g), "id")
  keep <- intersect(req_snps, snp_cols)
  if (!length(keep)) return(NULL)
  m <- as.matrix(g[, ..keep])
  storage.mode(m) <- "numeric"
  rownames(m) <- ids
  if (!is.null(req_ids)) {
    idx <- match(req_ids, rownames(m))
    k <- !is.na(idx)
    if (any(k)) m <- m[idx[k], , drop = FALSE]
  }
  m
}

#' Read a subset of SNPs via PLINK
#'
#' @param geno_file VCF or PLINK prefix.
#' @param plink_prefix PLINK prefix for `--bfile`.
#' @param plink_path PLINK executable path.
#' @param req_snps Required SNPs.
#' @param req_ids Required individuals.
#'
#' @return A numeric matrix or `NULL`.
#'
#' @noRd
.hbr_read_geno_via_plink <- function(geno_file, plink_prefix, plink_path,
                                      req_snps, req_ids = NULL) {
  snp_f <- tempfile(fileext = ".txt")
  writeLines(req_snps, snp_f)
  on.exit(unlink(snp_f), add = TRUE)
  out_pref <- tempfile("hbrsub_")
  on.exit({
    for (ext in c(".raw", ".bim", ".fam", ".bed", ".log", ".nosex"))
      if (file.exists(paste0(out_pref, ext)))
        unlink(paste0(out_pref, ext))
  }, add = TRUE)
  base_args <- if (!is.null(plink_prefix))
    c("--bfile", shQuote(plink_prefix))
  else c("--vcf", shQuote(geno_file))
  system2(plink_path,
          c(base_args, "--extract", shQuote(snp_f),
            "--recode", "A", "--out", shQuote(out_pref)),
          stdout = TRUE, stderr = TRUE)
  raw_f <- paste0(out_pref, ".raw")
  if (!file.exists(raw_f)) return(NULL)
  .hbr_read_geno_raw_subset(raw_f, req_snps, req_ids)
}

#' Dispatch to the right genotype subset reader
#'
#' @param geno_file Genotype path or prefix.
#' @param geno_format Format label.
#' @param req_snps Required SNPs.
#' @param req_ids Required individuals.
#' @param plink_path PLINK executable path.
#' @param verbose Print progress.
#'
#' @return A numeric matrix or `NULL`.
#'
#' @noRd
.hbr_read_geno_subset <- function(geno_file, geno_format,
                                   req_snps, req_ids = NULL,
                                   plink_path = "plink",
                                   verbose = FALSE) {
  geno_file <- .hbr_normalize_geno_path(geno_file, geno_format)
  if (verbose) {
    message("  resolved genotype file: ", basename(geno_file))
    message("  exists? ", file.exists(geno_file))
  }
  switch(geno_format,
         raw = {
           bim <- sub("\\.raw$", ".bim", geno_file)
           .hbr_read_geno_raw_subset(
             geno_file, req_snps, req_ids,
             bim_path = if (file.exists(bim)) bim else NULL)
         },
         csv = .hbr_read_geno_csv_subset(geno_file, req_snps, req_ids),
         vcf = .hbr_read_geno_via_plink(geno_file, NULL, plink_path,
                                        req_snps, req_ids),
         plink = .hbr_read_geno_via_plink(NULL, geno_file, plink_path,
                                          req_snps, req_ids),
         hapmap = {
           g <- data.table::fread(geno_file, header = TRUE)
           if (!"rs#" %in% names(g)) return(NULL)
           iupac <- c(R = "AG", Y = "CT", S = "GC", W = "AT", K = "GT",
                      M = "AC", A = "AA", C = "CC", G = "GG", T = "TT")
           dec <- function(call, ref, alt) {
             call <- toupper(trimws(call))
             two <- if (nchar(call) == 2L) call
             else if (call %in% names(iupac)) iupac[[call]]
             else NA_character_
             if (is.na(two) || grepl("[^ACGT]", two)) return(NA_real_)
             b <- strsplit(two, "")[[1]]
             if (!all(b %in% c(ref, alt))) return(NA_real_)
             sum(b == alt)
           }
           core <- c("rs#", "alleles", "chrom", "pos", "strand",
                     "assembly#", "center", "protLSID", "assayLSID",
                     "panelLSID", "QCcode")
           ind_cols <- setdiff(names(g), core)
           g <- g[g[["rs#"]] %in% req_snps]
           if (!nrow(g)) return(NULL)
           ref <- sub("/.*$", "", g$alleles)
           alt <- sub("^[^/]*/", "", g$alleles)
           m <- matrix(NA_real_, length(ind_cols), nrow(g),
                       dimnames = list(ind_cols, g[["rs#"]]))
           for (j in seq_along(ind_cols)) {
             col <- g[[ind_cols[j]]]
             m[j, ] <- mapply(dec, col, ref, alt, USE.NAMES = FALSE)
           }
           if (!is.null(req_ids)) {
             idx <- match(req_ids, rownames(m))
             k <- !is.na(idx)
             if (any(k)) m <- m[idx[k], , drop = FALSE]
           }
           m
         },
         stop("Unknown geno_format: ", geno_format, call. = FALSE))
}


# ----------------------------------------------------------------------------
# 4. Per-gene loader
# ----------------------------------------------------------------------------

#' Load per-gene files
#'
#' Loads `_hap_map.csv`, `_groups.csv`, `_functional_snps.csv`, and
#' `_individuals.csv` for a gene and returns them as a list.
#'
#' @param gene_id Gene ID.
#' @param per_gene_dir Directory with per-gene files.
#' @param summary_csv Summary CSV path.
#' @param genotype_m Optional genotype matrix.
#' @param geno_file Optional genotype file path.
#' @param geno_format Genotype file format.
#' @param plink_path PLINK executable path.
#' @param trait_col Trait column name.
#' @param verbose Print progress.
#'
#' @return A list with `hap_map`, `groups`, `func_snps`,
#'   `individuals_df`, `geno`, `meta`.
#'
#' @noRd
.hbr_load_gene_pieces <- function(gene_id,
                                   per_gene_dir,
                                   summary_csv = NULL,
                                   genotype_m  = NULL,
                                   geno_file   = NULL,
                                   geno_format = NULL,
                                   plink_path  = "plink",
                                   trait_col   = "NT",
                                   verbose     = FALSE) {

  find_f <- function(suffix, required = TRUE) {
    f <- file.path(per_gene_dir, paste0(gene_id, suffix))
    if (!file.exists(f)) {
      if (required) stop("Missing per-gene file: ", f, call. = FALSE)
      return(NULL)
    }
    f
  }

  hap_map   <- data.table::fread(find_f("_hap_map.csv"))
  groups    <- data.table::fread(find_f("_groups.csv"))
  func_snps <- data.table::fread(find_f("_functional_snps.csv"))

  ind_file <- find_f("_individuals.csv", required = FALSE)
  individuals_df <- NULL
  if (!is.null(ind_file)) {
    ind <- data.table::fread(ind_file, na.strings = c("", "NA"))
    id_col  <- intersect(c("ID", "id", "IID"), names(ind))[1]
    hap_col <- intersect(c("Haplotype", "haplotype"), names(ind))[1]
    tr_col  <- intersect(c(trait_col, "trait_value", "Trait"),
                         names(ind))[1]
    if (!is.na(id_col) && !is.na(hap_col)) {
      individuals_df <- data.table::data.table(
        ID        = as.character(ind[[id_col]]),
        Haplotype = as.character(ind[[hap_col]]),
        .trait    = if (!is.na(tr_col)) as.numeric(ind[[tr_col]])
                    else NA_real_)
      data.table::setnames(individuals_df, ".trait", trait_col)
    }
  }

  meta <- list()
  if (!is.null(summary_csv) && file.exists(summary_csv)) {
    s <- tryCatch(data.table::fread(summary_csv), error = function(e) NULL)
    if (!is.null(s)) {
      gc <- .hbr_resolve_col(s, "gene_id", required = FALSE)
      if (!is.na(gc)) {
        target_id <- as.character(gene_id)
        row <- s[as.character(s[[gc]]) == target_id]
        if (nrow(row)) {
          test_col <- .hbr_resolve_col(s, "test",            required = FALSE)
          p_col    <- .hbr_resolve_col(s, "p_value",         required = FALSE)
          pp_col   <- .hbr_resolve_col(s, "permutation_p",   required = FALSE)
          ap_col   <- .hbr_resolve_col(s, "adjusted_p",      required = FALSE)
          cd_col   <- .hbr_resolve_col(s, "cohens_d",        required = FALSE)
          imp_col  <- .hbr_resolve_col(s, "improvement_pct", required = FALSE)
          sup_col  <- .hbr_resolve_col(s, "superior_haplotype",
                                        required = FALSE)
          chr_col  <- .hbr_resolve_col(s, "chromosome",      required = FALSE)
          get1 <- function(cn, coerce = identity) {
            if (is.na(cn)) return(NA)
            tryCatch(coerce(row[[cn]][1]), error = function(e) NA)
          }
          meta$test          <- get1(test_col, as.character)
          meta$p_value       <- get1(p_col,    as.numeric)
          meta$permutation_p <- get1(pp_col,   as.numeric)
          meta$adjusted_p    <- get1(ap_col,   as.numeric)
          meta$cohens_d      <- get1(cd_col,   as.numeric)
          meta$improvement   <- get1(imp_col,  as.numeric)
          meta$superior_hap  <- get1(sup_col,  as.character)
          meta$chromosome    <- get1(chr_col,  as.character)
        }
      }
    }
  }

  geno_mat <- NULL
  if (!is.null(genotype_m)) {
    m <- genotype_m
    keep <- intersect(as.character(func_snps$SNP), colnames(m))
    if (length(keep)) {
      m <- m[, keep, drop = FALSE]
      if (!is.null(individuals_df)) {
        idx <- match(individuals_df$ID, rownames(m))
        k <- !is.na(idx)
        if (any(k)) m <- m[idx[k], , drop = FALSE]
      }
      geno_mat <- m
    }
  } else if (!is.null(geno_file)) {
    fmt <- tryCatch(geno_format %||% .hbr_detect_geno_format(geno_file),
                    error = function(e) NA_character_)
    if (!is.na(fmt)) {
      req_ids <- if (!is.null(individuals_df)) individuals_df$ID else NULL
      geno_mat <- tryCatch(
        .hbr_read_geno_subset(geno_file, fmt,
                              as.character(func_snps$SNP),
                              req_ids, plink_path, verbose = verbose),
        error = function(e) NULL)
    }
  }

  list(hap_map        = hap_map,
       groups         = groups,
       func_snps      = func_snps,
       individuals_df = individuals_df,
       geno           = geno_mat,
       meta           = meta)
}


# ----------------------------------------------------------------------------
# 5. Amino-acid label expansion
# ----------------------------------------------------------------------------

#' Expand one-letter amino-acid labels to three-letter
#'
#' @param labels Character vector of labels like `"R65G"`.
#'
#' @return Character vector, same length.
#'
#' @noRd
.hbr_expand_aa_labels <- function(labels) {
  map <- c(A = "Ala", R = "Arg", N = "Asn", D = "Asp", C = "Cys",
           Q = "Gln", E = "Glu", G = "Gly", H = "His", I = "Ile",
           L = "Leu", K = "Lys", M = "Met", F = "Phe", P = "Pro",
           S = "Ser", T = "Thr", W = "Trp", Y = "Tyr", V = "Val",
           `*` = "Ter")
  vapply(labels, function(s) {
    if (is.na(s) || !nzchar(s)) return(s)
    m <- regmatches(s, regexec("^([A-Za-z*])([0-9]+)([A-Za-z*])$",
                                s))[[1]]
    if (length(m) == 4) {
      a <- map[toupper(m[2])]; b <- map[toupper(m[4])]
      if (!is.na(a) && !is.na(b)) return(paste0(a, m[3], b))
    }
    s
  }, character(1), USE.NAMES = FALSE)
}

#' Human-readable trait label
#'
#' @param trait_col Trait code.
#'
#' @return A display label.
#'
#' @noRd
.hbr_trait_label <- function(trait_col) {
  map <- c(
    NT  = "Tiller number (NT)",  TN  = "Tiller number (TN)",
    PH  = "Plant height (PH)",   GY  = "Grain yield (GY)",
    GW  = "Grain weight (GW)",   TGW = "1000-grain weight (TGW)",
    DTF = "Days to flowering (DTF)", DTH = "Days to heading (DTH)",
    HD  = "Heading date (HD)",   SL  = "Seedling length (SL)",
    RL  = "Root length (RL)")
  if (trait_col %in% names(map)) unname(map[[trait_col]])
  else trait_col
}


# ----------------------------------------------------------------------------
# 6. Haplotype summary table builders
# ----------------------------------------------------------------------------

#' Build the haplotype summary table
#'
#' @param pieces Gene pieces from [`.hbr_load_gene_pieces()`].
#' @param highlight_superior Mark the superior haplotype with a symbol.
#' @param include_nt Include NT sequence column.
#' @param include_aa Include AA sequence column.
#' @param include_per_snp Include one column per SNP.
#' @param include_mean_se Include Mean and SE columns.
#' @param superior_symbol Symbol used to mark the superior row.
#' @param palette Palette name.
#'
#' @return A data.table.
#'
#' @noRd
.hbr_build_hap_table <- function(pieces,
                                  highlight_superior = TRUE,
                                  include_nt         = FALSE,
                                  include_aa         = FALSE,
                                  include_per_snp    = TRUE,
                                  include_mean_se    = TRUE,
                                  superior_symbol    = "\u2020",
                                  palette            = "haplobreed") {

  hap_map   <- pieces$hap_map
  groups    <- pieces$groups
  func_snps <- pieces$func_snps

  sup_hap <- if ("Superior_Haplotype" %in% names(groups))
    as.character(groups$Superior_Haplotype[1])
  else as.character(groups$Haplotype[which.max(groups$Mean)])

  tbl <- data.table::data.table(
    Haplotype  = as.character(hap_map$Haplotype),
    N          = as.integer(hap_map$Count),
    `Freq (%)` = round(100 * as.numeric(hap_map$Frequency), 1))

  if (isTRUE(include_nt)) tbl[["NT"]] <- as.character(hap_map$NT_Sequence)
  if (isTRUE(include_aa)) tbl[["AA"]] <- as.character(hap_map$AA_Sequence)

  if (isTRUE(include_per_snp) && nrow(func_snps)) {
    alleles_list <- strsplit(as.character(hap_map$NT_Sequence), "-",
                             fixed = TRUE)
    for (j in seq_len(nrow(func_snps))) {
      snp <- func_snps[j]
      hdr <- .hbr_short_snp_header(snp, j)
      tbl[[hdr]] <- vapply(alleles_list,
                           function(a) if (length(a) >= j) a[j] else "",
                           character(1))
    }
  }

  tbl[["Superior"]] <- if (isTRUE(highlight_superior))
    ifelse(as.character(tbl$Haplotype) == sup_hap, superior_symbol, "")
  else ""

  if (isTRUE(include_mean_se)) {
    gd <- data.table::data.table(
      Haplotype = as.character(groups$Haplotype),
      Mean      = round(as.numeric(groups$Mean), 2),
      SE        = round(as.numeric(groups$SE), 2))
    tbl <- merge(tbl, gd, by = "Haplotype", all.x = TRUE)
  }

  ord <- groups$Haplotype[order(-as.numeric(groups$Mean))]
  tbl[match(ord, tbl$Haplotype)]
}

#' Render a haplotype table with ggplot2
#'
#' @param tbl_df Table data.
#' @param fill_mat Cell colour matrix.
#' @param fontface_mat Font face matrix.
#' @param title Optional title.
#' @param palette Palette name.
#' @param theme Theme name.
#'
#' @return A ggplot object.
#'
#' @noRd
.hbr_build_hap_table_gg <- function(tbl_df, fill_mat, fontface_mat = NULL,
                                    title   = NULL,
                                    palette = "haplobreed",
                                    theme   = "haplobreed") {

  srf <- .hbr_theme_surface(theme)

  n_row <- nrow(tbl_df)
  n_col <- ncol(tbl_df)
  col_names <- names(tbl_df)

  cell_df <- data.table::data.table(
    row   = rep(seq_len(n_row), times = n_col),
    col   = rep(seq_len(n_col), each  = n_row),
    value = as.vector(as.matrix(tbl_df)),
    fill  = as.vector(fill_mat),
    face  = if (is.null(fontface_mat)) "plain"
            else as.vector(fontface_mat))

  header_df <- data.table::data.table(
    row   = 0L,
    col   = seq_len(n_col),
    value = col_names,
    fill  = srf$table_header)

  cell_plain <- cell_df[cell_df$face == "plain"]
  cell_bold  <- cell_df[cell_df$face == "bold"]

  rules <- data.table::data.table(
    y   = c(0.5, -0.5, -n_row - 0.5),
    lwd = c(1.0, 0.7, 1.0))

  ggplot2::ggplot() +
    ggplot2::geom_tile(data = header_df,
                       ggplot2::aes(x = col, y = -row, fill = fill),
                       colour = NA, width = 1, height = 1) +
    ggplot2::geom_tile(data = cell_df,
                       ggplot2::aes(x = col, y = -row, fill = fill),
                       colour = NA, width = 1, height = 1) +
    ggplot2::scale_fill_identity() +
    ggplot2::geom_text(data = header_df,
                       ggplot2::aes(x = col, y = -row, label = value),
                       fontface = "bold", size = 2.9,
                       colour = srf$text_primary) +
    ggplot2::geom_text(data = cell_plain,
                       ggplot2::aes(x = col, y = -row, label = value),
                       size = 2.7, colour = srf$text_primary) +
    ggplot2::geom_text(data = cell_bold,
                       ggplot2::aes(x = col, y = -row, label = value),
                       fontface = "bold", size = 2.7,
                       colour = srf$text_primary) +
    ggplot2::geom_segment(data = rules,
                          ggplot2::aes(x = 0.5, xend = n_col + 0.5,
                                       y = y, yend = y,
                                       linewidth = lwd),
                          colour = srf$rule_strong, lineend = "butt") +
    ggplot2::scale_linewidth_identity() +
    ggplot2::scale_x_continuous(expand = c(0, 0), breaks = NULL) +
    ggplot2::scale_y_continuous(expand = c(0, 0), breaks = NULL) +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(
        face = "bold", size = 11, hjust = 0,
        colour = srf$text_primary,
        margin = ggplot2::margin(b = 4)),
      plot.margin     = ggplot2::margin(6, 10, 6, 10),
      plot.background = ggplot2::element_rect(fill = "transparent",
                                              colour = NA))
}


# ----------------------------------------------------------------------------
# 7. gtable haplotype table renderer
# ----------------------------------------------------------------------------

#' Render a haplotype table as a grob
#'
#' @param gene_id Gene ID.
#' @param tbl Haplotype table.
#' @param highlight_superior Highlight the superior row.
#' @param groups Groups table.
#' @param palette Palette name.
#' @param theme Theme name.
#'
#' @return A gtable object.
#'
#' @noRd
.hbr_render_hap_table_grob <- function(gene_id, tbl, highlight_superior,
                                        groups,
                                        palette = "haplobreed",
                                        theme   = "haplobreed") {

  sup_hap <- if ("Superior_Haplotype" %in% names(groups))
    as.character(groups$Superior_Haplotype[1])
  else as.character(groups$Haplotype[which.max(groups$Mean)])

  col_names <- names(tbl)
  meta_cols <- c("Haplotype", "N", "Freq (%)", "Superior")
  snp_cols  <- setdiff(col_names, c(meta_cols, "NT_Seq", "AA_Seq"))

  srf   <- .hbr_theme_surface(theme)
  cells <- hbr_cell_colours(theme, palette)

  fill_mat <- .hbr_hap_table_fill(
    tbl, snp_cols,
    superior_haplotype = if (highlight_superior) sup_hap else NULL,
    cell_colours       = cells,
    neutral            = srf$table_bg,
    superior_row       = srf$table_superior_row)

  fnt_core <- if      (ncol(tbl) <=  5) 9.0
  else if (ncol(tbl) <=  8) 8.5
  else if (ncol(tbl) <= 12) 8.0
  else if (ncol(tbl) <= 16) 7.0
  else if (ncol(tbl) <= 20) 6.2
  else 5.5

  tg <- gridExtra::tableGrob(
    tbl, rows = NULL,
    theme = gridExtra::ttheme_minimal(
      core = list(
        fg_params = list(fontsize = fnt_core, col = srf$text_primary),
        bg_params = list(fill = fill_mat, col = NA)),
      colhead = list(
        fg_params = list(fontface = "bold", fontsize = fnt_core + 0.5,
                         col = srf$text_primary),
        bg_params = list(fill = srf$table_header, col = NA))))

  title_grob <- grid::textGrob(
    paste0(gene_id, " \u2014 haplotype structure"),
    gp = grid::gpar(fontface = "bold", fontsize = 12,
                    col = srf$text_primary), vjust = 0)
  tg <- gtable::gtable_add_rows(tg, grid::unit(1.6, "lines"), 0)
  tg <- gtable::gtable_add_grob(tg, title_grob, 1, 1, 1, ncol(tg))
  tg
}


# ----------------------------------------------------------------------------
# 8. LD heatmap + gene track
# ----------------------------------------------------------------------------

#' Pairwise LD matrix
#'
#' @param geno Genotype matrix.
#' @param stat `"r2"`.
#'
#' @return Upper-triangular r-squared matrix.
#'
#' @noRd
.hbr_compute_ld_matrix <- function(geno, stat = "r2") {
  n <- ncol(geno)
  out <- matrix(NA_real_, n, n,
                dimnames = list(colnames(geno), colnames(geno)))
  for (i in seq_len(n)) {
    out[i, i] <- 1
    if (i == n) next
    for (j in (i + 1):n) {
      xi <- geno[, i]; xj <- geno[, j]
      ok <- !is.na(xi) & !is.na(xj)
      if (sum(ok) < 5) next
      r <- suppressWarnings(stats::cor(xi[ok], xj[ok]))
      out[i, j] <- if (is.na(r)) NA_real_ else r^2
    }
  }
  out
}

#' Assemble the LD heatmap figure
#'
#' @param gene_id Gene ID.
#' @param ld_mat LD matrix.
#' @param func_snps Functional SNP table.
#' @param stat `"r2"`.
#' @param palette Palette name.
#' @param theme Theme name.
#' @param reference Optional reference for the gene model.
#'
#' @return A grob.
#'
#' @noRd
.hbr_build_ld_figure <- function(gene_id, ld_mat, func_snps, stat,
                                  palette   = "haplobreed",
                                  theme     = "haplobreed",
                                  reference = NULL) {

  th       <- hbr_theme(theme, base_size = 10)
  srf      <- .hbr_theme_surface(theme)
  accent   <- hbr_semantic("superior",      palette, theme)
  accent_d <- hbr_semantic("superior_dark", palette, theme)
  ld_pal   <- hbr_semantic("ld_gradient",   palette, theme)

  n <- nrow(ld_mat)
  if (n < 2) return(NULL)

  df <- data.frame(
    i = rep(seq_len(n), times = n),
    j = rep(seq_len(n), each  = n),
    value = as.vector(ld_mat))
  df <- df[df$i <= df$j, ]
  df$i <- factor(df$i, levels = seq_len(n))
  df$j <- factor(df$j, levels = seq_len(n))

  labels <- if ("Codon_Num" %in% names(func_snps) &&
                !all(is.na(func_snps$Codon_Num)))
    ifelse(is.na(func_snps$Codon_Num), paste0("p", func_snps$POS),
           paste0("c.", func_snps$Codon_Num))
  else paste0("p", func_snps$POS)

  p_ld <- ggplot2::ggplot(df, ggplot2::aes(x = j, y = i, fill = value)) +
    ggplot2::geom_tile(colour = srf$panel_bg, linewidth = 0.3) +
    ggplot2::scale_fill_gradientn(
      colours = ld_pal,
      limits  = c(0, 1), na.value = srf$grid_major,
      name    = if (stat == "r2") expression(r^2) else "D'") +
    ggplot2::scale_x_discrete(labels = labels, expand = c(0, 0)) +
    ggplot2::scale_y_discrete(labels = rev(labels), expand = c(0, 0)) +
    ggplot2::coord_fixed() +
    ggplot2::labs(title = paste0(gene_id, " \u2014 LD heatmap"),
                  x = NULL, y = NULL) +
    th +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.ticks = ggplot2::element_blank())

  p_gene <- NULL
  if (!is.null(reference)) {
    gene_gff_id <- .hbr_resolve_gene_gff_id(gene_id, func_snps, reference)
    if (!is.null(gene_gff_id)) {
      p_gene <- tryCatch(
        .hbr_gs_build_one(
          gene_id          = gene_id,
          gene_gff_id      = gene_gff_id,
          reference        = reference,
          variants         = func_snps,
          show_variants    = TRUE,
          variant_shape_by_region = TRUE,
          aa_full_names    = TRUE,
          palette          = palette,
          theme            = theme),
        error = function(e) NULL)

      if (!is.null(p_gene)) {
        p_gene <- p_gene +
          ggplot2::labs(title = NULL, subtitle = NULL) +
          ggplot2::theme(
            plot.title    = ggplot2::element_blank(),
            plot.subtitle = ggplot2::element_blank(),
            axis.text.x   = ggplot2::element_blank(),
            axis.ticks.x  = ggplot2::element_blank(),
            axis.title.x  = ggplot2::element_blank(),
            plot.margin   = ggplot2::margin(2, 4, 2, 4))
      }
    }
  }

  if (is.null(p_gene)) {
    pos <- as.numeric(func_snps$POS); pos <- pos[!is.na(pos)]
    if (!length(pos)) return(p_ld)
    x0 <- min(pos) - 200; x1 <- max(pos) + 200
    p_gene <- ggplot2::ggplot() +
      ggplot2::annotate("segment", x = x0, xend = x1, y = 0, yend = 0,
                        colour = srf$text_secondary, linewidth = 1) +
      ggplot2::geom_point(data = data.frame(x = pos, y = 0),
                          ggplot2::aes(x = x, y = y),
                          inherit.aes = FALSE, shape = 21, size = 2.5,
                          fill = accent, colour = accent_d, stroke = 0.5) +
      ggplot2::coord_cartesian(xlim = c(x0, x1), ylim = c(-0.5, 0.5)) +
      ggplot2::theme_void()
  }

  gridExtra::arrangeGrob(
    p_gene, p_ld, ncol = 1,
    heights = if (is.null(reference)) c(0.15, 0.85) else c(0.32, 0.68))
}


# ----------------------------------------------------------------------------
# 9. Haplotype network builder
# ----------------------------------------------------------------------------

#' Build a haplotype network figure
#'
#' @param gene_id Gene ID.
#' @param hap_map Haplotype map.
#' @param groups Groups table.
#' @param colour_by_superior Colour the superior node.
#' @param label_superior Label the superior node.
#' @param layout `"circle"` or `"grid"`.
#' @param node_scale Node radius multiplier.
#' @param palette Palette name.
#' @param theme Theme name.
#'
#' @return A ggplot object or `NULL`.
#'
#' @noRd
.hbr_build_network_figure <- function(gene_id, hap_map, groups,
                                       colour_by_superior = TRUE,
                                       label_superior     = TRUE,
                                       layout             = "circle",
                                       node_scale         = 1.0,
                                       palette            = "haplobreed",
                                       theme              = "haplobreed") {

  n_hap <- nrow(hap_map)
  if (n_hap < 2) return(NULL)

  srf <- .hbr_theme_surface(theme)
  pal <- .hbr_palette_anchors(palette)

  accent   <- pal$superior
  accent_d <- pal$superior_dark
  other    <- pal$other
  other_d  <- pal$other_dark

  seqs   <- as.character(hap_map$NT_Sequence)
  parts  <- strsplit(seqs, "-", fixed = TRUE)
  n_snps <- max(lengths(parts))
  mat <- t(vapply(parts,
                  function(x) c(x, rep("N", n_snps - length(x))),
                  character(n_snps)))

  edges <- matrix(integer(0), ncol = 2, nrow = 0)
  for (i in seq_len(n_hap - 1)) {
    for (j in (i + 1):n_hap) {
      if (sum(mat[i, ] != mat[j, ]) == 1)
        edges <- rbind(edges, c(i, j))
    }
  }

  if (layout == "circle") {
    angles <- seq(0, 2 * pi, length.out = n_hap + 1)[seq_len(n_hap)]
    coords <- cbind(x = cos(angles), y = sin(angles))
  } else {
    ncol <- ceiling(sqrt(n_hap)); nrow <- ceiling(n_hap / ncol)
    grid <- expand.grid(x = seq_len(ncol), y = seq_len(nrow))
    coords <- as.matrix(grid[seq_len(n_hap), ])
  }

  freq   <- as.numeric(hap_map$Frequency)
  max_f  <- max(freq, na.rm = TRUE)
  if (!is.finite(max_f) || max_f <= 0) max_f <- 1
  radius <- 3 * node_scale * (0.5 + 0.5 * freq / max_f)

  sup_hap <- if ("Superior_Haplotype" %in% names(groups))
    as.character(groups$Superior_Haplotype[1])
  else as.character(groups$Haplotype[which.max(groups$Mean)])
  is_sup <- as.character(hap_map$Haplotype) == sup_hap

  fill_col <- if (isTRUE(colour_by_superior))
    ifelse(is_sup, accent, other) else rep(other, n_hap)
  edge_col <- ifelse(is_sup, accent_d, other_d)

  nodes <- data.frame(x = coords[, 1], y = coords[, 2],
                      radius = radius,
                      label  = as.character(hap_map$Haplotype),
                      sup    = is_sup, stringsAsFactors = FALSE)
  nodes$label <- ifelse(nodes$sup & isTRUE(label_superior),
                        paste0(nodes$label, "*"), nodes$label)

  p <- ggplot2::ggplot()
  if (nrow(edges)) {
    edf <- data.frame(x = nodes$x[edges[, 1]], y = nodes$y[edges[, 1]],
                      xend = nodes$x[edges[, 2]],
                      yend = nodes$y[edges[, 2]])
    p <- p + ggplot2::geom_segment(
      data = edf,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      colour = srf$rule_light, linewidth = 0.6)
  }

  p +
    ggplot2::geom_point(data = nodes,
                        ggplot2::aes(x = x, y = y, size = radius),
                        shape = 21, fill = fill_col,
                        colour = edge_col, stroke = 1.2) +
    ggplot2::scale_size_identity() +
    ggplot2::geom_text(data = nodes,
                       ggplot2::aes(x = x, y = y, label = label),
                       nudge_y = 0.15, size = 3, fontface = "bold",
                       colour = srf$text_primary) +
    ggplot2::coord_equal() +
    ggplot2::theme_void(base_size = 11) +
    ggplot2::labs(title = paste0(gene_id, " \u2014 haplotype network"),
                  subtitle = paste0(n_hap, " haplotypes \u2022 superior = ",
                                    sup_hap)) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", hjust = 0.5,
                                            colour = srf$text_primary),
      plot.subtitle = ggplot2::element_text(hjust = 0.5,
                                            colour = srf$text_secondary))
}


# ----------------------------------------------------------------------------
# 10. GWAS normalisation
# ----------------------------------------------------------------------------

#' Normalise a GWAS table
#'
#' @param gwas_data A data.frame with GWAS results.
#' @param p_col,chr_col,pos_col,snp_col Column names. `NULL` auto-detects.
#'
#' @return A data.table with `SNP`, `Chr`, `POS`, `P`, or `NULL`.
#'
#' @noRd
.hbr_normalize_gwas_data <- function(gwas_data,
                                      p_col   = NULL,
                                      chr_col = NULL,
                                      pos_col = NULL,
                                      snp_col = NULL) {
  if (!is.data.frame(gwas_data)) return(NULL)
  gd <- data.table::as.data.table(gwas_data)
  p_col   <- p_col   %||% intersect(c("P", "p", "P_value", "Pvalue",
                                      "PVAL", "pval"), names(gd))[1]
  chr_col <- chr_col %||% intersect(c("Chromosome", "CHR", "Chr", "chr"),
                                    names(gd))[1]
  pos_col <- pos_col %||% intersect(c("Position", "POS", "Pos", "pos",
                                      "BP", "bp"), names(gd))[1]
  snp_col <- snp_col %||% intersect(c("SNP", "Snp", "snp",
                                      "Marker", "ID", "id"), names(gd))[1]
  if (is.na(p_col) || is.na(chr_col) || is.na(pos_col)) return(NULL)
  snp_ids <- if (!is.na(snp_col)) as.character(gd[[snp_col]])
  else paste0(gd[[chr_col]], ":", gd[[pos_col]])
  out <- data.table::data.table(
    SNP = snp_ids,
    Chr = as.character(gd[[chr_col]]),
    POS = suppressWarnings(as.numeric(gd[[pos_col]])),
    P   = suppressWarnings(as.numeric(gd[[p_col]])))
  out <- out[!is.na(POS) & !is.na(P) & P > 0 & P <= 1 & !is.na(Chr)]
  if (!nrow(out)) return(NULL)
  out[]
}


# ----------------------------------------------------------------------------
# 11. Reference cache and gene structure
# ----------------------------------------------------------------------------

#' Session reference cache
#'
#' Environment holding the GFF3-keyed cache used by
#' [`.hbr_resolve_reference()`].
#'
#' @noRd
.hbr_reference_cache <- new.env(parent = emptyenv())

#' Cache key for a GFF3
#'
#' @param gff3_path Path to a GFF3 file.
#'
#' @return A string, or `NA_character_` if the file does not exist.
#'
#' @noRd
.hbr_cache_key <- function(gff3_path) {
  if (is.null(gff3_path) || !file.exists(gff3_path)) return(NA_character_)
  info <- file.info(gff3_path)
  paste0(normalizePath(gff3_path, winslash = "/", mustWork = FALSE),
         "|", as.numeric(info$mtime), "|", info$size)
}

#' Look up a cached reference
#'
#' @param gff3_path Path to a GFF3 file.
#'
#' @return A reference object or `NULL`.
#'
#' @noRd
.hbr_get_cached_reference <- function(gff3_path) {
  key <- .hbr_cache_key(gff3_path)
  if (is.na(key)) return(NULL)
  if (exists(key, envir = .hbr_reference_cache, inherits = FALSE))
    return(get(key, envir = .hbr_reference_cache, inherits = FALSE))
  NULL
}

#' Store a reference in the cache
#'
#' @param gff3_path Path to a GFF3 file.
#' @param reference Reference object.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.hbr_set_cached_reference <- function(gff3_path, reference) {
  key <- .hbr_cache_key(gff3_path)
  if (is.na(key)) return(invisible(NULL))
  assign(key, reference, envir = .hbr_reference_cache)
  invisible(NULL)
}

#' Clear the session GFF3 cache
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
hbr_clear_reference_cache <- function() {
  rm(list = ls(envir = .hbr_reference_cache, all.names = TRUE),
     envir = .hbr_reference_cache)
  invisible(NULL)
}

#' Resolve an hbr_reference
#'
#' Accepts a reference, a result object, or a GFF3 path.
#'
#' @param reference Optional reference.
#' @param result Optional pipeline result.
#' @param gff3_path Optional GFF3 path.
#' @param verbose Print progress.
#'
#' @return A reference object or `NULL`.
#'
#' @noRd
.hbr_resolve_reference <- function(reference = NULL, result = NULL,
                                    gff3_path = NULL, verbose = FALSE) {

  if (!is.null(reference) && inherits(reference, "hbr_reference"))
    return(reference)

  if (!is.null(result) && is.list(result) &&
      !is.null(result$reference) &&
      inherits(result$reference, "hbr_reference"))
    return(result$reference)

  if (!is.null(gff3_path) && file.exists(gff3_path)) {
    cached <- .hbr_get_cached_reference(gff3_path)
    if (!is.null(cached)) return(cached)
    if (verbose)
      message("  building hbr_reference from GFF3: ", basename(gff3_path),
              " (cached for this session)")
    ref <- tryCatch(
      hbr_load_reference(genome = basename(gff3_path),
                         gff3   = gff3_path,
                         annotation = "CUSTOM"),
      error = function(e) {
        if (verbose)
          message("  reference build failed: ", conditionMessage(e))
        NULL
      })
    if (!is.null(ref)) .hbr_set_cached_reference(gff3_path, ref)
    return(ref)
  }
  NULL
}

#' Resolve a GFF3 gene_id
#'
#' @param gene_id Gene ID.
#' @param func_snps Functional SNP table.
#' @param reference Reference object.
#'
#' @return A GFF3 gene_id or `NULL`.
#'
#' @noRd
.hbr_resolve_gene_gff_id <- function(gene_id, func_snps, reference) {
  if (is.null(reference)) return(NULL)

  if ("Transcript_ID" %in% names(func_snps)) {
    tx_ids <- unique(as.character(func_snps$Transcript_ID))
    tx_ids <- tx_ids[!is.na(tx_ids) & nzchar(tx_ids)]
    if (length(tx_ids)) {
      tx_tbl <- data.table::as.data.table(reference$transcripts)
      tx_gene <- as.character(tx_tbl$gene_id)
      tx_id   <- as.character(tx_tbl$transcript_id)
      m <- match(tx_ids, tx_id)
      m <- m[!is.na(m)]
      if (length(m)) return(tx_gene[m[1]])
    }
  }

  gene_ids <- as.character(reference$genes$gene_id)
  m <- match(gene_id, gene_ids)
  if (!is.na(m)) return(gene_ids[m])

  m <- match(tolower(gene_id), tolower(gene_ids))
  if (!is.na(m)) return(gene_ids[m])

  NULL
}

#' Format a genomic span
#'
#' Auto Mb / kb / bp with thousands separators.
#'
#' @param x_start,x_end Start and end.
#'
#' @return A single string.
#'
#' @noRd
.hbr_format_span <- function(x_start, x_end) {
  span <- x_end - x_start
  if (span >= 1e6) {
    sprintf("%.2f\u2013%.2f Mb (%.2f Mb)",
            x_start / 1e6, x_end / 1e6, span / 1e6)
  } else if (span >= 1e3) {
    sprintf("%s\u2013%s kb (%.1f kb)",
            format(round(x_start / 1e3), big.mark = ","),
            format(round(x_end / 1e3), big.mark = ","),
            span / 1e3)
  } else {
    sprintf("%s\u2013%s bp",
            format(round(x_start), big.mark = ","),
            format(round(x_end), big.mark = ","))
  }
}

#' Gene-model colour list
#'
#' @param palette Palette name.
#' @param theme Theme name.
#'
#' @return A named list of colours.
#'
#' @noRd
.hbr_gene_model_colours <- function(palette = "haplobreed",
                                     theme   = "scientific") {

  srf     <- .hbr_theme_surface(theme)
  is_dark <- srf$is_dark
  pal     <- hbr_palette(palette)

  adjust <- function(hex, amount) {
    rgb <- grDevices::col2rgb(hex) / 255
    if (amount < 0) {
      f <- 1 + amount
      grDevices::rgb(rgb[1] * f, rgb[2] * f, rgb[3] * f)
    } else {
      grDevices::rgb(rgb[1] + (1 - rgb[1]) * amount,
                     rgb[2] + (1 - rgb[2]) * amount,
                     rgb[3] + (1 - rgb[3]) * amount)
    }
  }

  if (identical(theme, "bw")) {
    return(list(
      utr_fill      = "#E0E0E0",
      utr_border    = "#9E9E9E",
      exon_fill     = "#424242",
      exon_border   = "#000000",
      intron        = "#757575",
      marker_fill   = "#000000",
      marker_border = "#000000",
      star          = "#000000",
      panel_bg      = "#F5F5F5",
      panel_border  = "#BDBDBD",
      text          = "#000000"))
  }

  exon_base   <- pal[1]
  marker_base <- if (length(pal) >= 2) pal[2] else adjust(pal[1], 0.5)

  if (is_dark) {
    list(
      utr_fill      = adjust(exon_base, -0.45),
      utr_border    = adjust(exon_base, -0.25),
      exon_fill     = adjust(exon_base, -0.15),
      exon_border   = adjust(exon_base,  0.15),
      intron        = srf$text_muted,
      marker_fill   = marker_base,
      marker_border = adjust(marker_base, -0.25),
      star          = marker_base,
      panel_bg      = srf$table_bg,
      panel_border  = srf$rule_light,
      text          = srf$text_primary)
  } else {
    list(
      utr_fill      = adjust(exon_base,  0.60),
      utr_border    = adjust(exon_base,  0.25),
      exon_fill     = adjust(exon_base, -0.30),
      exon_border   = adjust(exon_base, -0.55),
      intron        = srf$text_secondary,
      marker_fill   = marker_base,
      marker_border = adjust(marker_base, -0.35),
      star          = marker_base,
      panel_bg      = srf$table_bg,
      panel_border  = srf$rule_light,
      text          = srf$text_primary)
  }
}


# ----------------------------------------------------------------------------
# 12. Gene-structure plot
# ----------------------------------------------------------------------------

#' Resolve the variants argument
#'
#' @param genes Gene IDs.
#' @param fun_variants A list, a data.frame, or `NULL`.
#' @param variant_panel_dir Directory with per-gene variant files.
#'
#' @return A named list of variant tables.
#'
#' @noRd
.hbr_gs_resolve_variants <- function(genes, fun_variants, variant_panel_dir) {
  out <- stats::setNames(vector("list", length(genes)), genes)

  if (is.list(fun_variants) && !is.data.frame(fun_variants) &&
      !is.null(names(fun_variants))) {
    for (g in genes)
      if (!is.null(fun_variants[[g]]))
        out[[g]] <- data.table::as.data.table(fun_variants[[g]])
    return(out)
  }

  if (is.data.frame(fun_variants)) {
    dt <- data.table::as.data.table(fun_variants)
    gc <- intersect(c("Gene_ID", "RAPDB_Gene_ID", "MSU_Gene_ID", "Gene"),
                    names(dt))[1]
    if (!is.na(gc)) {
      for (g in genes) {
        sub <- dt[as.character(get(gc)) == g]
        if (nrow(sub)) out[[g]] <- sub
      }
    }
    return(out)
  }

  if (!is.null(variant_panel_dir) && dir.exists(variant_panel_dir)) {
    for (g in genes) {
      f <- file.path(variant_panel_dir, paste0(g, "_functional_snps.csv"))
      if (file.exists(f)) {
        out[[g]] <- tryCatch(
          data.table::fread(f),
          error = function(e) NULL)
      }
    }
    return(out)
  }

  out
}

#' Draw a single gene-structure panel
#'
#' @param gene_id Gene ID.
#' @param gene_gff_id GFF3 gene_id.
#' @param reference Reference object.
#' @param variants Variant table.
#' @param show_variants Draw variant markers.
#' @param variant_shape_by_region Shape by region.
#' @param aa_full_names Use three-letter amino-acid codes.
#' @param palette Palette name.
#' @param theme Theme name.
#'
#' @return A ggplot object.
#'
#' @noRd
.hbr_gs_build_one <- function(gene_id, gene_gff_id, reference, variants,
                               show_variants = TRUE,
                               variant_shape_by_region = TRUE,
                               aa_full_names = TRUE,
                               palette = "haplobreed",
                               theme   = "scientific") {

  srf <- .hbr_theme_surface(theme)
  pal <- .hbr_palette_anchors(palette)
  th  <- hbr_theme(theme, base_size = 11)
  gm  <- .hbr_gene_model_colours(palette, theme)

  accent   <- pal$superior
  accent_d <- pal$superior_dark

  genes <- data.table::as.data.table(reference$genes)
  row <- genes[genes$gene_id == gene_gff_id, , drop = FALSE]
  if (!nrow(row))
    row <- genes[tolower(genes$gene_id) == tolower(gene_gff_id), ,
                 drop = FALSE]
  if (!nrow(row)) stop("Gene not found in reference: ", gene_id)

  g_start <- as.numeric(row$Start[1])
  g_end   <- as.numeric(row$End[1])
  strand  <- as.character(row$Strand[1])
  chr     <- as.character(row$Chr[1])
  if (!is.finite(g_start) || !is.finite(g_end) || g_end <= g_start)
    stop("Invalid gene coordinates for ", gene_id)

  tx_tbl <- data.table::as.data.table(reference$transcripts)
  tx_ids <- tx_tbl$transcript_id[tx_tbl$gene_id == gene_gff_id]
  if (!length(tx_ids)) stop("No transcripts for ", gene_id)

  cds <- if (!is.null(reference$cds))
    data.table::as.data.table(reference$cds)[
      data.table::as.data.table(reference$cds)$transcript_id %in% tx_ids, ,
      drop = FALSE] else NULL
  utr <- if (!is.null(reference$utr))
    data.table::as.data.table(reference$utr)[
      data.table::as.data.table(reference$utr)$transcript_id %in% tx_ids, ,
      drop = FALSE] else NULL

  collapse <- function(dt) {
    if (is.null(dt) || !nrow(dt)) return(NULL)
    d <- unique(dt[, c("start", "end")])
    data.table::setnames(d, c("start", "end"))
    d[, start := as.numeric(start)]
    d[, end   := as.numeric(end)]
    d <- d[!is.na(start) & !is.na(end)]
    if (!nrow(d)) return(NULL)
    data.table::setorder(d, start, end)
    d[, grp := cumsum(c(TRUE, start[-1] > cummax(end[-.N])))]
    d[, .(start = min(start), end = max(end)), by = grp][, grp := NULL][]
  }
  cds_u <- collapse(cds)
  utr_u <- collapse(utr)

  span     <- g_end - g_start
  intron_y <- 0.50
  exon_y   <- c(0.30, 0.70)
  utr_y    <- c(0.38, 0.62)

  p <- ggplot2::ggplot() +
    ggplot2::annotate("rect",
                      xmin = g_start - span * 0.03,
                      xmax = g_end   + span * 0.03,
                      ymin = 0.10, ymax = 0.90,
                      fill = gm$panel_bg, colour = gm$panel_border,
                      linewidth = 0.3) +
    ggplot2::annotate("segment",
                      x = g_start, xend = g_end,
                      y = intron_y, yend = intron_y,
                      colour = gm$intron, linewidth = 0.5)

  if (!is.null(utr_u) && nrow(utr_u))
    p <- p + ggplot2::geom_rect(
      data = utr_u,
      ggplot2::aes(xmin = start, xmax = end,
                   ymin = utr_y[1], ymax = utr_y[2]),
      inherit.aes = FALSE,
      fill = gm$utr_fill, colour = gm$utr_border, linewidth = 0.18)

  if (!is.null(cds_u) && nrow(cds_u))
    p <- p + ggplot2::geom_rect(
      data = cds_u,
      ggplot2::aes(xmin = start, xmax = end,
                   ymin = exon_y[1], ymax = exon_y[2]),
      inherit.aes = FALSE,
      fill = gm$exon_fill, colour = gm$exon_border, linewidth = 0.25)

  n_chev <- max(2L, min(8L, as.integer(round(span / 600))))
  chev_x <- seq(g_start + span * 0.10, g_end - span * 0.10,
                length.out = n_chev)
  chev_dx <- span * 0.015 * if (identical(strand, "-")) -1 else 1
  p <- p + ggplot2::geom_segment(
    data = data.frame(x = chev_x),
    ggplot2::aes(x = x, xend = x + chev_dx,
                 y = intron_y, yend = intron_y),
    inherit.aes = FALSE,
    arrow = ggplot2::arrow(length = grid::unit(0.10, "cm"),
                           type = "open"),
    colour = gm$intron, linewidth = 0.5)

  v <- NULL
  if (isTRUE(show_variants) && !is.null(variants) &&
      nrow(variants) && "POS" %in% names(variants)) {

    v <- data.table::as.data.table(variants)
    v <- v[!is.na(POS) & POS >= g_start & POS <= g_end]

    if (nrow(v)) {
      if (isTRUE(variant_shape_by_region) && "Region" %in% names(v)) {
        v[, vshape := data.table::fcase(
          Region == "CDS_nonsynonymous", 21L,
          Region == "CDS_synonymous",    1L,
          Region %in% c("5UTR", "3UTR"), 22L,
          default = 21L)]
      } else {
        v[, vshape := 21L]
      }

      aa_lab <- if ("AA_Change" %in% names(v)) as.character(v$AA_Change)
      else rep("", nrow(v))
      if (isTRUE(aa_full_names)) aa_lab <- .hbr_expand_aa_labels(aa_lab)
      v[, vlabel := aa_lab]

      p <- p +
        ggplot2::geom_segment(
          data = v,
          ggplot2::aes(x = POS, xend = POS, y = exon_y[2], yend = 1.00),
          inherit.aes = FALSE,
          colour = gm$marker_border, linewidth = 0.55) +
        ggplot2::geom_point(
          data = v,
          ggplot2::aes(x = POS, y = 1.08, shape = I(vshape)),
          inherit.aes = FALSE,
          fill = gm$marker_fill, colour = gm$marker_border,
          size = 3.6, stroke = 0.6)

      lbl <- v[!is.na(vlabel) & nzchar(vlabel)]
      if (nrow(lbl)) {
        if (requireNamespace("ggrepel", quietly = TRUE)) {
          p <- p + ggrepel::geom_text_repel(
            data = lbl,
            ggplot2::aes(x = POS, y = 1.08, label = vlabel),
            inherit.aes = FALSE,
            size = 3.0, fontface = "bold",
            colour = gm$text,
            direction    = "y",
            nudge_y      = 0.22,
            segment.size = 0.25,
            segment.colour = srf$text_muted,
            max.overlaps = Inf,
            seed         = 42,
            box.padding  = 0.4,
            point.padding = 0.25)
        } else {
          p <- p + ggplot2::geom_text(
            data = lbl,
            ggplot2::aes(x = POS, y = 1.35, label = vlabel),
            inherit.aes = FALSE,
            size = 3.0, fontface = "bold",
            colour = gm$text, hjust = 0.5, vjust = 0)
        }
      }
    }
  }

  five_end  <- if (identical(strand, "-")) g_end   else g_start
  three_end <- if (identical(strand, "-")) g_start else g_end
  pad <- span * 0.03
  x0  <- g_start - pad
  x1  <- g_end   + pad

  p <- p +
    ggplot2::annotate("text", x = five_end, y = 0.95,
                      label = "5'", size = 3.0, fontface = "bold",
                      colour = gm$text, vjust = 0) +
    ggplot2::annotate("text", x = three_end, y = 0.95,
                      label = "3'", size = 3.0, fontface = "bold",
                      colour = gm$text, vjust = 0)

  leg_y <- 0.02
  leg_x <- x0
  leg_step <- span * 0.16

  p <- p +
    ggplot2::annotate("rect",
                      xmin = leg_x, xmax = leg_x + span * 0.035,
                      ymin = leg_y - 0.05, ymax = leg_y + 0.05,
                      fill = gm$utr_fill, colour = gm$utr_border,
                      linewidth = 0.3) +
    ggplot2::annotate("text", x = leg_x + span * 0.045, y = leg_y,
                      label = "UTR", hjust = 0, size = 2.7,
                      colour = gm$text) +
    ggplot2::annotate("rect",
                      xmin = leg_x + leg_step,
                      xmax = leg_x + leg_step + span * 0.035,
                      ymin = leg_y - 0.05, ymax = leg_y + 0.05,
                      fill = gm$exon_fill, colour = gm$exon_border,
                      linewidth = 0.3) +
    ggplot2::annotate("text",
                      x = leg_x + leg_step + span * 0.045, y = leg_y,
                      label = "Exon", hjust = 0, size = 2.7,
                      colour = gm$text) +
    ggplot2::annotate("segment",
                      x = leg_x + 2 * leg_step,
                      xend = leg_x + 2 * leg_step + span * 0.035,
                      y = leg_y, yend = leg_y,
                      colour = gm$intron, linewidth = 0.5) +
    ggplot2::annotate("text",
                      x = leg_x + 2 * leg_step + span * 0.045, y = leg_y,
                      label = "Intron", hjust = 0, size = 2.7,
                      colour = gm$text) +
    ggplot2::annotate("point",
                      x = leg_x + 3 * leg_step, y = leg_y,
                      shape = 21, size = 3.0,
                      fill = gm$marker_fill, colour = gm$marker_border,
                      stroke = 0.5) +
    ggplot2::annotate("text",
                      x = leg_x + 3 * leg_step + span * 0.045, y = leg_y,
                      label = "Variant", hjust = 0, size = 2.7,
                      colour = gm$text)

  span_label <- .hbr_format_span(g_start, g_end)
  ymax <- if (!is.null(v) && nrow(v)) 1.60 else 1.10
  chr_prefix <- if (grepl("^chr", chr, ignore.case = TRUE)) "" else "Chr "

  p +
    ggplot2::scale_x_continuous(
      labels = function(val) {
        if (span >= 1e6)      sprintf("%.3f Mb", val / 1e6)
        else if (span >= 1e3) sprintf("%s kb",
                                      format(round(val / 1e3),
                                             big.mark = ","))
        else                  format(round(val), big.mark = ",")
      }) +
    ggplot2::coord_cartesian(xlim = c(x0, x1), ylim = c(-0.12, ymax),
                             clip = "off") +
    ggplot2::labs(
      x = paste0("Position on chromosome ", chr, " (bp)"),
      y = NULL,
      title    = gene_id,
      subtitle = paste0(chr_prefix, chr, "  \u2022  ", span_label,
                        "  \u2022  ", strand, " strand",
                        if (!is.null(row$gene_function) &&
                            !is.na(row$gene_function[1]))
                          paste0("\n", row$gene_function[1]) else "")) +
    th +
    ggplot2::theme(
      axis.text.y   = ggplot2::element_blank(),
      axis.ticks.y  = ggplot2::element_blank(),
      panel.grid    = ggplot2::element_blank(),
      plot.title    = ggplot2::element_text(face = "bold", size = 13,
                                            colour = gm$text),
      plot.subtitle = ggplot2::element_text(size = 9.5,
                                            colour = srf$text_secondary),
      plot.margin   = ggplot2::margin(6, 8, 4, 8))
}

#' Build a combined-axis gene-structure facet figure
#'
#' @param genes Gene IDs.
#' @param reference Reference object.
#' @param variant_source Named list of variant tables.
#' @param palette Palette name.
#' @param theme Theme name.
#' @param aa_full_names Use three-letter amino-acid codes.
#'
#' @return A patchwork object or `NULL`.
#'
#' @noRd
.hbr_gs_build_facet <- function(genes, reference, variant_source,
                                 palette, theme, aa_full_names) {
  if (!requireNamespace("patchwork", quietly = TRUE)) return(NULL)
  panels <- lapply(genes, function(g) {
    gene_gff_id <- .hbr_resolve_gene_gff_id(
      g, data.table::data.table(Transcript_ID = character(0)), reference)
    if (is.null(gene_gff_id)) gene_gff_id <- g
    tryCatch(
      .hbr_gs_build_one(
        gene_id = g, gene_gff_id = gene_gff_id, reference = reference,
        variants = variant_source[[g]], show_variants = TRUE,
        variant_shape_by_region = TRUE, aa_full_names = aa_full_names,
        palette = palette, theme = theme),
      error = function(e) NULL)
  })
  panels <- panels[!vapply(panels, is.null, logical(1))]
  if (!length(panels)) return(NULL)
  do.call(patchwork::wrap_plots, c(panels, list(ncol = 1)))
}

#' Add a bold panel title
#'
#' @param p A ggplot object.
#' @param title Title.
#' @param colour Title colour.
#'
#' @return The modified ggplot.
#'
#' @noRd
.hbr_panel_title <- function(p, title, colour = "#222222") {
  p + ggplot2::labs(title = title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face     = "bold",
        size     = 11,
        hjust    = 0,
        margin   = ggplot2::margin(b = 4),
        colour   = colour))
}


# ----------------------------------------------------------------------------
# 13. Trait histogram
# ----------------------------------------------------------------------------

#' Trait histogram panel
#'
#' @param values Numeric vector.
#' @param trait_label Trait label.
#' @param bins Number of histogram bins.
#' @param palette Palette name.
#' @param theme Theme name.
#'
#' @return A ggplot object.
#'
#' @noRd
.hbr_build_trait_histogram <- function(values,
                                        trait_label = "Trait",
                                        bins        = 30L,
                                        palette     = "haplobreed",
                                        theme       = "scientific") {

  srf <- .hbr_theme_surface(theme)
  pal <- hbr_palette(palette)
  th  <- hbr_theme(theme, base_size = 10)

  df <- data.frame(x = suppressWarnings(as.numeric(values)))
  df <- df[is.finite(df$x), , drop = FALSE]

  if (!nrow(df)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::annotate("text", x = 0, y = 0,
                               label = "No trait data",
                               size = 3, colour = srf$text_muted))
  }

  m  <- mean(df$x, na.rm = TRUE)
  sd <- stats::sd(df$x, na.rm = TRUE)

  rng <- range(df$x, na.rm = TRUE)
  if (!is.finite(diff(rng)) || diff(rng) <= 0)
    rng <- c(rng[1] - 0.5, rng[2] + 0.5)

  bar_cols <- grDevices::colorRampPalette(pal)(max(2L,
                                                    as.integer(bins)))

  ggplot2::ggplot(df, ggplot2::aes(x = x)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density),
                   fill = ggplot2::after_stat(x)),
      bins       = as.integer(bins),
      colour     = srf$panel_bg,
      linewidth  = 0.25,
      alpha      = 0.9) +
    ggplot2::scale_fill_gradientn(colours = bar_cols, guide = "none") +
    ggplot2::geom_density(colour = srf$text_primary, linewidth = 0.9) +
    ggplot2::annotate(
      "text",
      x = min(rng), y = Inf,
      label = sprintf("Mean: %.2f\nSd: %.2f\nN = %d", m, sd, nrow(df)),
      hjust = 0, vjust = 1.15,
      size = 3.0, fontface = "bold",
      colour = srf$text_primary) +
    ggplot2::labs(x = trait_label, y = "Density", title = NULL) +
    th +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      axis.text        = ggplot2::element_text(size = 9),
      axis.title       = ggplot2::element_text(face = "bold", size = 10),
      plot.margin      = ggplot2::margin(4, 6, 4, 4))
}


# ----------------------------------------------------------------------------
# 14. Four-panel assembler
# ----------------------------------------------------------------------------

#' Assemble the four-panel figure
#'
#' @param gene_id Gene ID.
#' @param pieces Gene pieces.
#' @param panels Panels to include.
#' @param trait_col Trait column name.
#' @param trait_label Trait label.
#' @param palette Palette name.
#' @param theme Theme name.
#' @param table_style `"color_coded"` or `"row_highlight"`.
#' @param include_nt,include_aa Include sequence columns.
#' @param aa_full_names Expand amino-acid codes.
#' @param reference Reference object.
#'
#' @return A composed figure.
#'
#' @noRd
.hbr_build_4panel <- function(gene_id, pieces, panels, trait_col,
                               trait_label   = NULL,
                               palette       = "haplobreed",
                               theme         = "haplobreed",
                               table_style   = c("color_coded",
                                                 "row_highlight"),
                               include_nt    = FALSE,
                               include_aa    = FALSE,
                               aa_full_names = TRUE,
                               reference     = NULL) {

  table_style <- match.arg(table_style)
  .hbr_check_palette_theme(palette, theme)
  th <- hbr_theme(theme, base_size = 11)

  pal_anchors <- .hbr_palette_anchors(palette)
  srf         <- .hbr_theme_surface(theme)
  cells       <- hbr_cell_colours(theme, palette)

  accent      <- pal_anchors$superior
  accent_dark <- pal_anchors$superior_dark
  txt_primary <- srf$text_primary
  txt_second  <- srf$text_secondary
  txt_muted   <- srf$text_muted

  uniform_margin <- ggplot2::margin(6, 10, 6, 10)

  if (is.null(trait_label) || !nzchar(trait_label))
    trait_label <- .hbr_trait_label(trait_col)

  hap_map   <- pieces$hap_map
  groups    <- pieces$groups
  func_snps <- pieces$func_snps
  ind_df    <- pieces$individuals_df
  meta      <- pieces$meta %||% list()

  gene_chr <- meta$chromosome
  if (is.null(gene_chr) || is.na(gene_chr)) {
    if ("Chr" %in% names(func_snps)) {
      chrs <- unique(as.character(func_snps$Chr))
      chrs <- chrs[!is.na(chrs) & nzchar(chrs)]
      if (length(chrs)) gene_chr <- chrs[1]
    }
  }
  chr_label <- if (!is.null(gene_chr) && !is.na(gene_chr))
    as.character(gene_chr) else "\u2014"

  superior_hap <- if ("Superior_Haplotype" %in% names(groups))
    as.character(groups$Superior_Haplotype[1])
  else if (!is.null(meta$superior_hap) && !is.na(meta$superior_hap))
    meta$superior_hap
  else as.character(groups$Haplotype[which.max(groups$Mean)])

  p_gene_track <- if ("gene_model" %in% panels) {

    p_struct <- NULL
    if (!is.null(reference)) {
      gene_gff_id <- .hbr_resolve_gene_gff_id(gene_id, func_snps, reference)
      if (!is.null(gene_gff_id)) {
        p_struct <- tryCatch(
          .hbr_gs_build_one(
            gene_id          = gene_id,
            gene_gff_id      = gene_gff_id,
            reference        = reference,
            variants         = func_snps,
            show_variants    = TRUE,
            variant_shape_by_region = TRUE,
            aa_full_names    = aa_full_names,
            palette          = palette,
            theme            = theme),
          error = function(e) NULL)
      }
    }

    if (!is.null(p_struct)) {
      p_struct <- .hbr_panel_title(p_struct, "A.i   Functional variants",
                                    colour = txt_primary) +
        ggplot2::labs(subtitle = NULL) +
        ggplot2::theme(plot.subtitle = ggplot2::element_blank())
      p_struct
    } else {
      pos <- as.numeric(func_snps$POS); pos <- pos[!is.na(pos)]
      if (length(pos)) {
        g_start <- min(pos); g_end <- max(pos)
        x0 <- g_start - 200; x1 <- g_end + 200

        g <- ggplot2::ggplot() +
          ggplot2::annotate("segment", x = x0, xend = x1, y = 0, yend = 0,
                            colour = txt_primary, linewidth = 0.7) +
          ggplot2::geom_segment(
            data = data.frame(x = pos, y = 0, xend = pos, yend = 0.55),
            ggplot2::aes(x = x, xend = xend, y = y, yend = yend),
            inherit.aes = FALSE, colour = accent, linewidth = 0.55) +
          ggplot2::geom_point(data = data.frame(x = pos, y = 0.55),
                              ggplot2::aes(x = x, y = y),
                              inherit.aes = FALSE, shape = 21, size = 3.0,
                              fill = accent, colour = accent_dark,
                              stroke = 0.55)

        ns <- func_snps[
          func_snps$Region == "CDS_nonsynonymous" &
            !is.na(func_snps$AA_Change) &
            nzchar(func_snps$AA_Change), ]
        if (nrow(ns)) {
          aa_lab <- as.character(ns$AA_Change)
          if (isTRUE(aa_full_names)) aa_lab <- .hbr_expand_aa_labels(aa_lab)
          lbl_df <- data.frame(x = as.numeric(ns$POS),
                               y = 0.72, label = aa_lab)
          if (requireNamespace("ggrepel", quietly = TRUE)) {
            g <- g + ggrepel::geom_text_repel(
              data = lbl_df,
              ggplot2::aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 2.9, fontface = "bold",
              colour = txt_primary, angle = 45,
              hjust = 0, vjust = 0, direction = "y", nudge_y = 0.1,
              segment.size = 0.2, segment.colour = txt_muted,
              max.overlaps = Inf, seed = 42)
          } else {
            g <- g + ggplot2::geom_text(
              data = lbl_df,
              ggplot2::aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 2.9, fontface = "bold",
              colour = txt_primary, angle = 45, hjust = 0, vjust = 0)
          }
        }

        span    <- x1 - x0
        x_scale <- if (span >= 1e6)
          ggplot2::scale_x_continuous(
            labels = function(v) sprintf("%.2f Mb", v / 1e6))
        else if (span >= 1e3)
          ggplot2::scale_x_continuous(
            labels = function(v) sprintf("%s kb",
                                         format(round(v / 1e3),
                                                big.mark = ",")))
        else
          ggplot2::scale_x_continuous(
            labels = function(v) format(round(v), big.mark = ","))

        g2 <- g + x_scale +
          ggplot2::coord_cartesian(xlim = c(x0, x1),
                                   ylim = c(-0.25, 1.15)) +
          ggplot2::labs(
            x = paste0("Position on chromosome ", chr_label, " (bp)"),
            y = NULL,
            subtitle = paste0("Chr ", chr_label, "  \u2022  ",
                              .hbr_format_span(g_start, g_end))) +
          th +
          ggplot2::theme(
            axis.text.y   = ggplot2::element_blank(),
            axis.ticks.y  = ggplot2::element_blank(),
            panel.grid    = ggplot2::element_blank(),
            plot.margin   = uniform_margin,
            plot.subtitle = ggplot2::element_text(
              size = 8.5, colour = txt_second, hjust = 0))
        .hbr_panel_title(g2, "A.i   Functional variants",
                         colour = txt_primary)
      } else {
        .hbr_panel_title(
          ggplot2::ggplot() + ggplot2::theme_void() +
            ggplot2::annotate("text", x = 0, y = 0,
                              label = "No functional variants",
                              size = 4, colour = txt_muted),
          "A.i   Functional variants", colour = txt_primary)
      }
    }
  } else NULL

  p_network <- if ("haplotype_network" %in% panels &&
                   !is.null(hap_map) && nrow(hap_map) >= 2) {

    net_p <- tryCatch(
      .hbr_build_network_figure(
        gene_id            = gene_id,
        hap_map            = hap_map,
        groups             = groups,
        colour_by_superior = TRUE,
        label_superior     = TRUE,
        layout             = "circle",
        node_scale         = 1.0,
        palette            = palette,
        theme              = theme),
      error = function(e) NULL)

    if (!is.null(net_p)) {
      net_p <- net_p +
        ggplot2::labs(title = "A.ii   Haplotype network",
                      subtitle = NULL) +
        ggplot2::theme(
          plot.title    = ggplot2::element_text(
            face = "bold", size = 11, hjust = 0,
            margin = ggplot2::margin(b = 4),
            colour = txt_primary),
          plot.subtitle = ggplot2::element_blank(),
          plot.margin   = ggplot2::margin(4, 6, 4, 6))
      net_p
    } else NULL
  } else NULL

  tbl_df <- .hbr_build_hap_table(
    pieces             = pieces,
    highlight_superior = TRUE,
    include_nt         = include_nt,
    include_aa         = include_aa,
    include_per_snp    = TRUE,
    include_mean_se    = TRUE,
    superior_symbol    = "\u2020",
    palette            = palette)

  n_col     <- ncol(tbl_df)
  col_names <- names(tbl_df)
  seq_cols  <- intersect(c("NT", "AA"), col_names)
  meta_cols <- c("Haplotype", "N", "Freq (%)", "Superior", "Mean", "SE")
  snp_cols  <- setdiff(col_names, c(meta_cols, seq_cols))

  fill_mat <- if (table_style == "color_coded" && length(snp_cols)) {
    .hbr_hap_table_fill(tbl_df, snp_cols,
                        superior_haplotype = superior_hap,
                        cell_colours       = cells,
                        neutral            = srf$table_bg,
                        superior_row       = srf$table_superior_row)
  } else {
    m <- matrix(srf$table_bg, nrow(tbl_df), n_col,
                dimnames = list(NULL, col_names))
    sup_i <- which(tbl_df$Haplotype == superior_hap)
    if (length(sup_i)) m[sup_i, ] <- srf$table_superior_row
    m
  }

  sup_row <- which(tbl_df$Haplotype == superior_hap)
  fg_mat  <- matrix("plain", nrow(tbl_df), n_col,
                    dimnames = list(NULL, col_names))
  if (length(sup_row)) fg_mat[sup_row, ] <- "bold"

  p_table <- .hbr_build_hap_table_gg(
    tbl_df       = tbl_df,
    fill_mat     = fill_mat,
    fontface_mat = fg_mat,
    title        = "B.   Functional haplotypes",
    palette      = palette,
    theme        = theme)

  p_bp <- if ("trait_boxplot" %in% panels &&
              !is.null(ind_df) && nrow(ind_df) &&
              trait_col %in% names(ind_df)) {

    df <- data.table::as.data.table(ind_df)
    df[, Haplotype := as.character(Haplotype)]
    df[, trait_val := as.numeric(get(trait_col))]
    df <- df[!is.na(trait_val)]

    ord <- groups$Haplotype[order(-as.numeric(groups$Mean))]
    ord <- ord[!is.na(ord)]
    df <- df[Haplotype %in% ord]
    df[, Haplotype := factor(Haplotype, levels = ord)]

    base_cols <- rep(srf$text_muted, length(ord))
    names(base_cols) <- ord
    if (superior_hap %in% names(base_cols))
      base_cols[superior_hap] <- accent

    y_max <- max(df$trait_val, na.rm = TRUE)
    y_min <- min(df$trait_val, na.rm = TRUE)
    y_rng <- y_max - y_min
    if (!is.finite(y_rng) || y_rng <= 0) y_rng <- 1

    letter_col <- intersect(c("Letter", "letter",
                              "Significance_Group",
                              "Duncan_Group", "Tukey_Group"),
                            names(groups))[1]

    lbl_df <- NULL
    if (!is.na(letter_col)) {
      grp_max <- df[, .(y = max(trait_val, na.rm = TRUE) + y_rng * 0.05),
                    by = Haplotype]
      lbl_df <- merge(
        data.table::data.table(
          Haplotype = as.character(groups$Haplotype),
          Letter    = as.character(groups[[letter_col]])),
        grp_max, by = "Haplotype", all.x = TRUE)
      lbl_df[, Haplotype := factor(Haplotype, levels = ord)]
      lbl_df <- lbl_df[!is.na(y)]
    }

    n_by_hap <- df[, .N, by = Haplotype]
    x_labels <- stats::setNames(
      sprintf("%s\n(N = %d)", n_by_hap$Haplotype, n_by_hap$N),
      as.character(n_by_hap$Haplotype))

    p <- ggplot2::ggplot(df, ggplot2::aes(x = Haplotype, y = trait_val,
                                          fill = Haplotype)) +
      ggplot2::geom_boxplot(width = 0.58, linewidth = 0.55,
                            alpha = 0.85, outlier.shape = NA,
                            colour = txt_primary)

    if (requireNamespace("ggbeeswarm", quietly = TRUE)) {
      p <- p + ggbeeswarm::geom_quasirandom(
        width = 0.22, size = 1.35, alpha = 0.55,
        colour = txt_primary, shape = 16)
    } else {
      p <- p + ggplot2::geom_jitter(
        width = 0.12, size = 1.35, alpha = 0.45,
        colour = txt_primary, seed = 42)
    }

    p <- p +
      ggplot2::stat_summary(fun = mean, geom = "point", shape = 23,
                            fill = "white", colour = txt_primary,
                            size = 2.6, stroke = 0.55) +
      ggplot2::scale_fill_manual(values = base_cols) +
      ggplot2::scale_x_discrete(labels = x_labels) +
      ggplot2::coord_cartesian(
        ylim = c(y_min - y_rng * 0.05, y_max + y_rng * 0.20),
        clip = "off")

    if (!is.null(lbl_df) && nrow(lbl_df)) {
      p <- p + ggplot2::geom_text(
        data = lbl_df,
        ggplot2::aes(x = Haplotype, y = y, label = Letter),
        inherit.aes = FALSE, fontface = "bold",
        size = 4.5, colour = txt_primary, na.rm = TRUE)
    }

    sub_parts <- character(0)
    if (!is.null(meta$test) && !is.na(meta$test))
      sub_parts <- c(sub_parts, as.character(meta$test))
    if (!is.null(meta$p_value) && is.finite(meta$p_value))
      sub_parts <- c(sub_parts, paste0("P = ", hbr_format_p(meta$p_value)))
    if (!is.null(meta$permutation_p) && is.finite(meta$permutation_p))
      sub_parts <- c(sub_parts,
                     paste0("P_perm = ",
                            hbr_format_p(meta$permutation_p)))
    if (!is.null(meta$adjusted_p) && is.finite(meta$adjusted_p))
      sub_parts <- c(sub_parts,
                     paste0("FDR = ", hbr_format_p(meta$adjusted_p)))
    if (!is.null(meta$cohens_d) && is.finite(meta$cohens_d))
      sub_parts <- c(sub_parts, sprintf("d = %.2f", meta$cohens_d))
    stat_txt <- if (length(sub_parts))
      paste(sub_parts, collapse = "\n") else NULL

    if (!is.null(stat_txt)) {
      p <- p + ggplot2::annotate(
        "label", x = Inf, y = Inf, label = stat_txt,
        hjust = 1.05, vjust = 1.10, size = 2.7,
        linewidth = 0.25,
        label.padding = grid::unit(0.3, "lines"),
        fill = srf$table_bg, colour = txt_primary)
    }

    p <- p +
      ggplot2::labs(x = NULL, y = trait_label, subtitle = NULL) +
      th +
      ggplot2::theme(
        legend.position    = "none",
        panel.grid.minor   = ggplot2::element_blank(),
        panel.grid.major.x = ggplot2::element_blank(),
        plot.margin        = uniform_margin,
        axis.text.x        = ggplot2::element_text(size = 10,
                                                   colour = txt_primary),
        axis.title.y       = ggplot2::element_text(face = "bold",
                                                   size = 11))

    .hbr_panel_title(p, "C.i   Distribution by haplotype",
                     colour = txt_primary)
  } else if ("trait_boxplot" %in% panels) {
    .hbr_panel_title(
      ggplot2::ggplot() + ggplot2::theme_void() +
        ggplot2::annotate("text", x = 0, y = 0,
                          label = "No individual-level trait data",
                          size = 4, colour = txt_muted),
      "C.i   Distribution by haplotype", colour = txt_primary)
  } else NULL

  p_hist <- if ("trait_boxplot" %in% panels &&
                !is.null(ind_df) && nrow(ind_df) &&
                trait_col %in% names(ind_df)) {

    hist_p <- .hbr_build_trait_histogram(
      values      = ind_df[[trait_col]],
      trait_label = trait_label,
      bins        = 30L,
      palette     = palette,
      theme       = theme)

    .hbr_panel_title(hist_p, "C.ii   Overall distribution",
                     colour = txt_primary)
  } else NULL

  has_pw <- requireNamespace("patchwork", quietly = TRUE)

  panel_A <- if (!is.null(p_gene_track) && !is.null(p_network) && has_pw) {
    patchwork::wrap_plots(p_gene_track, p_network,
                          ncol = 2, widths = c(0.6, 0.4))
  } else if (!is.null(p_gene_track) && !is.null(p_network) && !has_pw) {
    gridExtra::arrangeGrob(
      ggplot2::ggplotGrob(p_gene_track),
      ggplot2::ggplotGrob(p_network),
      ncol = 2, widths = c(0.6, 0.4))
  } else if (!is.null(p_gene_track)) {
    p_gene_track
  } else {
    p_network
  }

  panel_C <- if (!is.null(p_bp) && !is.null(p_hist) && has_pw) {
    patchwork::wrap_plots(p_bp, p_hist,
                          ncol = 2, widths = c(0.6, 0.4))
  } else if (!is.null(p_bp) && !is.null(p_hist) && !has_pw) {
    gridExtra::arrangeGrob(
      ggplot2::ggplotGrob(p_bp),
      ggplot2::ggplotGrob(p_hist),
      ncol = 2, widths = c(0.6, 0.4))
  } else if (!is.null(p_bp)) {
    p_bp
  } else {
    p_hist
  }

  h_A <- if (!is.null(p_network)) 1.10 else 0.90
  h_B <- 0.42 + nrow(tbl_df) * 0.030
  h_C <- max(0.9, 0.55 + nrow(tbl_df) * 0.05)

  all_panels  <- list(panel_A, p_table, panel_C)
  all_heights <- c(h_A, h_B, h_C)

  keep        <- !vapply(all_panels, is.null, logical(1))
  all_panels  <- all_panels[keep]
  all_heights <- all_heights[keep]

  if (!length(all_panels)) {
    return(grid::textGrob(
      paste0(gene_id, " \u2014 no panels to draw"),
      gp = grid::gpar(fontface = "italic", fontsize = 11,
                      col = txt_primary)))
  }

  master_title <- paste0(gene_id, " \u2014 functional haplotype summary")

  if (has_pw) {
    pw <- do.call(
      patchwork::wrap_plots,
      c(all_panels, list(ncol = 1L, heights = all_heights)))

    pw <- pw + patchwork::plot_annotation(
      title = master_title,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold", size = 13, hjust = 0,
          colour = txt_primary,
          margin = ggplot2::margin(b = 6))))

    return(pw)
  }

  gt_list <- lapply(all_panels, function(g) {
    if (inherits(g, "gtable")) return(g)
    ggplot2::ggplotGrob(g)
  })

  do.call(gridExtra::arrangeGrob,
          c(gt_list,
            list(ncol    = 1,
                 heights = all_heights,
                 top     = grid::textGrob(
                   master_title,
                   gp = grid::gpar(fontface = "bold", fontsize = 13,
                                   col = txt_primary),
                   x = grid::unit(10, "pt"), hjust = 0, vjust = 1),
                 padding = grid::unit(2, "lines"))))
}


# ----------------------------------------------------------------------------
# 15. Progressive-zoom tier
# ----------------------------------------------------------------------------

#' Progressive-zoom tier
#'
#' @param gwas_norm Normalised GWAS table.
#' @param xlim X-axis limits.
#' @param chr Chromosome label.
#' @param title Tier title.
#' @param highlight_range Highlight range.
#' @param box_range Box range.
#' @param y_max Y-axis maximum.
#' @param threshold Horizontal threshold line.
#' @param palette Palette name.
#' @param theme Theme name.
#' @param ld_to_lead LD to lead SNP.
#' @param lead Lead SNP info.
#' @param stats_box Optional stats box text.
#'
#' @return A ggplot object.
#'
#' @noRd
.hbr_build_zoom_tier <- function(gwas_norm, xlim,
                                  chr = NULL, title = NULL,
                                  highlight_range = NULL,
                                  box_range       = NULL,
                                  y_max           = NULL,
                                  threshold       = NULL,
                                  palette         = "haplobreed",
                                  theme           = "haplobreed",
                                  ld_to_lead      = NULL,
                                  lead            = NULL,
                                  stats_box       = NULL) {

  th     <- hbr_theme(theme, base_size = 10)
  srf    <- .hbr_theme_surface(theme)
  gwas_l <- hbr_semantic("gwas_lead",   palette, theme)
  gwas_w <- hbr_semantic("gwas_window", palette, theme)

  gd <- data.table::as.data.table(gwas_norm)[POS >= xlim[1] & POS <= xlim[2]]
  if (!nrow(gd)) {
    return(ggplot2::ggplot() + ggplot2::theme_void() +
             ggplot2::annotate("text", x = 0.5, y = 0.5,
                               label = "No SNPs in range",
                               size = 4, colour = srf$text_muted))
  }

  gd <- data.table::copy(gd)
  gd[, neglogP := -log10(pmax(P, 1e-300))]
  if (is.null(y_max)) y_max <- max(gd$neglogP, na.rm = TRUE) * 1.15
  if (!is.finite(y_max) || y_max <= 0) y_max <- 1

  use_ld <- !is.null(ld_to_lead) && any(names(ld_to_lead) %in% gd$SNP)

  if (use_ld) {
    gd[, r2 := as.numeric(ld_to_lead[SNP])]
    gd[is.na(r2), r2 := 0]
    ld_pal <- c("#3B7DD8", "#4CAF50", "#FFC107", "#FF9800", "#D62728")
    p <- ggplot2::ggplot(gd, ggplot2::aes(x = POS, y = neglogP))
    if (!is.null(highlight_range)) {
      p <- p + ggplot2::annotate("rect",
                                 xmin = highlight_range[1],
                                 xmax = highlight_range[2],
                                 ymin = -Inf, ymax = Inf,
                                 fill = gwas_w, alpha = 0.35)
    }
    p <- p + ggplot2::geom_point(
      ggplot2::aes(colour = r2), size = 1.6, alpha = 0.9) +
      ggplot2::scale_colour_gradientn(
        colours = ld_pal, limits = c(0, 1),
        name    = expression(r^2),
        guide   = ggplot2::guide_colourbar(
          barwidth  = grid::unit(0.5, "cm"),
          barheight = grid::unit(2.0, "cm")))
  } else {
    gd[, in_highlight := if (!is.null(highlight_range))
      POS >= highlight_range[1] & POS <= highlight_range[2] else FALSE]
    p <- ggplot2::ggplot(gd, ggplot2::aes(x = POS, y = neglogP))
    if (!is.null(highlight_range)) {
      p <- p + ggplot2::annotate("rect",
                                 xmin = highlight_range[1],
                                 xmax = highlight_range[2],
                                 ymin = -Inf, ymax = Inf,
                                 fill = gwas_w, alpha = 0.55)
    }
    p <- p + ggplot2::geom_point(
      ggplot2::aes(colour = in_highlight), size = 1.4, alpha = 0.85) +
      ggplot2::scale_colour_manual(
        values = c(`TRUE` = gwas_l, `FALSE` = srf$text_muted),
        guide = "none")
  }

  if (!is.null(threshold) && is.finite(threshold)) {
    p <- p + ggplot2::geom_hline(yintercept = threshold,
                                 linetype   = "dashed",
                                 colour     = srf$text_muted,
                                 linewidth  = 0.4)
  }

  if (!is.null(box_range)) {
    p <- p + ggplot2::annotate("rect",
                               xmin = box_range[1], xmax = box_range[2],
                               ymin = -Inf, ymax = Inf,
                               fill = NA, colour = gwas_l,
                               linetype = "dashed", linewidth = 0.6)
  }

  if (!is.null(lead) && is.list(lead) && !is.null(lead$pos) &&
      is.finite(lead$pos) &&
      lead$pos >= xlim[1] && lead$pos <= xlim[2]) {

    y_lead <- lead$neglogP %||% y_max * 0.95

    p <- p +
      ggplot2::annotate("point", x = lead$pos, y = y_lead,
                        shape = 21, size = 7,
                        fill = NA, colour = gwas_l, stroke = 1.2) +
      ggplot2::annotate("point", x = lead$pos, y = y_lead,
                        shape = 8, size = 4.5,
                        colour = gwas_l, stroke = 1.4) +
      ggplot2::annotate("text", x = lead$pos, y = y_lead,
                        label = lead$label %||% "",
                        vjust = -1.6, hjust = 0.5, size = 3.0,
                        fontface = "bold", colour = gwas_l)
  }

  span <- diff(xlim)
  x_scale <- if (span > 1e6) {
    ggplot2::scale_x_continuous(
      labels = function(v) sprintf("%.1f Mb", v / 1e6))
  } else if (span > 1e3) {
    ggplot2::scale_x_continuous(
      labels = function(v) sprintf("%.0f kb", v / 1e3))
  } else {
    ggplot2::scale_x_continuous(
      labels = function(v) sprintf("%d bp", round(v)))
  }

  sub_parts <- c(
    if (!is.null(chr))   paste0("Chr ", chr),
    if (!is.null(title)) title,
    paste0("window: ",
           if (span > 1e6) sprintf("%.2f Mb", span / 1e6)
           else if (span > 1e3) sprintf("%.0f kb", span / 1e3)
           else sprintf("%d bp", round(span))))

  p <- p + x_scale +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, y_max * 1.10)) +
    ggplot2::labs(x = NULL, y = expression(-log[10](P)),
                  subtitle = paste(sub_parts, collapse = "  \u2022  ")) +
    th +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(colour = srf$grid_major,
                                               linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.y     = ggplot2::element_text(face = "bold", size = 10,
                                               colour = srf$text_primary),
      plot.subtitle    = ggplot2::element_text(size = 9,
                                               colour = srf$text_secondary),
      legend.position  = if (use_ld) "right" else "none",
      legend.title     = ggplot2::element_text(size = 8),
      legend.text      = ggplot2::element_text(size = 7))

  if (!is.null(stats_box) && nzchar(stats_box)) {
    p <- p + ggplot2::annotate(
      "label",
      x = Inf, y = -Inf,
      label = stats_box,
      hjust = 1.05, vjust = -0.20,
      size = 2.6, linewidth = 0.25,
      label.padding = grid::unit(0.3, "lines"),
      fill = srf$table_bg, colour = srf$text_primary)
  }

  p
}


# ----------------------------------------------------------------------------
# 16. Shared utilities
# ----------------------------------------------------------------------------

#' Add a letter tag above a panel
#'
#' @param g A grob or ggplot.
#' @param label Tag label.
#' @param size Font size.
#' @param pad_pt Left padding.
#' @param colour Text colour.
#'
#' @return A gtable.
#'
#' @noRd
.hbr_add_panel_tag <- function(g, label, size = 11, pad_pt = 4,
                                colour = "#222222") {
  if (inherits(g, "ggplot") || inherits(g, "gg"))
    g <- ggplot2::ggplotGrob(g)
  if (!inherits(g, "gtable"))
    g <- gtable::gtable_matrix("panel",
                               matrix(list(g), 1, 1),
                               grid::unit(1, "npc"),
                               grid::unit(1, "npc"))
  g <- gtable::gtable_add_rows(g, grid::unit(2.0, "lines"), 0)
  tag <- grid::textGrob(label,
                        x = grid::unit(pad_pt, "pt"),
                        y = grid::unit(2, "pt"),
                        hjust = 0, vjust = 1,
                        gp = grid::gpar(fontface = "bold",
                                        fontsize = size, col = colour))
  gtable::gtable_add_grob(g, tag, 1, 1, 1, ncol(g), clip = "off")
}

#' Align panel widths
#'
#' @param grobs List of grobs.
#'
#' @return A list of aligned gtables.
#'
#' @noRd
.hbr_align_panels <- function(grobs) {
  gt <- lapply(grobs, function(g) {
    if (inherits(g, "gtable")) return(g)
    if (inherits(g, "ggplot") || inherits(g, "gg"))
      return(ggplot2::ggplotGrob(g))
    gtable::gtable_matrix("wrap", matrix(list(g), 1, 1),
                          grid::unit(1, "npc"), grid::unit(1, "npc"))
  })
  find_panel_col <- function(g) {
    if (!"layout" %in% names(g)) return(NA_integer_)
    hit <- which(g$layout$name == "panel")
    if (!length(hit)) return(NA_integer_)
    g$layout$l[hit[1]]
  }
  panel_cols <- vapply(gt, find_panel_col, numeric(1))
  have       <- !is.na(panel_cols)
  if (sum(have) >= 2) {
    widths <- lapply(gt[have], function(g) g$widths[find_panel_col(g)])
    max_w  <- do.call(grid::unit.pmax, widths)
    for (i in which(have)) {
      pc <- panel_cols[i]
      gt[[i]]$widths[pc] <- max_w
    }
  }
  gt
}

#' Resolve the figures directory
#'
#' @param output_dir Output directory.
#' @param result Optional result object.
#' @param figures_dir Subdirectory name.
#'
#' @return A directory path.
#'
#' @noRd
.hbr_resolve_figures_dir <- function(output_dir, result = NULL,
                                      figures_dir = "figures") {
  if (is.null(output_dir) && !is.null(result)) output_dir <- getwd()
  if (is.null(output_dir))
    stop("Supply `output_dir` or a `result` object.", call. = FALSE)
  d <- if (!is.null(figures_dir) && nzchar(figures_dir))
    file.path(output_dir, figures_dir) else output_dir
  if (!dir.exists(d))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Save a figure in multiple formats
#'
#' @param fig The grob or ggplot.
#' @param gene_id Gene ID.
#' @param suffix File suffix.
#' @param fig_dir Output directory.
#' @param width,height Dimensions in inches.
#' @param formats Output formats.
#'
#' @return A named list of file paths.
#'
#' @noRd
.hbr_save_fig <- function(fig, gene_id, suffix, fig_dir,
                           width, height, formats) {

  if (is.null(fig_dir) || length(fig_dir) != 1L ||
      is.na(fig_dir) || !nzchar(fig_dir))
    fig_dir <- getwd()

  if (!dir.exists(fig_dir))
    dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  files <- character(0)
  for (fmt in formats) {
    fn <- file.path(fig_dir, paste0(gene_id, "_", suffix, ".", fmt))
    ok <- tryCatch({
      .hbr_draw_and_save(fig, fn, fmt, width, height)
      TRUE
    }, error = function(e) {
      warning("Could not write ", fmt, " for ", gene_id, ": ",
              conditionMessage(e), call. = FALSE)
      FALSE
    })
    if (ok) files[[fmt]] <- fn
  }
  files
}

#' Draw and save one figure file
#'
#' @param fig The grob or ggplot.
#' @param fn Output path.
#' @param fmt Format.
#' @param width,height Dimensions in inches.
#'
#' @return `NULL`, invisibly.
#'
#' @noRd
.hbr_draw_and_save <- function(fig, fn, fmt, width, height) {

  if (!is.character(fn) || length(fn) != 1L || is.na(fn) || !nzchar(fn))
    stop("Invalid output path: ", deparse(fn))
  if (!is.numeric(width)  || length(width)  != 1L || !is.finite(width))
    stop("Invalid width: ",  deparse(width))
  if (!is.numeric(height) || length(height) != 1L || !is.finite(height))
    stop("Invalid height: ", deparse(height))

  d <- dirname(fn)
  if (!dir.exists(d))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

  if (fmt == "pdf") {
    grDevices::pdf(fn, width = width, height = height)
  } else if (fmt == "png") {
    grDevices::png(fn, width = width, height = height,
                   units = "in", res = 300)
  } else if (fmt == "tiff") {
    grDevices::tiff(fn, width = width, height = height,
                    units = "in", res = 600, compression = "lzw")
  } else if (fmt == "svg") {
    grDevices::svg(fn, width = width, height = height)
  } else {
    stop("Unsupported format: ", fmt)
  }
  on.exit(grDevices::dev.off(), add = TRUE)

  if (inherits(fig, "gtable") || inherits(fig, "grob") ||
      inherits(fig, "gTree")) {
    grid::grid.newpage()
    grid::grid.draw(fig)
  } else {
    print(fig)
  }
  invisible(NULL)
}

#' Multi-page combined PDF
#'
#' @param plots Named list of figures.
#' @param fig_dir Output directory.
#' @param combined_pdf Write the combined PDF.
#' @param width,height Page dimensions in inches.
#' @param verbose Print progress.
#'
#' @return The output path or `NULL`.
#'
#' @noRd
.hbr_write_combined <- function(plots, fig_dir, combined_pdf,
                                 width, height, verbose = TRUE) {
  if (!isTRUE(combined_pdf) || length(plots) < 2) return(NULL)
  out <- file.path(fig_dir, "combined_gene_haplotypes.pdf")
  tryCatch({
    grDevices::pdf(out, width = width, height = height, onefile = TRUE)
    for (p in plots) {
      grid::grid.newpage()
      if (inherits(p, "gtable") || inherits(p, "grob") ||
          inherits(p, "gTree"))
        grid::grid.draw(p)
      else
        print(p)
    }
    grDevices::dev.off()
    if (verbose) message("\nCombined PDF written: ", out)
    out
  }, error = function(e) {
    if (grDevices::dev.cur() > 1) grDevices::dev.off()
    warning("Combined PDF failed: ", conditionMessage(e), call. = FALSE)
    NULL
  })
}

#' Short SNP column header
#'
#' @param snp SNP row.
#' @param index SNP index.
#'
#' @return A short header string.
#'
#' @noRd
.hbr_short_snp_header <- function(snp, index) {
  region <- if ("Region" %in% names(snp)) as.character(snp$Region) else NA
  codon  <- if ("Codon_Num" %in% names(snp)) snp$Codon_Num else NA
  if (is.na(region)) return(paste0("S", index))
  switch(region,
         CDS_nonsynonymous = if (!is.na(codon)) paste0("c.", codon)
         else paste0("cS", index),
         CDS_synonymous    = if (!is.na(codon)) paste0("s.", codon)
         else paste0("sS", index),
         `5UTR`            = paste0("5U.", index),
         `3UTR`            = paste0("3U.", index),
         paste0("S", index))
}

#' Zoom-gap trapezoid connector
#'
#' @param colour Line colour.
#' @param label Optional label.
#'
#' @return A grob.
#'
#' @noRd
.hbr_zoom_gap_grob <- function(colour = "#B2182B", label = NULL) {
  parts <- list(
    grid::polygonGrob(
      x = grid::unit(c(0.40, 0.60, 0.55, 0.45), "npc"),
      y = grid::unit(c(1, 1, 0, 0), "npc"),
      gp = grid::gpar(col = colour, lty = "dashed", lwd = 0.7, fill = NA)))
  if (!is.null(label) && nzchar(label)) {
    parts <- c(parts, list(
      grid::textGrob(label,
                     x = grid::unit(0.50, "npc"),
                     y = grid::unit(0.50, "npc"),
                     gp = grid::gpar(col = colour, fontsize = 8,
                                     fontface = "italic"))))
  }
  do.call(grid::grobTree, parts)
}

#' Haplotype-table cell colour matrix
#'
#' @param tbl_df Table data.
#' @param snp_cols SNP column names.
#' @param superior_haplotype Superior haplotype.
#' @param cell_colours Named vector of A/C/G/T colours.
#' @param neutral Neutral cell colour.
#' @param superior_row Superior row highlight.
#' @param ... Ignored.
#'
#' @return A character matrix of fill colours.
#'
#' @noRd
.hbr_hap_table_fill <- function(tbl_df, snp_cols,
                                 superior_haplotype = NULL,
                                 cell_colours       = NULL,
                                 neutral            = "#FAFAFA",
                                 superior_row       = "#FFF3E0",
                                 ...) {
  if (is.null(cell_colours))
    cell_colours <- c(A = "#D7EDD8", C = "#CFE2F3",
                      G = "#FCE5CD", T = "#F4CCCC")
  n_row <- nrow(tbl_df)
  n_col <- ncol(tbl_df)
  fill  <- matrix(neutral, n_row, n_col,
                  dimnames = list(NULL, names(tbl_df)))
  if (!is.null(superior_haplotype) && "Haplotype" %in% names(tbl_df)) {
    sup_i <- which(as.character(tbl_df$Haplotype) == superior_haplotype)
    if (length(sup_i)) fill[sup_i, ] <- superior_row
  }
  for (sc in snp_cols) {
    if (!sc %in% names(tbl_df)) next
    j <- which(names(tbl_df) == sc)
    vals <- toupper(as.character(tbl_df[[sc]]))
    for (i in seq_len(n_row)) {
      if (vals[i] %in% names(cell_colours))
        fill[i, j] <- cell_colours[[vals[i]]]
    }
  }
  fill
}

#' Letter and subtitle panel title
#'
#' @param p A ggplot object.
#' @param letter Letter tag.
#' @param subtitle Optional subtitle.
#' @param size Font size.
#' @param colour Text colour.
#'
#' @return The modified ggplot.
#'
#' @noRd
.hbr_letter_panel <- function(p, letter, subtitle = NULL,
                               size = 11, colour = "#222222") {
  title <- if (is.null(subtitle) || !nzchar(subtitle)) letter
  else paste0(letter, "   ", subtitle)
  p + ggplot2::labs(title = title) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold", size = size, hjust = 0,
        margin = ggplot2::margin(b = 4),
        colour = colour))
}

#' Bin a Manhattan table
#'
#' @param gd Normalised GWAS table.
#' @param bin_size Bin size in bp.
#'
#' @return A binned data.table.
#'
#' @noRd
.hbr_manhattan_binned <- function(gd, bin_size = 50000L) {
  gd <- data.table::copy(gd)
  gd[, bin := floor(GenomicPos / bin_size) * bin_size]
  gd[, .(
    GenomicPos = mean(GenomicPos),
    neglogP    = max(-log10(P), na.rm = TRUE),
    n_snps     = .N,
    top_snp    = SNP[which.max(-log10(P))]
  ), by = .(Chr, bin)][]
}

#' SNP caption strip
#'
#' @param func_snps Functional SNP table.
#' @param palette Palette name.
#' @param theme Theme name.
#'
#' @return A ggplot object.
#'
#' @noRd
.hbr_snp_caption_strip <- function(func_snps,
                                    palette = "haplobreed",
                                    theme   = "scientific") {

  srf <- .hbr_theme_surface(theme)

  if (is.null(func_snps) || !nrow(func_snps)) {
    return(ggplot2::ggplot() + ggplot2::theme_void())
  }

  n <- nrow(func_snps)

  items <- vapply(seq_len(n), function(j) {
    s <- func_snps[j]
    id <- paste0("S", j)

    change <- if ("AA_Change" %in% names(s) &&
                  !is.na(s$AA_Change) &&
                  nzchar(as.character(s$AA_Change))) {
      as.character(s$AA_Change)
    } else if ("Variant_Label" %in% names(s) &&
               !is.na(s$Variant_Label) &&
               nzchar(as.character(s$Variant_Label))) {
      as.character(s$Variant_Label)
    } else {
      ref <- if ("REF" %in% names(s)) as.character(s$REF) else "?"
      alt <- if ("ALT" %in% names(s)) sub(",.*$", "", as.character(s$ALT))
             else "?"
      paste0(ref, ">", alt)
    }

    pos <- if ("POS" %in% names(s))
      sprintf("%.2f", as.numeric(s$POS) / 1000) else ""
    pos_txt <- if (nzchar(pos)) paste0(" (", pos, " kb)") else ""

    paste0(id, " = ", change, pos_txt)
  }, character(1))

  per_row <- if (n <= 4) n else 4L
  rows <- split(items, ceiling(seq_along(items) / per_row))
  caption <- paste(
    vapply(rows, function(r) paste(r, collapse = "     "), character(1)),
    collapse = "\n")

  ggplot2::ggplot() +
    ggplot2::annotate("text",
                      x = 0.5, y = 0.5,
                      label = caption,
                      hjust = 0.5, vjust = 0.5,
                      size = 2.6,
                      colour = srf$text_secondary,
                      lineheight = 1.5) +
    ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1),
                             clip = "off") +
    ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(4, 10, 6, 10))
}

#' LD to lead SNP
#'
#' @param geno Genotype matrix.
#' @param lead_col Lead SNP column name.
#'
#' @return A named numeric vector of r-squared values.
#'
#' @noRd
.hbr_compute_ld_to_lead <- function(geno, lead_col) {
  if (is.null(geno) || is.null(lead_col) || !lead_col %in% colnames(geno))
    return(NULL)
  lead <- geno[, lead_col]
  if (all(is.na(lead))) return(NULL)
  r2 <- vapply(seq_len(ncol(geno)), function(j) {
    x <- geno[, j]
    ok <- !is.na(x) & !is.na(lead)
    if (sum(ok) < 5L) return(NA_real_)
    r <- suppressWarnings(stats::cor(x[ok], lead[ok]))
    if (is.na(r)) NA_real_ else r^2
  }, numeric(1))
  names(r2) <- colnames(geno)
  r2
}

#' Cumulative genomic positions
#'
#' @param gd Normalised GWAS table.
#'
#' @return The table with `ChrPos_Start` and `GenomicPos`.
#'
#' @noRd
.hbr_build_cumulative_positions <- function(gd) {
  chr_order <- sort(unique(gd$Chr))
  chr_num <- suppressWarnings(as.numeric(chr_order))
  if (all(!is.na(chr_num))) chr_order <- chr_order[order(chr_num)]
  offset <- 0
  gd[, ChrPos_Start := NA_real_]
  for (i in seq_along(chr_order)) {
    rows <- which(gd$Chr == chr_order[i])
    if (!length(rows)) next
    pmin <- min(gd$POS[rows], na.rm = TRUE)
    pmax <- max(gd$POS[rows], na.rm = TRUE)
    gd[rows, ChrPos_Start := pmin]
    gd[rows, GenomicPos := offset + (POS - pmin)]
    offset <- offset + (pmax - pmin) + 1
  }
  gd[]
}

#' Nearest gene label
#'
#' @param pos Position.
#' @param chr Chromosome.
#' @param reference Reference object.
#' @param max_dist Maximum distance.
#'
#' @return A list with `gene_id` and `dist`, or `NULL`.
#'
#' @noRd
.hbr_nearest_gene_label <- function(pos, chr, reference,
                                     max_dist = 200000) {
  if (is.null(reference) || !inherits(reference, "hbr_reference")) return(NULL)
  genes <- data.table::as.data.table(reference$genes)
  genes <- genes[as.character(Chr) == as.character(chr)]
  if (!nrow(genes)) return(NULL)
  genes[, dist := pmin(abs(Start - pos), abs(End - pos))]
  genes <- genes[dist <= max_dist][order(dist)]
  if (!nrow(genes)) return(NULL)
  list(gene_id = genes$gene_id[1], dist = genes$dist[1])
}

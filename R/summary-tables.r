# ============================================================================
# R/summary-tables.r
# Publication summary row construction and gene-object loading from CSV.
#
# ----------------------------------------------------------------------------
# Summary row construction
# ----------------------------------------------------------------------------

#' Map a RAP-DB gene to its MSU7 orthologue
#'
#' @param chr Gene chromosome.
#' @param start,end Gene coordinates.
#' @param all_annotation MSU7 annotation table.
#'
#' @return A one-row data.table with `MSU_Gene_ID`,
#'   `MSU_Gene_Function`, `MSU_Start`, `MSU_End`.
#'
#' @noRd
.hbr_map_to_msu <- function(chr, start, end, all_annotation) {
  empty <- data.table::data.table(
    MSU_Gene_ID       = NA_character_,
    MSU_Gene_Function = NA_character_,
    MSU_Start         = NA_integer_,
    MSU_End           = NA_integer_)
  if (is.null(all_annotation) || !nrow(all_annotation)) return(empty)
  if (!"Type" %in% names(all_annotation)) return(empty)

  genes <- all_annotation[Type == "gene" &
                           as.character(Chr) == as.character(chr) &
                           Start <= end & End >= start]
  if (!nrow(genes)) return(empty)

  genes[, overlap := pmin(End, end) - pmax(Start, start) + 1L]
  data.table::setorder(genes, -overlap)
  best <- genes[1]

  data.table::data.table(
    MSU_Gene_ID       = best$Gene_ID,
    MSU_Gene_Function = if ("Gene_Function" %in% names(best))
                          best$Gene_Function else NA_character_,
    MSU_Start         = best$Start,
    MSU_End           = best$End)
}

#' Build the pipeline's publication summary row for one gene
#'
#' @param gene_row Gene row from the candidate table.
#' @param ns_res Output of [hbr_get_functional_variants()].
#' @param gd Output of [hbr_extract_genotypes()].
#' @param hap_mat Output of [hbr_build_haplotype_matrix()].
#' @param hap_data Output of [hbr_define_haplotypes()].
#' @param test_res Output of [hbr_test_association()].
#' @param sup Output of [hbr_identify_superior()].
#' @param phenotype_type Phenotype type label.
#' @param avail_pcs Available PC columns.
#' @param all_genes Reference gene table.
#' @param all_annotation Optional full annotation table.
#' @param gwas_snp_col,gwas_pos_col Optional GWAS column names.
#'
#' @return A one-row data.table.
#'
#' @noRd
.hbr_build_summary_row <- function(gene_row, ns_res, gd, hap_mat, hap_data,
                                    test_res, sup, phenotype_type, avail_pcs,
                                    all_genes = NULL, all_annotation = NULL,
                                    gwas_snp_col = NULL, gwas_pos_col = NULL) {

  target_gene <- gene_row$RAPDB_Gene_ID

  gene_function <- NA_character_
  if (!is.null(all_genes) && "gene_function" %in% names(all_genes)) {
    gf_row <- all_genes[gene_id == target_gene]
    if (nrow(gf_row)) gene_function <- gf_row$gene_function[1]
  }

  msu <- data.table::data.table(
    MSU_Gene_ID = NA_character_, MSU_Gene_Function = NA_character_,
    MSU_Start = NA_integer_, MSU_End = NA_integer_)
  if (!is.null(all_annotation) &&
      all(c("Chromosomes", "Gene_Start", "Gene_End") %in% names(gene_row))) {
    msu <- .hbr_map_to_msu(gene_row$Chromosomes, gene_row$Gene_Start,
                            gene_row$Gene_End, all_annotation)
  }

  gwas_snp <- NA_character_; gwas_pos <- NA_real_
  if (!is.null(gwas_snp_col) && gwas_snp_col %in% names(gene_row))
    gwas_snp <- as.character(gene_row[[gwas_snp_col]])
  if (!is.null(gwas_pos_col) && gwas_pos_col %in% names(gene_row))
    gwas_pos <- suppressWarnings(as.numeric(gene_row[[gwas_pos_col]]))

  core_cols  <- c("RAPDB_Gene_ID", "Chromosomes", "Gene_Start", "Gene_End",
                  gwas_snp_col, gwas_pos_col)
  extra_cols <- setdiff(names(gene_row), core_cols)
  extra_dt   <- if (length(extra_cols)) {
    ex <- as.list(gene_row)[extra_cols]
    names(ex) <- paste0("Candidate_", extra_cols)
    data.table::as.data.table(ex)
  } else data.table::data.table()

  sr <- data.table::data.table(
    Gene                     = target_gene,
    Gene_Function            = gene_function,
    GFF3_Gene_ID             = ns_res$gene_gff_id,
    Transcript_ID            = ns_res$transcript_id,
    N_Transcripts_Considered = ns_res$n_transcripts_considered,
    Chr                      = as.character(gene_row$Chromosomes),
    MSU_Gene_ID              = msu$MSU_Gene_ID,
    MSU_Gene_Function        = msu$MSU_Gene_Function,
    MSU_Start                = msu$MSU_Start,
    MSU_End                  = msu$MSU_End,
    GWAS_SNP                 = gwas_snp,
    GWAS_Position            = gwas_pos,
    CDS_bp                   = ns_res$total_cds_bp,
    N_CDS_Exons              = ns_res$n_cds_exons,
    UTR5_bp                  = ns_res$utr5_bp,
    UTR3_bp                  = ns_res$utr3_bp,
    N_CDS_SNPs               = nrow(ns_res$all_cds_snps),
    N_NS_CDS_SNPs            = sum(ns_res$functional_snps$Region ==
                                     "CDS_nonsynonymous"),
    N_5UTR_SNPs              = sum(ns_res$functional_snps$Region == "5UTR"),
    N_3UTR_SNPs              = sum(ns_res$functional_snps$Region == "3UTR"),
    N_REF_Mismatch           = ns_res$n_ref_mismatch,
    N_Functional_SNPs_Raw    = nrow(ns_res$functional_snps),
    N_Functional_SNPs_QC     = gd$n_snps,
    N_Functional_SNPs_Final  = hap_mat$n_snps,
    Functional_Changes       = paste(ns_res$functional_snps$Variant_Label,
                                       collapse = ";"),
    N_Individuals_Genotyped  = hap_data$n_individuals,
    N_Haplotypes             = hap_data$n_haplotypes,
    Nei_Hd                   = hap_data$diversity$Hd,
    Neff                     = hap_data$diversity$Neff,
    PIC                      = hap_data$diversity$PIC,
    Test                     = test_res$test,
    PC_Corrected             = isTRUE(test_res$pc_corrected),
    PCs_Used                 = paste(avail_pcs, collapse = ","),
    P_value                  = test_res$p_value,
    Permutation_P            = test_res$permutation_p,
    Permutation_Statistic    = test_res$permutation_statistic,
    Adjusted_P               = NA_real_,
    Omega_Squared            = hbr_safe_num(test_res$omega_squared),
    CV_Percent               = hbr_safe_num(test_res$cv_percent),
    Superior_Haplotype       = sup$haplotype,
    Superior_NT_Seq          = sup$nt_sequence,
    Superior_AA_Seq          = sup$aa_sequence,
    Superior_Mean            = hbr_safe_num(sup$mean),
    Superior_SE              = hbr_safe_num(sup$se),
    Superior_CI_Lower        = hbr_safe_num(sup$ci_lower),
    Superior_CI_Upper        = hbr_safe_num(sup$ci_upper),
    Superior_Freq_Analysed_Pct = hbr_safe_num(sup$freq_analysed_pct),
    Superior_Freq_Panel_Pct    = hbr_safe_num(sup$freq_panel_pct),
    Second_Haplotype         = sup$second_haplotype,
    Descriptive_Letter       = sup$letter,
    Is_Superior              = sup$is_superior,
    Gate_A_Sig               = sup$gate_a_sig,
    Gate_B_Perm_Sig          = sup$gate_b_perm_sig,
    Gate_C_Effect            = sup$gate_c_effect,
    Gate_D_Threshold         = sup$gate_d_threshold,
    Improvement_Pct          = hbr_safe_num(sup$improvement_pct),
    Cohens_d                 = hbr_safe_num(sup$cohens_d),
    Diff_Mean                = hbr_safe_num(sup$diff_mean),
    Diff_CI_Lower            = hbr_safe_num(sup$diff_ci_lower),
    Diff_CI_Upper            = hbr_safe_num(sup$diff_ci_upper),
    Phenotype_Type           = phenotype_type)

  if (nrow(extra_dt)) sr <- cbind(sr, extra_dt)
  sr
}


# ----------------------------------------------------------------------------
# Gene object construction from CSV
# ----------------------------------------------------------------------------

#' Build a gene object for visualization
#'
#' @param gene_csv Path to `<gene>_viz.csv`.
#' @param gene_data In-memory data.frame/data.table.
#' @param individuals_csv Path to `<gene>_individuals.csv`.
#' @param individuals_data In-memory data.frame/data.table.
#' @param geno_file Path to a genotype file.
#' @param geno_matrix In-memory numeric matrix.
#' @param geno_format One of `"csv"`, `"raw"`, `"vcf"`, `"plink"`, `"hapmap"`.
#' @param trait_col Internal trait column name.
#' @param plink_path Path to PLINK.
#' @param verbose Print progress.
#'
#' @return A list structured like the pipeline's per-gene result.
#'
#' @noRd
hbr_build_gene_object <- function(
    gene_csv         = NULL,
    gene_data        = NULL,
    individuals_csv  = NULL,
    individuals_data = NULL,
    geno_file        = NULL,
    geno_matrix      = NULL,
    geno_format      = NULL,
    trait_col        = "trait",
    plink_path       = "plink",
    verbose          = TRUE) {

  if (!is.null(gene_data)) {
    dt <- data.table::as.data.table(gene_data)
  } else if (!is.null(gene_csv)) {
    if (!file.exists(gene_csv))
      stop("Gene CSV not found: ", gene_csv, call. = FALSE)
    dt <- data.table::fread(gene_csv, na.strings = c("", "NA"))
  } else {
    stop("Supply either `gene_data` or `gene_csv`.", call. = FALSE)
  }

  required <- c("gene_id", "chromosome", "gene_start", "gene_end",
                "haplotype", "n", "mean", "se", "ci_lower", "ci_upper",
                "nt_sequence", "aa_sequence")
  missing <- setdiff(required, names(dt))
  if (length(missing))
    stop("Missing required columns: ",
         paste(missing, collapse = ", "), call. = FALSE)

  meta <- list(
    gene_id            = as.character(dt$gene_id[1]),
    chromosome         = as.character(dt$chromosome[1]),
    gene_start         = as.numeric(dt$gene_start[1]),
    gene_end           = as.numeric(dt$gene_end[1]),
    p_value            = as.numeric(dt$p_value[1]),
    permutation_p      = as.numeric(dt$permutation_p[1]),
    superior_haplotype = as.character(dt$superior_haplotype[1]),
    is_superior        = as.logical(dt$is_superior[1]),
    cohens_d           = as.numeric(dt$cohens_d[1]),
    improvement_pct    = as.numeric(dt$improvement_pct[1]))

  hap_map <- data.table::data.table(
    Haplotype   = as.character(dt$haplotype),
    Count       = as.integer(dt$n),
    Frequency   = if ("frequency" %in% names(dt))
                    as.numeric(dt$frequency) else as.numeric(dt$n) / sum(dt$n),
    NT_Sequence = as.character(dt$nt_sequence),
    AA_Sequence = as.character(dt$aa_sequence))

  groups <- data.table::data.table(
    Haplotype = as.character(dt$haplotype),
    N         = as.integer(dt$n),
    Mean      = as.numeric(dt$mean),
    SE        = as.numeric(dt$se),
    CI_Lower  = as.numeric(dt$ci_lower),
    CI_Upper  = as.numeric(dt$ci_upper),
    Letter    = if ("letter" %in% names(dt))
                  as.character(dt$letter) else "a")

  snp_id_cols <- grep("^snp_[0-9]+_id$", names(dt), value = TRUE)
  if (!length(snp_id_cols))
    stop("No `snp_N_id` columns found.", call. = FALSE)

  snp_idx <- sort(as.integer(sub("^snp_([0-9]+)_id$", "\\1", snp_id_cols)))

  functional_snps <- data.table::rbindlist(lapply(snp_idx, function(i) {
    prefix  <- paste0("snp_", i, "_")
    get_col <- function(suffix, default = NA) {
      col <- paste0(prefix, suffix)
      if (col %in% names(dt)) dt[[col]][1] else default
    }
    region <- as.character(get_col("region"))
    data.table::data.table(
      SNP              = as.character(get_col("id")),
      POS              = as.numeric(get_col("pos")),
      REF              = as.character(get_col("ref")),
      ALT              = as.character(get_col("alt")),
      Region           = region,
      AA_Change        = as.character(get_col("aa_change")),
      Codon_Num        = as.numeric(get_col("codon_num")),
      Is_Nonsynonymous = !is.na(region) &
                           region == "CDS_nonsynonymous",
      Variant_Label    = if (!is.na(region) &&
                             region == "CDS_nonsynonymous")
                            as.character(get_col("aa_change"))
                          else paste0(region, ":",
                                       as.character(get_col("ref")), ">",
                                       as.character(get_col("alt"))))
  }))

  required_snps <- as.character(functional_snps$SNP)

  if (!is.null(individuals_data)) {
    ind <- data.table::as.data.table(individuals_data)
    data.table::setnames(ind, tolower(names(ind)))
    ind_dt <- data.table::data.table(
      ID        = as.character(ind$id),
      Haplotype = as.character(ind$haplotype))
    ind_dt[[trait_col]] <- as.numeric(ind$trait_value)
    required_ids <- as.character(ind$id)
  } else if (!is.null(individuals_csv) && file.exists(individuals_csv)) {
    ind <- data.table::fread(individuals_csv, na.strings = c("", "NA"))
    data.table::setnames(ind, tolower(names(ind)))
    ind_dt <- data.table::data.table(
      ID        = as.character(ind$id),
      Haplotype = as.character(ind$haplotype))
    ind_dt[[trait_col]] <- as.numeric(ind$trait_value)
    required_ids <- as.character(ind$id)
  } else {
    ind_dt <- data.table::rbindlist(lapply(seq_len(nrow(groups)), function(i) {
      d <- data.table::data.table(
        ID        = paste0("ind_", seq_len(groups$N[i]), "_",
                            groups$Haplotype[i]),
        Haplotype = groups$Haplotype[i])
      d[[trait_col]] <- stats::rnorm(groups$N[i],
                                       mean = groups$Mean[i],
                                       sd   = groups$SE[i] *
                                                sqrt(groups$N[i]))
      d
    }))
    required_ids <- NULL
  }

  geno_mat <- NULL

  if (!is.null(geno_matrix)) {
    m <- if (is.matrix(geno_matrix)) geno_matrix
         else as.matrix(data.table::as.data.table(geno_matrix)[, -1])
    storage.mode(m) <- "numeric"
    if (is.null(rownames(m)) && is.data.frame(geno_matrix))
      rownames(m) <- as.character(geno_matrix[[1]])

    keep <- intersect(required_snps, colnames(m))
    if (length(keep)) {
      m <- m[, keep, drop = FALSE]
      if (!is.null(required_ids)) {
        idx <- match(required_ids, rownames(m))
        keep_rows <- !is.na(idx)
        if (any(keep_rows)) m <- m[idx[keep_rows], , drop = FALSE]
      }
      geno_mat <- m
    } else if (verbose) {
      message("  In-memory genotype matrix has no columns matching ",
              "required SNPs -- skipping.")
    }
  } else if (!is.null(geno_file)) {
    fmt <- geno_format %||% .hbr_detect_geno_format(geno_file)
    if (verbose)
      message("  Genotype input: ", geno_file, " (format: ", fmt, ")")
    geno_mat <- tryCatch(
      .hbr_read_geno_subset(
        geno_file     = geno_file,
        geno_format   = fmt,
        req_snps    = required_snps,
        req_ids     = required_ids,
        plink_path    = plink_path,
        verbose       = verbose),
      error = function(e) {
        warning("Genotype loading failed for '", geno_file, "': ",
                conditionMessage(e), call. = FALSE)
        NULL
      })
  }

  if (!is.null(geno_mat) && verbose)
    message("  Genotype matrix: ", nrow(geno_mat), " x ", ncol(geno_mat))

  sup_row <- groups[Haplotype %in% meta$superior_haplotype]
  superior <- list(
    haplotype = meta$superior_haplotype,
    is_superior = meta$is_superior,
    mean       = if (nrow(sup_row)) sup_row$Mean[1] else NA_real_,
    se         = if (nrow(sup_row)) sup_row$SE[1] else NA_real_,
    ci_lower   = if (nrow(sup_row)) sup_row$CI_Lower[1] else NA_real_,
    ci_upper   = if (nrow(sup_row)) sup_row$CI_Upper[1] else NA_real_,
    n          = if (nrow(sup_row)) sup_row$N[1] else NA_integer_,
    letter     = if (nrow(sup_row)) sup_row$Letter[1] else NA_character_,
    cohens_d   = meta$cohens_d,
    improvement_pct = meta$improvement_pct,
    nt_sequence = hap_map$NT_Sequence[
      match(meta$superior_haplotype, hap_map$Haplotype)],
    aa_sequence = hap_map$AA_Sequence[
      match(meta$superior_haplotype, hap_map$Haplotype)],
    freq_analysed_pct = if (nrow(sup_row))
      round(sup_row$N[1] / sum(groups$N) * 100, 2) else NA_real_)

  list(
    source          = "csv",
    meta            = meta,
    trait_col       = trait_col,
    hap_map         = hap_map,
    groups          = groups,
    functional_snps = functional_snps,
    individuals_df  = ind_dt,
    geno            = geno_mat,
    superior        = superior,
    gene_annotation = data.table::data.table(
      RAPDB_Gene_ID = meta$gene_id,
      Chromosomes   = meta$chromosome,
      Gene_Start    = meta$gene_start,
      Gene_End      = meta$gene_end),
    p_value         = meta$p_value,
    permutation_p   = meta$permutation_p,
    is_region       = FALSE)
}

#' Load a gene visualization object from a CSV triplet
#'
#' Thin wrapper around [hbr_build_gene_object()] that reads the per-gene
#' `_viz.csv`, `_individuals.csv`, and optional genotype file, and returns
#' the standard gene object used by the plotting functions.
#'
#' @inheritParams hbr_build_gene_object
#'
#' @return A gene object.
#'
#' @noRd
hbr_load_gene_csv <- function(gene_csv,
                               individuals_csv = NULL,
                               geno_file       = NULL,
                               geno_format     = NULL,
                               trait_col       = "trait",
                               plink_path      = "plink",
                               verbose         = TRUE) {
  hbr_build_gene_object(
    gene_csv        = gene_csv,
    individuals_csv = individuals_csv,
    geno_file       = geno_file,
    geno_format     = geno_format,
    trait_col       = trait_col,
    plink_path      = plink_path,
    verbose         = verbose)
}

#' Resolve a summary CSV path
#'
#' @param summary_csv Directory or file path.
#'
#' @return The resolved file path.
#'
#' @noRd
.hbr_resolve_summary_csv_path <- function(summary_csv) {
  if (length(summary_csv) != 1L)
    stop("`summary_csv` must be a single path.", call. = FALSE)

  if (dir.exists(summary_csv)) {
    candidates <- c("publication_table_functional.csv",
                    "haplotype_summary_functional.csv",
                    "superior_genes_summary.csv")
    for (f in candidates) {
      p <- file.path(summary_csv, f)
      if (file.exists(p)) return(p)
    }
    stop("No summary CSV found in directory: ", summary_csv, call. = FALSE)
  }
  if (!file.exists(summary_csv))
    stop("Summary CSV not found: ", summary_csv, call. = FALSE)
  summary_csv
}

#' Convert a gene object to the internal ctx used by draw helpers
#'
#' @param res Gene object.
#'
#' @return A list for the plotting helpers.
#'
#' @noRd
.hbr_ctx_from_object <- function(res) {
  list(
    func_snps        = res$functional_snps,
    func_snps_used   = res$functional_snps,
    cds_snps         = res$functional_snps[
      res$functional_snps$Region %in%
        c("CDS_nonsynonymous", "CDS_synonymous")],
    utr_snps         = res$functional_snps[
      res$functional_snps$Region %in% c("5UTR", "3UTR")],
    hap_map          = res$hap_map,
    geno             = res$geno,
    test_data        = res$individuals_df,
    test_groups      = res$groups,
    superior         = res$superior,
    gene_annotation  = res$gene_annotation,
    network_data     = NULL)
}

#' Resolve a gene_list argument
#'
#' @param gene_list IDs, a CSV path, a directory, or a data.table.
#' @param summary_csv Optional summary CSV path.
#' @param per_gene_dir Override for the per-gene directory.
#' @param output_dir Fallback output directory.
#' @param verbose Print progress.
#'
#' @return A list with `genes`, `summary_csv`, `per_gene_dir`.
#'
#' @noRd
.hbr_resolve_gene_list <- function(gene_list,
                                   summary_csv  = NULL,
                                   per_gene_dir = NULL,
                                   output_dir   = NULL,
                                   verbose      = TRUE) {

  gene_col_candidates <- c("MSU7_Gene_ID", "RAPDB_Gene_ID",
                           "MSU_Gene_ID", "Gene_ID", "Gene", "gene_id")

  resolve_summary_file <- function(x) {
    if (dir.exists(x)) {
      for (f in c("publication_table_functional.csv",
                  "haplotype_summary_functional.csv",
                  "superior_genes_summary.csv")) {
        p <- file.path(x, f)
        if (file.exists(p)) return(p)
      }
      return(NULL)
    }
    if (file.exists(x) && grepl("\\.csv$", x, ignore.case = TRUE)) return(x)
    NULL
  }

  extract_ids_from_table <- function(dt) {
    gc <- intersect(gene_col_candidates, names(dt))[1]
    if (is.na(gc)) return(character(0))
    if ("Is_Superior" %in% names(dt)) {
      sup <- dt[dt$Is_Superior == TRUE, , drop = FALSE]
      if (nrow(sup))
        return(unique(as.character(sup[[gc]])))
    }
    unique(as.character(dt[[gc]]))
  }

  if (is.data.frame(gene_list)) {
    ids <- extract_ids_from_table(data.table::as.data.table(gene_list))
    return(list(genes = ids,
                summary_csv = summary_csv,
                per_gene_dir = per_gene_dir))
  }

  if (is.list(gene_list) && !is.data.frame(gene_list)) {
    return(list(genes = names(gene_list),
                summary_csv = summary_csv,
                per_gene_dir = per_gene_dir))
  }

  if (is.character(gene_list)) {

    if (length(gene_list) == 1L) {
      sf <- resolve_summary_file(gene_list)
      if (!is.null(sf)) {
        if (is.null(summary_csv)) summary_csv <- sf
        if (is.null(per_gene_dir)) per_gene_dir <- dirname(sf)
        dt  <- data.table::fread(sf)
        ids <- extract_ids_from_table(dt)
        if (verbose)
          message("`gene_list` is a source path -> ", length(ids),
                  " gene(s) extracted from ", basename(sf))
        return(list(genes = ids,
                    summary_csv = summary_csv,
                    per_gene_dir = per_gene_dir))
      }
      return(list(genes = gene_list,
                  summary_csv = summary_csv,
                  per_gene_dir = per_gene_dir))
    }

    is_path_like <- grepl("\\.csv$", gene_list, ignore.case = TRUE) |
      dir.exists(gene_list)
    if (any(is_path_like)) {
      ids <- character(0)
      for (el in gene_list) {
        if (dir.exists(el) || (file.exists(el) &&
                               grepl("\\.csv$", el, ignore.case = TRUE))) {
          sf2 <- resolve_summary_file(el)
          if (!is.null(sf2)) {
            if (is.null(summary_csv)) summary_csv <- sf2
            if (is.null(per_gene_dir)) per_gene_dir <- dirname(sf2)
            ids <- c(ids, extract_ids_from_table(data.table::fread(sf2)))
          }
        } else {
          ids <- c(ids, el)
        }
      }
      return(list(genes = unique(ids),
                  summary_csv = summary_csv,
                  per_gene_dir = per_gene_dir))
    }

    return(list(genes = gene_list,
                summary_csv = summary_csv,
                per_gene_dir = per_gene_dir))
  }

  list(genes = character(0),
       summary_csv = summary_csv,
       per_gene_dir = per_gene_dir)
}

#' Build a gene object from a summary CSV
#'
#' @param gene_id Gene ID.
#' @param summary_csv Summary CSV path.
#' @param per_gene_dir Directory with per-gene files.
#' @param trait_col Trait column name.
#' @param geno_file Optional genotype file.
#' @param geno_matrix Optional genotype matrix.
#' @param geno_format Optional format override.
#' @param plink_path Path to PLINK.
#' @param verbose Print progress.
#'
#' @return A gene object.
#'
#' @noRd
hbr_build_gene_object_from_summary <- function(
    gene_id,
    summary_csv,
    per_gene_dir = NULL,
    trait_col    = "NT",
    geno_file    = NULL,
    geno_matrix  = NULL,
    geno_format  = NULL,
    plink_path   = "plink",
    verbose      = TRUE) {

  if (!file.exists(summary_csv))
    stop("Summary CSV not found: ", summary_csv, call. = FALSE)

  if (is.null(per_gene_dir)) per_gene_dir <- dirname(summary_csv)

  s <- data.table::fread(summary_csv)
  gene_col <- intersect(
    c("RAPDB_Gene_ID", "MSU7_Gene_ID", "MSU_Gene_ID", "Gene_ID", "Gene"),
    names(s))[1]
  if (is.na(gene_col))
    stop("Cannot detect gene-ID column in ", summary_csv, call. = FALSE)

  row <- s[as.character(get(gene_col)) == gene_id]
  if (!nrow(row))
    stop("Gene '", gene_id, "' not in ", summary_csv, call. = FALSE)

  find_file <- function(suffix, required = TRUE) {
    f <- file.path(per_gene_dir, paste0(gene_id, suffix))
    if (!file.exists(f)) {
      if (required)
        stop("Missing per-gene file: ", f, call. = FALSE)
      return(NULL)
    }
    f
  }

  hap_map_file   <- find_file("_hap_map.csv")
  groups_file    <- find_file("_groups.csv")
  func_snps_file <- find_file("_functional_snps.csv")
  ind_file       <- find_file("_individuals.csv", required = FALSE)

  hap_map   <- data.table::fread(hap_map_file)
  groups    <- data.table::fread(groups_file)
  func_snps <- data.table::fread(func_snps_file)

  ind <- NULL
  if (!is.null(ind_file)) {
    ind <- data.table::fread(ind_file)
    if (!"Haplotype" %in% names(ind))
      stop("No Haplotype column in ", ind_file, call. = FALSE)
    if (!trait_col %in% names(ind))
      stop("Trait column '", trait_col, "' not in ", ind_file,
           ". Available: ", paste(names(ind), collapse = ", "),
           call. = FALSE)
  }

  n_hap <- nrow(hap_map)
  get_meta <- function(col, default = NA) {
    if (col %in% names(row)) row[[col]][1] else default
  }

  letter_col <- if ("Letter" %in% names(groups)) "Letter"
                else if ("Significance_Group" %in% names(groups))
                  "Significance_Group"
                else NULL

  gene_viz <- data.table::data.table(
    gene_id            = rep(gene_id, n_hap),
    chromosome         = rep(as.character(get_meta(
      "Chr", get_meta("Chromosome", NA_character_))), n_hap),
    gene_start         = rep(as.numeric(get_meta("Gene_Start", NA_real_)),
                              n_hap),
    gene_end           = rep(as.numeric(get_meta("Gene_End",   NA_real_)),
                              n_hap),
    p_value            = rep(as.numeric(get_meta("P_value", NA_real_)),
                              n_hap),
    permutation_p      = rep(as.numeric(get_meta("Permutation_P",
                                                  NA_real_)), n_hap),
    superior_haplotype = rep(as.character(get_meta("Superior_Haplotype",
                                                    NA_character_)), n_hap),
    is_superior        = rep(as.logical(get_meta("Is_Superior", TRUE)), n_hap),
    cohens_d           = rep(as.numeric(get_meta("Cohens_d", NA_real_)),
                              n_hap),
    improvement_pct    = rep(as.numeric(get_meta("Improvement_Pct",
                                                  NA_real_)), n_hap),
    haplotype   = as.character(hap_map$Haplotype),
    n           = as.integer(hap_map$Count),
    frequency   = as.numeric(hap_map$Frequency),
    nt_sequence = as.character(hap_map$NT_Sequence),
    aa_sequence = as.character(hap_map$AA_Sequence))

  groups_sel <- groups[, c("Haplotype", "Mean", "SE", "CI_Lower",
                           "CI_Upper", letter_col), with = FALSE]
  if (!is.null(letter_col) && letter_col != "Letter")
    data.table::setnames(groups_sel, letter_col, "Letter")
  if (is.null(letter_col))
    groups_sel[, Letter := "a"]

  gene_viz <- merge(gene_viz, groups_sel, by = "Haplotype", all.x = TRUE)

  if ("Genomic_Order" %in% names(func_snps))
    data.table::setorder(func_snps, Genomic_Order)

  for (j in seq_len(nrow(func_snps))) {
    snp <- func_snps[j]
    prefix <- paste0("snp_", j, "_")
    gene_viz[[paste0(prefix, "id")]]        <- as.character(snp$SNP)
    gene_viz[[paste0(prefix, "pos")]]       <- as.numeric(snp$POS)
    gene_viz[[paste0(prefix, "ref")]]       <- as.character(snp$REF)
    gene_viz[[paste0(prefix, "alt")]]       <- as.character(snp$ALT)
    gene_viz[[paste0(prefix, "region")]]    <- as.character(snp$Region)
    gene_viz[[paste0(prefix, "aa_change")]] <- if ("AA_Change" %in% names(snp))
      as.character(snp$AA_Change) else NA_character_
    gene_viz[[paste0(prefix, "codon_num")]] <- if ("Codon_Num" %in% names(snp))
      as.numeric(snp$Codon_Num) else NA_real_
  }

  individuals_data <- if (!is.null(ind)) {
    data.table::data.table(
      id          = as.character(ind$ID),
      haplotype   = as.character(ind$Haplotype),
      trait_value = as.numeric(ind[[trait_col]]))
  } else NULL

  hbr_build_gene_object(
    gene_data        = gene_viz,
    individuals_data = individuals_data,
    geno_file        = geno_file,
    geno_matrix      = geno_matrix,
    geno_format      = geno_format,
    trait_col        = trait_col,
    plink_path       = plink_path,
    verbose          = verbose)
}

# R/hgsi-io.R
# Reads HaploBreedR output (in memory or from the trait output directory) and
# rebuilds per-gene haplotype codes for all tested genes.

#' Normalise a trait direction
#'
#' Maps pipeline direction labels to `"higher"` or `"lower"`; stops if the value is missing or unknown.
#' @param x Direction as stored in `pipeline_summary$Direction`.
#' @return `"higher"` or `"lower"`.
#' @noRd
.hgsi_norm_direction <- function(x) {
  x <- tolower(trimws(as.character(x)[1]))
  if (!is.na(x) && x %in% c("higher", "high", "max", "increase", "up")) return("higher")
  if (!is.na(x) && x %in% c("lower", "low", "min", "decrease", "down")) return("lower")
  stop("pipeline direction missing or unrecognised ('", x,
       "'); expected higher/lower in pipeline_summary$Direction", call. = FALSE)
}

#' Split a delimited string into unique tokens
#'
#' Drops empty, `NA`, `"none"` and `"null"` entries.
#' @param x Character vector (comma, semicolon, plus or space separated).
#' @return Character vector of unique tokens.
#' @noRd
.hgsi_tokens <- function(x) {
  x <- unlist(strsplit(as.character(x), "[,;+ ]+"))
  x <- x[!is.na(x) & nzchar(x) & !tolower(x) %in% c("na", "none", "null")]
  unique(x)
}

#' Read a CSV from a directory if present
#'
#' @param dir Directory, or `NULL`.
#' @param file File name.
#' @return A data.frame, or `NULL` if the file does not exist.
#' @noRd
.hgsi_read_csv <- function(dir, file) {
  p <- file.path(dir, file)
  if (!is.null(dir) && !is.na(dir) && file.exists(p))
    data.table::fread(p, data.table = FALSE) else NULL
}

# ---- pipeline summary and gate parameters -----------------------------------

#' Get the pipeline summary row
#'
#' Uses `pipeline_result$pipeline_summary`, else `pipeline_summary.csv` in `pipeline_dir`.
#' @param pipeline_result In-memory pipeline result or `NULL`.
#' @param pipeline_dir Pipeline output directory or `NULL`.
#' @return One-row data.frame.
#' @noRd
.hgsi_read_summary <- function(pipeline_result, pipeline_dir) {
  ps <- pipeline_result$pipeline_summary
  if (is.null(ps)) ps <- .hgsi_read_csv(pipeline_dir, "pipeline_summary.csv")
  if (is.null(ps) || !nrow(ps))
    stop("pipeline_summary not found (result or ", pipeline_dir,
         "/pipeline_summary.csv)", call. = FALSE)
  as.data.frame(ps)
}

#' Superior-haplotype gate settings of a pipeline run
#'
#' Reads alpha, minimum effect size, minimum group size, minimum haplotype frequency, MAF, trait threshold and FDR requirement from the summary row; unavailable values get the pipeline defaults and are reported. `trait_threshold` is read per trait.
#' @param ps Summary row from [.hgsi_read_summary()].
#' @param direction `"higher"` or `"lower"`.
#' @param user Optional override list (`alpha`, `min_d`/`min_cohens_d`, `min_ind`/`min_ind_per_hap`, `min_hap_freq`, `maf_min`, `trait_threshold`, `require_fdr`, `hedges`).
#' @return List with `gate` (settings plus `direction` and `hedges`) and `defaulted` (names filled with defaults).
#' @noRd
.hgsi_gate_from_summary <- function(ps, direction, user = NULL) {
  num <- function(col) {
    v <- suppressWarnings(as.numeric(ps[[col]][1]))
    if (length(v) && is.finite(v)) v else NA_real_
  }
  rf <- as.logical(ps[["Require_FDR_Gate"]][1])
  g <- list(alpha = num("Alpha_Gate"), min_d = num("Min_Cohens_D_Gate"),
            min_ind = num("Min_Ind_Per_Hap_Gate"),
            min_hap_freq = num("Min_Hap_Freq_Gate"),
            maf_min = num("MAF_Min_Gate"),
            trait_threshold = num("Trait_Threshold_Gate"),
            require_fdr = if (length(rf) && !is.na(rf)) rf else TRUE,
            hedges = TRUE)   # Hedges' g; use gate_params = list(hedges = FALSE) for older runs
  if (!is.null(user)) {
    if (is.null(user$min_d) && !is.null(user$min_cohens_d)) user$min_d <- user$min_cohens_d
    if (is.null(user$min_ind) && !is.null(user$min_ind_per_hap))
      user$min_ind <- user$min_ind_per_hap
    nm <- intersect(names(user), names(g))
    g[nm] <- user[nm]
  }
  notes <- character()
  fill <- function(name, value) {
    if (!is.finite(g[[name]])) { g[[name]] <<- value; notes <<- c(notes, name) }
  }
  fill("alpha", 0.05); fill("min_d", 0.3); fill("min_ind", 3)
  fill("min_hap_freq", 0); fill("maf_min", 0.05)
  g$min_ind <- as.integer(g$min_ind)
  g$direction <- direction
  list(gate = g, defaulted = notes)
}

# ---- tested genes, superior genes, superior haplotypes ---------------------

#' Tested-gene table of a pipeline run
#'
#' Uses `pipeline_result$summary`, else `publication_table_functional.csv`, else `haplotype_summary_functional.csv`. Adds `gene` and logical `is_superior`.
#' @param pipeline_result In-memory pipeline result or `NULL`.
#' @param out_dir Pipeline output directory.
#' @param ps Summary row (gives `Annotation_Source`, which names the gene-ID column).
#' @return List with `table` and `source`.
#' @noRd
.hgsi_read_tested <- function(pipeline_result, out_dir, ps) {
  tb <- pipeline_result$summary; src <- "pipeline_result$summary"
  if (is.null(tb) || !NROW(tb)) {
    tb <- NULL
    for (f in c("publication_table_functional.csv",
                "haplotype_summary_functional.csv")) {
      tb <- .hgsi_read_csv(out_dir, f)
      if (!is.null(tb) && nrow(tb)) { src <- f; break }
    }
  }
  if (is.null(tb) || !nrow(tb))
    stop("tested-gene table not found (result$summary or ", out_dir,
         "/publication_table_functional.csv)", call. = FALSE)
  tb <- as.data.frame(tb)
  annot <- as.character(ps$Annotation_Source[1])
  gc <- intersect(c(if (!is.na(annot) && annot != "CUSTOM")
                      paste0(annot, "_Gene_ID") else "Gene_ID",
                    "Gene_ID", "MSU7_Gene_ID", "RAPDB_Gene_ID", "Gene"),
                  names(tb))[1]
  if (is.na(gc) || !"Is_Superior" %in% names(tb))
    stop("tested-gene table needs a gene-ID column and Is_Superior; found: ",
         paste(names(tb), collapse = ", "), call. = FALSE)
  tb$gene <- as.character(tb[[gc]])
  tb$is_superior <- as.logical(tb$Is_Superior) %in% TRUE
  list(table = tb[!is.na(tb$gene) & nzchar(tb$gene), , drop = FALSE], source = src)
}

#' Functional SNPs per gene
#'
#' @param fs Functional SNP table with `gene`, `SNP` and optionally `Genomic_Order`.
#' @return Named list: gene -> SNP IDs in genomic order.
#' @noRd
.hgsi_gene_snps <- function(fs) {
  if (!all(c("gene", "SNP") %in% names(fs)))
    stop("functional SNP table needs `gene` and `SNP` columns", call. = FALSE)
  fs <- as.data.frame(fs)
  if ("Genomic_Order" %in% names(fs)) fs <- fs[order(fs$Genomic_Order), , drop = FALSE]
  split(as.character(fs$SNP), as.character(fs$gene))
}


#' Union of phenotypes from the per-gene individuals CSVs
#'
#' Fallback used when neither the `phenotypes` argument, nor
#' `pipeline_result$phenotype_data`, nor `phenotype_data_used.csv` is
#' available. Reads every `<gene>_individuals.csv` in the trait output
#' directory (and in `superior_genes/`), pulls the ID and the trait
#' column, and returns one row per individual.
#'
#' @param out_dir Trait output directory.
#' @param trait_col Trait column name as written by the pipeline.
#' @param verbose Print the number of files and individuals.
#'
#' @return A data.frame with columns `ID` and `trait_col`, or `NULL`.
#'
#' @noRd
.hgsi_phenotypes_from_individuals <- function(out_dir, trait_col,
                                              verbose = TRUE) {
  if (is.null(out_dir) || !dir.exists(out_dir)) return(NULL)
  dirs <- unique(c(file.path(out_dir, "superior_genes"), out_dir))
  files <- unique(unlist(lapply(dirs, function(d) {
    if (!dir.exists(d)) return(character(0))
    list.files(d, pattern = "_individuals\\.csv$", full.names = TRUE)
  })))
  if (!length(files)) return(NULL)
  
  parts <- lapply(files, function(f) {
    d <- tryCatch(data.table::fread(f, data.table = FALSE),
                  error = function(e) NULL)
    if (is.null(d)) return(NULL)
    idc <- intersect(c("ID", "id", "IID"), names(d))[1]
    yc  <- intersect(c(trait_col, "trait_value", "Trait"), names(d))[1]
    if (is.na(idc) || is.na(yc)) return(NULL)
    data.frame(ID = hbr_clean_id(d[[idc]]),
               y  = suppressWarnings(as.numeric(d[[yc]])),
               stringsAsFactors = FALSE)
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) return(NULL)
  
  out <- data.table::rbindlist(parts, fill = TRUE)
  out <- out[!is.na(ID) & !is.na(y)]
  # If a gene's haplotype assignment puts an individual in more than one
  # row (unlikely but possible after ID cleaning), keep the first.
  out <- out[!duplicated(ID)]
  
  if (verbose)
    message(sprintf(
      "    phenotypes from %d per-gene individuals CSV(s): %d individuals",
      length(parts), nrow(out)))
  
  df <- data.frame(ID = out$ID, stringsAsFactors = FALSE)
  df[[trait_col]] <- out$y
  df
}
# ---- main loader ------------------------------------------------------------

#' Load one trait from HaploBreedR output
#'
#' Collects direction, gate settings, tested and superior genes, functional SNPs, phenotypes and the PCs the pipeline used (none if it ran with `pc_cols = NULL`).
#' @details
#' Sources, in priority order: pipeline summary (`pipeline_summary` | `pipeline_summary.csv`); tested genes and superior haplotypes (`summary` | `publication_table_functional.csv`); functional SNPs (`functional_snps_all` | `functional_snps_all_tested_genes.csv`); phenotypes and PCs (`phenotypes` argument | `phenotype_data` | `phenotype_data_used.csv`, which has no PCs).
#' @param trait_name Trait label.
#' @param pipeline_result In-memory pipeline result or `NULL`.
#' @param pipeline_dir Pipeline output directory or `NULL`.
#' @param phenotypes Optional wide data.frame (ID, trait, PC columns).
#' @param phenotype_id_col ID column of `phenotypes`.
#' @param gate_params_user Optional override list for the gate settings.
#' @param verbose Print a one-line summary.
#' @return List: `direction`, `gate`, `phenotype`, `pcs`, `pc_cols`, `gene_snps`, `tested`, `superior`, `superior_label`, `plink_prefix`, `out_dir`, `is_gwas`, `notes`, `provenance`.
#' @noRd
.hgsi_load_trait <- function(trait_name, pipeline_result = NULL,
                             pipeline_dir = NULL, phenotypes = NULL,
                             phenotype_id_col = "ID", gate_params_user = NULL,
                             verbose = TRUE) {
  if (is.null(pipeline_result) && is.null(pipeline_dir))
    stop("Supply a pipeline result or its output directory for trait '",
         trait_name, "'.", call. = FALSE)
  ps <- .hgsi_read_summary(pipeline_result, pipeline_dir)
  out_dir <- pipeline_dir %||% as.character(ps$Output_Dir[1])
  notes <- character()

  direction <- .hgsi_norm_direction(ps$Direction[1])
  trait_col <- as.character(ps$Trait[1]) %||% trait_name
  mode      <- as.character(ps$Mode[1])
  plink     <- as.character(ps$PLINK_Prefix[1])

  gg   <- .hgsi_gate_from_summary(ps, direction, gate_params_user)
  gate <- gg$gate
  if (length(gg$defaulted))
    notes <- c(notes, sprintf(
      "%s: gate settings not stored in the pipeline summary, defaults used for: %s",
      trait_name, paste(gg$defaulted, collapse = ", ")))

  # ---- tested / superior genes ----------------------------------------------
  tt <- .hgsi_read_tested(pipeline_result, out_dir, ps)
  tb <- tt$table

  # ---- functional SNPs ------------------------------------------------------
  fs <- pipeline_result$functional_snps_all; fs_source <- "pipeline_result"
  if (is.null(fs) || !NROW(fs)) {
    fs <- .hgsi_read_csv(out_dir, "functional_snps_all_tested_genes.csv")
    fs_source <- "functional_snps_all_tested_genes.csv"
  }
  if (is.null(fs) || !NROW(fs))
    stop("functional SNPs of tested genes not found for trait '", trait_name,
         "' (functional_snps_all or functional_snps_all_tested_genes.csv)",
         call. = FALSE)
  gene_snps <- .hgsi_gene_snps(fs)
  gene_snps <- gene_snps[intersect(unique(tb$gene), names(gene_snps))]
  tested    <- names(gene_snps)
  sup_rows  <- tb[tb$is_superior & tb$gene %in% tested, , drop = FALSE]
  superior  <- unique(sup_rows$gene)
  lost <- setdiff(tb$gene[tb$is_superior], tested)
  if (length(lost))
    notes <- c(notes, sprintf("%s: %d superior gene(s) have no functional SNPs on file and are ignored",
                              trait_name, length(lost)))
  sup_label <- if ("Superior_Haplotype" %in% names(sup_rows))
    stats::setNames(as.character(sup_rows$Superior_Haplotype), sup_rows$gene)
  else character(0)

  # ---- phenotypes and PCs ---------------------------------------------------
  pdf <- NULL; ph_source <- NA_character_
  if (!is.null(phenotypes)) {
    pdf <- as.data.frame(phenotypes)
    ph_source <- "phenotypes argument"
  } else if (!is.null(pipeline_result$phenotype_data)) {
    pdf <- as.data.frame(pipeline_result$phenotype_data)
    ph_source <- "pipeline_result$phenotype_data"
  } else {
    pdf <- .hgsi_read_csv(out_dir, "phenotype_data_used.csv")
    if (!is.null(pdf)) ph_source <- "phenotype_data_used.csv"
  }
  # Fallback: union of per-gene individuals CSVs written by the pipeline
  if (is.null(pdf)) {
    pdf <- .hgsi_phenotypes_from_individuals(out_dir, trait_col,
                                             verbose = verbose)
    if (!is.null(pdf)) ph_source <- "per-gene individuals CSV"
  }
  if (is.null(pdf)) stop("No phenotypes found for trait '", trait_name,
                         "'. Supply `phenotypes`.", call. = FALSE)
  idc <- intersect(c(phenotype_id_col, "ID", "id"), names(pdf))[1]
  # phenotype_data_used.csv holds the pipeline's own ID column (ID, then trait)
  if (is.na(idc) && identical(ph_source, "phenotype_data_used.csv")) idc <- names(pdf)[1]
  yc  <- intersect(c(trait_col, trait_name), names(pdf))[1]
  if (is.na(idc) || is.na(yc))
    stop("phenotype table needs an ID column and a '", trait_col, "' column",
         call. = FALSE)
  ph <- data.frame(ID = hbr_clean_id(pdf[[idc]]), y = as.numeric(pdf[[yc]]),
                   stringsAsFactors = FALSE)
  keep <- !is.na(ph$ID) & !is.na(ph$y) & !duplicated(ph$ID)
  ph <- ph[keep, , drop = FALSE]

  # PCs: those the pipeline actually used (PCs_Used), else its PC_Cols_Gate.
  # NULL / empty -> no PC adjustment (as in the pipeline).
  used_tbl <- if ("PCs_Used" %in% names(tb)) .hgsi_tokens(tb$PCs_Used) else NULL
  want <- if (!is.null(used_tbl)) used_tbl else .hgsi_tokens(ps$PC_Cols_Gate[1])
  pcs <- NULL; pc_cols <- character(0)
  if (length(want)) {
    have <- intersect(want, names(pdf))
    if (length(have) < length(want) && !is.null(used_tbl))
      stop("the pipeline used PCs (", paste(want, collapse = ", "),
           ") but they are not in the phenotype source '", ph_source,
           "'. Supply `phenotypes` with ID, trait and these PC columns.",
           call. = FALSE)
    if (length(have)) {
      pc_cols <- have
      pcs <- data.frame(ID = hbr_clean_id(pdf[[idc]]),
                        as.data.frame(lapply(pdf[have], as.numeric)),
                        stringsAsFactors = FALSE)
      pcs <- pcs[!duplicated(pcs$ID), , drop = FALSE]
    } else notes <- c(notes, sprintf(
      "%s: PC columns (%s) not available; gate replayed without PCs",
      trait_name, paste(want, collapse = ", ")))
  }

  if (verbose)
    message(sprintf("  %s: direction=%s | %d tested, %d superior | %d phenotyped | PCs=%s | %s",
                    trait_name, direction, length(tested), length(superior),
                    nrow(ph), if (length(pc_cols)) paste(pc_cols, collapse = ",") else "none",
                    tt$source))

  list(trait = trait_name, trait_col = trait_col, direction = direction,
       mode = mode, is_gwas = isTRUE(grepl("gwas", mode, ignore.case = TRUE)),
       plink_prefix = plink, out_dir = out_dir,
       phenotype = ph, pcs = pcs, pc_cols = pc_cols,
       gene_snps = gene_snps, tested = tested, superior = superior,
       superior_label = sup_label,
       superior_results = pipeline_result$superior_results %||% list(),
       gate = gate, notes = notes,
       provenance = list(phenotypes = ph_source, tested_genes = tt$source,
                         functional_snps = fs_source))
}

# ---- haplotype codes --------------------------------------------------------

#' Haplotype codes for all tested genes
#'
#' Mirrors the pipeline's haplotype construction: heterozygous or missing calls are `"N"`; SNPs with MAF below `maf_min` or more than 5% `"N"` are dropped; individuals with more than 10% `"N"` get `NA`. The code is the `"-"`-joined dosage string; allele orientation does not change the partition.
#' @param plink_prefix PLINK prefix.
#' @param gene_snps Named list gene -> SNP IDs.
#' @param ids Cleaned individual IDs in `.fam` order.
#' @param snps SNP IDs in `.bim` order.
#' @param maf_min Minimum minor allele frequency.
#' @param verbose Print the SNP mapping rate.
#' @return Character matrix, individuals x genes, row names = `ids`.
#' @noRd
.hgsi_build_codes <- function(plink_prefix, gene_snps, ids, snps,
                              maf_min = 0.05, verbose = TRUE) {
  if (!requireNamespace("BEDMatrix", quietly = TRUE))
    stop("BEDMatrix required. install.packages('BEDMatrix')", call. = FALSE)
  G <- suppressMessages(BEDMatrix::BEDMatrix(plink_prefix, simple_names = TRUE))
  allsnp <- unique(unlist(gene_snps, use.names = FALSE))
  idx <- match(allsnp, snps); mapped <- mean(!is.na(idx))
  if (verbose)
    message(sprintf("    functional SNPs mapped to PLINK: %.1f%% (%d)",
                    100 * mapped, sum(!is.na(idx))))
  if (mapped < 0.5)
    stop("fewer than 50% of functional SNP IDs found in the PLINK file; ",
         "check ID format", call. = FALSE)
  keep_idx <- sort(unique(idx[!is.na(idx)]))
  X <- G[, keep_idx, drop = FALSE]; storage.mode(X) <- "double"
  colnames(X) <- snps[keep_idx]; rownames(X) <- ids
  cols <- list()
  for (g in names(gene_snps)) {
    s <- intersect(gene_snps[[g]], colnames(X))
    if (!length(s)) next
    m <- X[, s, drop = FALSE]
    p <- colMeans(m, na.rm = TRUE) / 2
    m[m == 1] <- NA
    ok <- is.finite(p) & pmin(p, 1 - p) >= maf_min & colMeans(is.na(m)) <= 0.05
    if (!any(ok)) next
    m <- m[, ok, drop = FALSE]
    parts <- lapply(seq_len(ncol(m)), function(j)
      ifelse(is.na(m[, j]), "N", as.character(m[, j])))
    key <- do.call(paste, c(parts, sep = "-"))
    key[rowMeans(is.na(m)) > 0.10] <- NA_character_
    cols[[g]] <- key
  }
  if (!length(cols)) stop("no haplotype codes could be built", call. = FALSE)
  codes <- do.call(cbind, cols)
  rownames(codes) <- ids
  codes
}

#' Relabel haplotype codes to the pipeline's HapN names
#'
#' The engine's codes are joined dosage patterns (e.g. `"0-2-0"`). This
#' helper maps each code to the pipeline's haplotype label (`Hap1`,
#' `Hap2`, ...) using the per-gene haplotype assignment stored in the
#' pipeline result or its `<gene>_individuals.csv` files. Genes for
#' which no pipeline assignment exists (tested but not superior) keep
#' their dosage codes. The partition of individuals is unchanged.
#'
#' @param codes Character matrix, individuals x genes, from
#'   [.hgsi_build_codes()].
#' @param ids Individual IDs (row names of `codes`).
#' @param ld Output of [.hgsi_load_trait()].
#' @return The relabelled codes matrix.
#' @noRd
.hgsi_relabel_codes <- function(codes, ids, ld) {
  for (g in colnames(codes)) {
    hp <- NULL
    # 1. Try the pipeline result held in memory
    if (!is.null(ld$superior_results) &&
        !is.null(ld$superior_results[[g]])) {
      hp <- ld$superior_results[[g]]$hap_data$haplotypes
    }
    # 2. Try the per-gene files written by the pipeline
    if (is.null(hp) || !all(c("ID", "Haplotype") %in% names(hp))) {
      for (d in c(file.path(ld$out_dir, "superior_genes"), ld$out_dir)) {
        f <- file.path(d, paste0(g, "_individuals.csv"))
        if (file.exists(f)) {
          hp <- tryCatch(data.table::fread(f, data.table = FALSE),
                         error = function(e) NULL)
          if (!is.null(hp) && all(c("ID", "Haplotype") %in% names(hp))) break
          hp <- NULL
        }
      }
    }
    if (is.null(hp) || !all(c("ID", "Haplotype") %in% names(hp))) next
    
    ref  <- as.character(hp$Haplotype)[match(ids, hbr_clean_id(hp$ID))]
    mine <- codes[, g]
    ok   <- !is.na(mine) & !is.na(ref)
    if (sum(ok) < 10L) next
    
    # For each unique code, use the modal HapN among the individuals the
    # pipeline assigned to that code. If a code maps to exactly one HapN,
    # the modal value is that HapN. Ties are broken by first occurrence.
    map <- tapply(ref[ok], mine[ok], function(v)
      names(which.max(table(v))))
    
    codes[, g] <- ifelse(is.na(mine), NA_character_,
                         ifelse(mine %in% names(map), map[mine], mine))
  }
  codes
}

# ---- fidelity of the rebuilt codes -------------------------------------------

#' Compare rebuilt haplotypes with the pipeline's
#'
#' For superior genes, reports partition agreement (purity both ways) and the Jaccard overlap between the pipeline's superior-haplotype carriers and the replayed best haplotype; the latter also checks that the same direction was used.
#' @param codes Matrix from [.hgsi_build_codes()].
#' @param ids Individual IDs.
#' @param ld Output of [.hgsi_load_trait()].
#' @param best_hap Named character: gene -> replayed best haplotype code.
#' @return Data.frame: `gene`, `n`, `agreement`, `superior_overlap`.
#' @noRd
.hgsi_check_fidelity <- function(codes, ids, ld, best_hap) {
  out <- list()
  for (g in intersect(ld$superior, colnames(codes))) {
    hp <- ld$superior_results[[g]]$hap_data$haplotypes
    if (is.null(hp))
      for (d in c(file.path(ld$out_dir, "superior_genes"), ld$out_dir)) {
        hp <- .hgsi_read_csv(d, paste0(g, "_individuals.csv"))
        if (!is.null(hp)) break
      }
    if (is.null(hp)) next
    hp <- as.data.frame(hp)
    if (!all(c("ID", "Haplotype") %in% names(hp))) next
    ref  <- as.character(hp$Haplotype)[match(ids, hbr_clean_id(hp$ID))]
    mine <- codes[, g]
    ok <- !is.na(ref) & !is.na(mine)
    if (sum(ok) < 10L) next
    tab <- table(mine[ok], ref[ok])
    agree <- min(sum(apply(tab, 1, max)), sum(apply(tab, 2, max))) / sum(tab)
    lab <- unname(ld$superior_label[g]); bh <- unname(best_hap[g])
    jac <- NA_real_
    if (length(lab) && length(bh) && !is.na(lab) && !is.na(bh)) {
      a <- which(!is.na(mine) & ref %in% lab)
      b <- which(!is.na(ref) & !is.na(mine) & mine == bh)
      jac <- length(intersect(a, b)) / max(1L, length(union(a, b)))
    }
    out[[length(out) + 1L]] <- data.frame(gene = g, n = sum(ok),
      agreement = agree, superior_overlap = jac, stringsAsFactors = FALSE)
  }
  if (!length(out)) return(data.frame(gene = character(), n = integer(),
                                       agreement = numeric(),
                                       superior_overlap = numeric()))
  do.call(rbind, out)
}

#' Label how honest the CV evaluation is
#'
#' @param is_gwas Pipeline ran in GWAS mode.
#' @param universe_all_tested Codes cover at least 80% of the tested genes.
#' @param concordance_ok Replayed and pipeline superior genes agree well enough.
#' @param fidelity_ok Rebuilt haplotypes match the pipeline's.
#' @return One of `"conditional_on_gwas_candidates"`, `"conditional_on_discovered_set"`, `"approximate_replay"`, `"honest_discovery_in_fold"`.
#' @noRd
.hgsi_eval_label <- function(is_gwas, universe_all_tested, concordance_ok,
                             fidelity_ok) {
  if (isTRUE(is_gwas)) return("conditional_on_gwas_candidates")
  if (!isTRUE(universe_all_tested)) return("conditional_on_discovered_set")
  if (!isTRUE(concordance_ok) || !isTRUE(fidelity_ok)) return("approximate_replay")
  "honest_discovery_in_fold"
}

#' Carrier status of the superior haplotype at a gene
#'
#' Reads the per-gene `<gene>_individuals.csv` written by
#' [haplobreed_analysis()] and returns, for every individual in the HGSI
#' result, whether they carry the pipeline's superior haplotype at that
#' gene. Used to colour per-individual figures by carrier status.
#'
#' @param x An `hgsi_advanced` object.
#' @param trait Trait label.
#' @param gene Gene ID.
#' @param trait_dir Path to the trait's pipeline output directory. If
#'   `NULL`, the helper looks in `x$trait_dirs[[trait]]`.
#' @param superior_haplotype Optional. When `NULL` (default), the
#'   helper reads it from `superior_genes_summary.csv` in `trait_dir`.
#'
#' @return A list with `carrier` (logical, length `length(x$ids)`),
#'   `haplotype` (character, the assigned haplotype per individual),
#'   `superior_haplotype` (character), and `gene` and `trait`.
#' @examples
#' toy <- system.file("extdata", "hgsi_toy", package = "HaploBreedR")
#' td  <- file.path(toy, "trait_YLD")
#' if (nzchar(toy) && dir.exists(td)) {
#'   x <- structure(
#'     list(ids = sprintf("L%03d", 1:150),
#'          traits = "YLD",
#'          trait_dirs = c(YLD = td)),
#'     class = "hgsi_advanced")
#'   res <- try(hbr_gene_carriers(
#'     x, trait = "YLD", gene = "GENE04", trait_dir = td))
#' }
#' @export
hbr_gene_carriers <- function(x, trait, gene,
                              trait_dir = NULL,
                              superior_haplotype = NULL) {
  
  if (!inherits(x, "hgsi_advanced"))
    stop("`x` must be an hgsi_advanced object.", call. = FALSE)
  if (!trait %in% x$traits)
    stop("Trait '", trait, "' not in x$traits.", call. = FALSE)
  
  if (is.null(trait_dir)) {
    if (!is.null(x$trait_dirs) && trait %in% names(x$trait_dirs))
      trait_dir <- x$trait_dirs[[trait]]
    else stop("Supply `trait_dir` or populate `x$trait_dirs`.", call. = FALSE)
  }
  
  # Locate the per-gene individuals file
  candidates <- c(
    file.path(trait_dir, "superior_genes", paste0(gene, "_individuals.csv")),
    file.path(trait_dir, paste0(gene, "_individuals.csv")))
  f <- candidates[file.exists(candidates)][1]
  if (is.na(f))
    stop("No individuals file for gene '", gene, "' under ", trait_dir,
         ". Confirm the gene is superior for this trait.", call. = FALSE)
  
  ind <- data.table::fread(f, data.table = FALSE)
  idc  <- intersect(c("ID", "id", "IID"), names(ind))[1]
  hapc <- intersect(c("Haplotype", "haplotype"), names(ind))[1]
  if (is.na(idc) || is.na(hapc))
    stop("Individuals file must have `ID` and `Haplotype` columns.",
         call. = FALSE)
  
  # Resolve the superior haplotype label
  if (is.null(superior_haplotype)) {
    s <- tryCatch(
      data.table::fread(file.path(trait_dir, "superior_genes_summary.csv"),
                        data.table = FALSE),
      error = function(e) NULL)
    if (!is.null(s) && "Superior_Haplotype" %in% names(s)) {
      gc <- intersect(c("Gene", "RAPDB_Gene_ID", "MSU7_Gene_ID",
                        "MSU_Gene_ID", "Gene_ID"),
                      names(s))[1]
      if (!is.na(gc)) {
        hit <- s[as.character(s[[gc]]) == gene, ]
        if (nrow(hit))
          superior_haplotype <- as.character(hit$Superior_Haplotype[1])
      }
    }
  }
  if (is.null(superior_haplotype) || is.na(superior_haplotype))
    stop("Could not determine superior_haplotype for gene '", gene,
         "'. Pass it explicitly.", call. = FALSE)
  
  # Match to the HGSI individual set
  idx <- match(hbr_clean_id(x$ids), hbr_clean_id(ind[[idc]]))
  hap <- as.character(ind[[hapc]])[idx]
  
  list(carrier            = !is.na(hap) & hap == superior_haplotype,
       haplotype          = hap,
       superior_haplotype = superior_haplotype,
       gene               = gene,
       trait              = trait)
}

#' Load one HaploBreedR run into an `hbr_run` object
#'
#' Internal worker shared by [hbr_read_results()]. Reads the pipeline
#' summary, the tested-gene table and the trait table for a single run
#' and assembles them into a classed object.
#'
#' @param x A [haplobreed_analysis()] result or a pipeline output
#'   directory.
#' @param trait Trait label; default the pipeline's `Trait`.
#' @param phenotypes Optional wide data.frame (ID, trait, PC columns).
#' @param phenotype_id_col ID column of `phenotypes`.
#' @param gate_params Optional override list for the gate settings.
#' @return An object of class `hbr_run` (a list).
#' @noRd
.hbr_read_one <- function(x, trait, phenotypes, phenotype_id_col,
                          gate_params) {
  is_dir <- is.character(x) && length(x) == 1L
  pr   <- if (is_dir) NULL else x
  pdir <- if (is_dir) x else NULL
  ps <- .hgsi_read_summary(pr, pdir)
  nm <- if (is.null(trait)) as.character(ps$Trait[1]) else trait
  ld <- .hgsi_load_trait(nm, pr, pdir, phenotypes, phenotype_id_col,
                         gate_params, verbose = FALSE)
  tt <- .hgsi_read_tested(pr, ld$out_dir, ps)
  ld$table            <- tt$table
  ld$table_source     <- tt$source
  ld$summary          <- ps
  ld$pipeline_summary <- ps
  class(ld) <- c("hbr_run", "list")
  ld
}

#' Normalise a run input to a named list of `hbr_run` objects
#'
#' Accepts a single `hbr_run`, a named list of them, or anything
#' [hbr_read_results()] takes. Used by [hbr_haplotype_blocks()] and
#' [hbr_regate()].
#'
#' @param x An `hbr_run`, a named list of them, or a loadable run input.
#' @param ... Passed to [hbr_read_results()].
#' @return A named list of `hbr_run` objects.
#' @noRd
.hbr_as_runs <- function(x, ...) {
  if (inherits(x, "hbr_run")) return(stats::setNames(list(x), x$trait))
  if (is.list(x) && length(x) &&
      all(vapply(x, inherits, logical(1), "hbr_run")))
    return(x)
  r <- hbr_read_results(x, ...)
  if (inherits(r, "hbr_run")) stats::setNames(list(r), r$trait) else r
}

#' Read finished HaploBreedR results from disk or memory
#'
#' Loads the output of [haplobreed_analysis()] — either the in-memory
#' result or its output directory — into one object that carries the
#' trait direction, the effective gate settings, the PCs the run
#' adjusted for (none if it used no PCs), the tested and superior gene
#' lists, the functional SNPs per gene, the phenotype table, and the
#' publication-style results table. [hbr_advanced_decision()],
#' [hbr_regate()] and [hbr_haplotype_blocks()] consume runs the same
#' way.
#'
#' A single pipeline directory or result returns one `hbr_run`. A
#' named character vector of directories or a named list of results
#' returns one `hbr_run` per element, classed as `hbr_run_list`.
#'
#' @param x A [haplobreed_analysis()] result, a pipeline output
#'   directory, a named character vector of directories, or a named
#'   list of results. **Alternatively**, the legacy arguments `path`
#'   and `traits` (see below) are still honoured.
#' @param trait Trait label; default the pipeline's stored `Trait`.
#' @param phenotypes Optional wide data.frame (ID column plus one
#'   column per trait, plus PC columns if the run used PCs). Needed
#'   only when the run used PCs and is read from files, because
#'   `phenotype_data_used.csv` does not carry the PC columns.
#' @param phenotype_id_col ID column of `phenotypes`. Default `"ID"`.
#' @param gate_params Optional named list overriding the stored gate
#'   settings (`alpha`, `min_d`, `min_ind`, `min_hap_freq`, `maf_min`,
#'   `trait_threshold`, `require_fdr`, `hedges`). Unspecified entries
#'   keep the values stored by the pipeline.
#' @param path Legacy alias: a single pipeline output directory.
#' @param traits Legacy alias: a named character vector of pipeline
#'   output directories, one element per trait.
#'
#' @return An object of class `hbr_run` (a list), or a named list of
#'   them classed as `hbr_run_list`. Elements of a run:
#'   \describe{
#'     \item{`trait`, `direction`, `mode`}{run identity and the trait
#'       direction (`"higher"` or `"lower"`)}
#'     \item{`gate`}{effective gate settings}
#'     \item{`phenotype`, `pcs`, `pc_cols`}{trait table and PC columns}
#'     \item{`gene_snps`}{functional SNPs per tested gene}
#'     \item{`tested`, `superior`, `superior_label`}{gene lists}
#'     \item{`table`, `table_source`}{publication-style results table
#'       and the file it came from}
#'     \item{`summary`, `pipeline_summary`}{the pipeline-summary row}
#'   }
#'
#' @seealso [hbr_regate()] to re-apply the gate with different
#'   settings, [hbr_haplotype_blocks()] to group superior haplotypes
#'   into blocks.
#'
#' @examples
#' toy <- system.file("extdata", "hgsi_toy", "trait_YLD",
#'                    package = "HaploBreedR")
#' if (nzchar(toy) && dir.exists(toy)) {
#'   run <- hbr_read_results(toy)
#'   run$trait
#'   length(run$superior)
#'   head(run$table)
#' }
#'
#' @export
hbr_read_results <- function(x = NULL, trait = NULL, phenotypes = NULL,
                             phenotype_id_col = "ID", gate_params = NULL,
                             path = NULL, traits = NULL) {
  
  # ---- Legacy aliases ---------------------------------------------------
  if (is.null(x)) {
    if (!is.null(traits)) x <- traits
    else if (!is.null(path)) x <- path
  }
  if (is.null(x))
    stop("supply a pipeline output directory, an in-memory result, or ",
         "the legacy `path`/`traits` argument", call. = FALSE)
  
  multi_chr <- is.character(x) && length(x) > 1L
  multi_lst <- is.list(x) && is.null(x$pipeline_summary) &&
    length(x) > 0L && !inherits(x, "hbr_run")
  
  if (multi_chr || multi_lst) {
    nms <- names(x)
    if (is.null(nms))
      nms <- if (multi_chr) basename(x) else NULL
    if (is.null(nms))
      stop("a list of results must be named", call. = FALSE)
    res <- stats::setNames(
      lapply(seq_along(x), function(i)
        .hbr_read_one(x[[i]], nms[i], phenotypes, phenotype_id_col,
                      gate_params)),
      nms)
    class(res) <- "hbr_run_list"
    return(res)
  }
  
  .hbr_read_one(x, trait, phenotypes, phenotype_id_col, gate_params)
}
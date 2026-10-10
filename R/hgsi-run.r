# ============================================================================
# R/hgsi-run.R
# Part 2 of the HGSI pipeline: haplotype-block calling, cross ranking, parent-
# set selection. Consumes an `hgsi_fit` object from hbr_fit_models().
#
#' Design crosses from a fitted model
#'
#' Takes the `hgsi_fit` object produced by [hbr_fit_models()], groups the
#' superior haplotypes into haplotype blocks, scores crosses by usefulness
#' (Schnell and Utz, 1975) and by block coverage, and selects a diverse
#' parent set. Returns the full `hgsi_advanced` object. Because the expensive
#' cross-validation is already done, this function reruns in minutes with a
#' different `haploblock_method` or a different `genetic_map`.
#'
#' @param fit An `hgsi_fit` object from [hbr_fit_models()].
#' @param output_dir Directory for the cross-design CSVs and the final RDS.
#'   `NULL` inherits `fit$output_dir`.
#' @param cM_per_Mb,genetic_map,bin_cM Genetic map. See
#'   [hbr_advanced_decision()].
#' @param top_crosses Number of crosses reported.
#' @param relatedness_max Maximum relatedness between parents of a cross.
#' @param mating `"DH"` or `"RIL"`.
#' @param sel_frac,h_sel Selection intensity and heritability for usefulness.
#' @param cross_rank How crosses are ordered: `"uc"`, `"blocks"` or
#'   `"balanced"`.
#' @param haploblock_method See [hbr_advanced_decision()].
#' @param haploblock_rmax,haploblock_r2,haploblock_ld_max_cM,haploblock_max_span_cM
#'   Map-based / LD block thresholds.
#' @param haploblock_flank_cM,haploblock_max_snps,haploblock_min_gamete
#'   SNP-level method settings.
#' @param run_set_selection,set_size,n_candidates,block_cM
#'   Set-selection settings.
#' @param verbose Print progress.
#'
#' @return An `hgsi_advanced` object (a superset of the `hgsi_fit`).
#'
#'@examples
#' toy <- system.file("extdata", "hgsi_toy", package = "HaploBreedR")
#' td  <- file.path(toy, "trait_YLD")
#' fit <- try(hbr_fit_models(
#'   trait_dirs   = c(YLD = td),
#'   plink_prefix = file.path(toy, "geno_toy"),
#'   output_dir   = tempdir(),
#'   cv_repeats   = 1L, cv_folds = 2L, strict = FALSE,
#'   verbose      = FALSE))
#' if (inherits(fit, "hgsi_fit")) {
#'   try(hbr_cross_design(
#'     fit,
#'     output_dir        = tempdir(),
#'     haploblock_method = "recomb",
#'     run_set_selection = FALSE,
#'     verbose           = FALSE))
#' }
#'
#' @export
hbr_cross_design <- function(
    fit,
    output_dir = NULL,
    cM_per_Mb = NULL,
    genetic_map = NULL,
    bin_cM = 1,
    top_crosses = 20L,
    relatedness_max = 0.30,
    mating = c("DH", "RIL"),
    sel_frac = 0.10,
    h_sel = 1,
    cross_rank = c("uc", "blocks", "balanced"),
    haploblock_method = c("recomb", "ld", "fourgamete", "gabriel",
                          "recomb+gabriel", "recomb+fourgamete"),
    haploblock_rmax = 0.01,
    haploblock_r2 = 0.5,
    haploblock_ld_max_cM = 3,
    haploblock_max_span_cM = 5,
    haploblock_flank_cM = 1,
    haploblock_max_snps = 400L,
    haploblock_min_gamete = 0.01,
    run_set_selection = TRUE,
    set_size = 8L,
    n_candidates = 60L,
    block_cM = 5,
    verbose = TRUE) {
  
  if (!inherits(fit, "hgsi_fit"))
    stop("`fit` must be an hgsi_fit object from hbr_fit_models().",
         call. = FALSE)
  mating      <- match.arg(mating)
  cross_rank  <- match.arg(cross_rank)
  haploblock_method <- match.arg(haploblock_method)
  
  warn <- fit$warnings
  addw <- function(msg) {
    warn <<- c(warn, msg)
    if (verbose) message("  WARNING: ", msg)
  }
  t_all <- Sys.time()
  
  # ------------------------------------------------------------------------
  # Unpack the fit
  # ------------------------------------------------------------------------
  traits      <- fit$traits
  ids         <- fit$ids
  n           <- length(ids)
  grm         <- fit$grm
  G_all       <- fit$G_all
  codes_l     <- fit$codes_l
  loaded      <- fit$loaded
  fits        <- fit$fits
  cv_l        <- fit$cv
  gebv_mat    <- fit$gebv_mat
  w_eff       <- fit$w_eff
  idx         <- fit$index
  weights     <- fit$weights
  plink_prefix <- fit$plink_prefix
  if (is.null(output_dir)) output_dir <- fit$output_dir
  ord <- order(idx, decreasing = TRUE)
  
  top_parents <- nrow(fit$parents)
  
  # ------------------------------------------------------------------------
  # Stage 1 - top parents and cross pairs
  # ------------------------------------------------------------------------
  if (verbose) message("[1/3] Cross pairs and usefulness")
  top_idx <- ord[seq_len(min(top_parents, n))]
  pairs <- t(utils::combn(seq_along(top_idx), 2))
  relv <- G_all[cbind(top_idx[pairs[, 1]], top_idx[pairs[, 2]])]
  ok_rel <- relv <= relatedness_max
  if (!any(ok_rel))
    stop("no parent pair satisfies relatedness_max = ",
         relatedness_max, call. = FALSE)
  
  bim <- .hgsi_read_bim(plink_prefix)
  gpos <- .hgsi_genetic_pos(bim, cM_per_Mb, genetic_map)
  if (!is.null(gpos) && !identical(bim$snp, grm$snps)) {
    addw("SNP IDs of the .bim and of BEDMatrix differ; genetic positions ",
         "cannot be used.")
    gpos <- NULL
  }
  crosses <- data.table::data.table()
  selected_set <- NULL
  blocks_tab <- NULL; block_set <- NULL; tgt <- NULL
  map_source <- NA_character_
  parents <- fit$parents
  
  if (!is.null(gpos)) {
    map_source <- gpos$source
    bins <- .hgsi_make_bins(bim$chr, gpos$cm, bin_cM)
    cand_n <- if (isTRUE(run_set_selection))
      max(n_candidates, top_parents) else top_parents
    cand_idx <- ord[seq_len(min(cand_n, n))]
    Tn <- length(traits)
    V <- matrix(0, n, Tn); scale_v <- numeric(Tn); excl <- vector("list", Tn)
    for (t in seq_len(Tn)) {
      f <- fits[[t]]; V[f$train, t] <- f$vir
      scale_v[t] <- f$s2R / f$kernel_scale; excl[[t]] <- f$exclude
    }
    U <- .hgsi_snp_effects(plink_prefix, grm, V, scale_v, excl)
    S <- .hgsi_parent_bin_sums(plink_prefix, grm, U, cand_idx,
                               bins$bin_of, bins$n)
    
    gene_first <- list()
    for (tr in traits) for (g in names(loaded[[tr]]$gene_snps))
      if (is.null(gene_first[[g]])) {
        ii <- stats::na.omit(match(loaded[[tr]]$gene_snps[[g]], grm$snps))
        if (length(ii)) gene_first[[g]] <- ii[1]
      }
    gene_bin <- vapply(gene_first, function(i) bins$bin_of[i], integer(1))
    gene_cm  <- vapply(gene_first, function(i) gpos$cm[i], numeric(1))
    gene_chr <- vapply(gene_first, function(i) bim$chr[i], character(1))
    for (t in seq_len(Tn)) {
      f <- fits[[t]]
      if (f$model == "two_kernel")
        S <- .hgsi_add_gene_effects(S, t, codes_l[[t]], ids[cand_idx],
                                    f$hap_info, f$hap_effects, gene_bin)
    }
    
    uc <- .hgsi_cross_uc(S, gebv_mat[cand_idx, , drop = FALSE], bins, w_eff,
                         pairs[ok_rel, , drop = FALSE], sel_frac, h_sel,
                         mating)
    
    # ------------------------------------------------------------------
    # Stage 2 - haplotype blocks
    # ------------------------------------------------------------------
    if (verbose) message("[2/3] Haplotype blocks (", haploblock_method, ")")
    sup_tabs <- lapply(stats::setNames(traits, traits), function(tr) {
      nd <- cv_l[[tr]]$naive_discovery; g <- loaded[[tr]]$superior
      data.frame(gene = g, hap = nd$best_hap[match(g, nd$gene)],
                 stringsAsFactors = FALSE)
    })
    tgt <- .hgsi_target_blocks(sup_tabs, codes_l, ids, gene_chr, gene_cm,
                               gene_first, plink_prefix, gpos$cm, bim$chr,
                               haploblock_method, haploblock_rmax,
                               haploblock_r2, haploblock_ld_max_cM,
                               haploblock_max_span_cM, haploblock_flank_cM,
                               haploblock_max_snps, haploblock_min_gamete)
    units <- NULL; block <- integer(0)
    ucarry <- matrix(FALSE, length(cand_idx), 0L)
    if (!is.null(tgt)) {
      units <- tgt$units; block <- tgt$block
      ucarry <- tgt$carry[cand_idx, , drop = FALSE]
      blocks_tab <- tgt$blocks
      blocks_tab$n_candidates_full <-
        colSums(tgt$full[cand_idx, , drop = FALSE])
      if (verbose)
        message(sprintf("    %d target haplotypes -> %d haplotype blocks",
                        nrow(units), nrow(blocks_tab)))
      parents <- cbind(parents, data.table::data.table(
        n_blocks_full    = rowSums(tgt$full[top_idx, , drop = FALSE]),
        n_blocks_partial = rowSums(tgt$partial[top_idx, , drop = FALSE]),
        n_blocks_target  = nrow(blocks_tab)))
    } else {
      addw("no superior haplotype with a map position: crosses are not ",
           "scored on haplotype blocks.")
    }
    
    pb <- .hgsi_cross_blocks(ucarry, block, cbind(uc$a, uc$b),
                             units$chr, units$cm, mating)
    
    crosses <- data.table::data.table(
      parent_a             = ids[cand_idx[uc$a]],
      parent_b             = ids[cand_idx[uc$b]],
      relatedness          = relv[ok_rel],
      mid_index            = uc$mid_index,
      sd_index             = uc$sd_index,
      uc                   = uc$uc,
      n_blocks_fixed       = pb[, "n_blocks_fixed"],
      n_blocks_segregating = pb[, "n_blocks_segregating"],
      n_blocks_split       = pb[, "n_blocks_split"],
      n_blocks_incomplete  = pb[, "n_blocks_incomplete"],
      n_blocks_covered     = pb[, "n_blocks_covered"],
      n_blocks_target      = length(unique(block)),
      p_pyramid_DH         = pb[, "p_pyramid_DH"],
      p_pyramid_F2         = pb[, "p_pyramid_F2"],
      p_pyramid_F2_hom     = pb[, "p_pyramid_F2_hom"],
      n_DH_95              = pb[, "n_DH_95"],
      n_F2_95              = .hgsi_n_for_confidence(pb[, "p_pyramid_F2"]),
      n_F2hom_95           = .hgsi_n_for_confidence(pb[, "p_pyramid_F2_hom"]))
    
    for (tr in traits) {
      mid_tr <- (gebv_mat[cand_idx[uc$a], tr] +
                   gebv_mat[cand_idx[uc$b], tr]) / 2
      data.table::set(crosses, j = paste0(tr, "_mid_gebv"), value = mid_tr)
    }
    crosses <- .hgsi_rank_crosses(crosses, cross_rank)
    crosses <- crosses[seq_len(min(top_crosses, nrow(crosses)))]
    
    # ------------------------------------------------------------------
    # Stage 3 - set selection
    # ------------------------------------------------------------------
    if (verbose) message("[3/3] Set selection")
    if (isTRUE(run_set_selection) && ncol(ucarry) > 0L &&
        length(cand_idx) > set_size) {
      bs <- .hgsi_select_set_blocks(ucarry, block, set_size,
                                    K = G_all[cand_idx, cand_idx],
                                    max_rel = relatedness_max,
                                    score = idx[cand_idx])
      block_set <- list(ids = ids[cand_idx][bs$set],
                        index = idx[cand_idx][bs$set],
                        unique_blocks = bs$unique_blocks,
                        n_blocks_covered = bs$n_blocks_covered,
                        n_blocks = bs$n_blocks,
                        coverage = bs$coverage,
                        truncation_coverage = bs$trunc_coverage,
                        mean_relatedness = bs$mean_relatedness,
                        feasible = bs$feasible,
                        uncovered = blocks_tab$block[bs$uncovered])
    }
    if (isTRUE(run_set_selection) && length(cand_idx) > set_size) {
      H <- .hgsi_block_values(S, w_eff, bins, block_cM)
      sel <- .hgsi_select_set(H, set_size, K = G_all[cand_idx, cand_idx],
                              max_rel = relatedness_max, seed = fit$settings$seed)
      selected_set <- list(set = sel$set,
                           ids = ids[cand_idx][sel$set],
                           value = sel$value,
                           trunc_value = sel$trunc_value,
                           coverage = sel$coverage,
                           mean_relatedness = sel$mean_relatedness,
                           feasible = sel$feasible,
                           block_cM = block_cM)
    }
  } else {
    addw(paste0("No genetic map available (informative .bim cM column, ",
                "`genetic_map` or `cM_per_Mb`): progeny variance, ",
                "usefulness and set selection are not computed; crosses ",
                "are ranked by the mid-parent index only."))
    mid <- vapply(seq_len(nrow(pairs)),
                  function(k) mean(idx[top_idx[pairs[k, ]]]), numeric(1))
    crosses <- data.table::data.table(
      parent_a = ids[top_idx[pairs[, 1]]],
      parent_b = ids[top_idx[pairs[, 2]]],
      relatedness = relv, mid_index = mid)[ok_rel]
    crosses <- crosses[order(-crosses$mid_index)][
      seq_len(min(top_crosses, nrow(crosses)))]
  }
  
  # ------------------------------------------------------------------------
  # Reassemble the hgsi_advanced object
  # ------------------------------------------------------------------------
  result <- list(
    traits = fit$traits, directions = fit$directions,
    ids = fit$ids, pheno_mat = fit$pheno_mat,
    gebv_mat = fit$gebv_mat, rel_mat = fit$rel_mat,
    sd_g = fit$sd_g, w_eff = fit$w_eff, index = fit$index,
    fits = fit$fits, model_used = fit$model_used,
    cv = fit$cv, cv_summary = fit$cv_summary, cv_compare = fit$cv_compare,
    concordance = fit$concordance, fidelity = fit$fidelity,
    evaluation_label = fit$evaluation_label,
    parents = parents, parents_all = fit$parents_all,
    crosses = crosses, selected_set = selected_set,
    blocks = blocks_tab, block_set = block_set,
    target = if (!is.null(tgt))
      tgt[c("units", "block", "carry", "full", "partial")] else NULL,
    conflicts = fit$conflicts,
    weights = fit$weights, plink_prefix = fit$plink_prefix,
    output_dir = output_dir,
    provenance = fit$provenance,
    warnings = warn,
    map_source = map_source,
    settings = c(fit$settings, list(
      mating = mating, sel_frac = sel_frac, h_sel = h_sel,
      cross_rank = cross_rank,
      top_crosses = top_crosses, relatedness_max = relatedness_max,
      run_set_selection = run_set_selection,
      set_size = set_size, n_candidates = n_candidates,
      block_cM = block_cM,
      haploblocks = list(method = haploblock_method,
                         rmax = haploblock_rmax, r2 = haploblock_r2,
                         ld_max_cM = haploblock_ld_max_cM,
                         max_span_cM = haploblock_max_span_cM,
                         flank_cM = haploblock_flank_cM,
                         max_snps = haploblock_max_snps,
                         min_gamete = haploblock_min_gamete))))
  class(result) <- c("hgsi_advanced", "list")
  
  if (!is.null(output_dir) && !is.na(output_dir)) {
    .hgsi_write_cross_outputs(result, output_dir, verbose)
  }
  
  .hbr_log_time(t_all, "cross design", verbose)
  invisible(result)
}


# ============================================================================
# One-shot wrapper
# ============================================================================

#' Breeding decisions from HaploBreedR output
#'
#' Convenience wrapper that runs [hbr_fit_models()] followed by
#' [hbr_cross_design()] and returns the combined `hgsi_advanced` object.
#' See the two component functions for the full argument descriptions.
#'
#' @inheritParams hbr_fit_models
#' @inheritParams hbr_cross_design
#'
#' @return An `hgsi_advanced` object.
#' @examples
#' toy <- system.file("extdata", "hgsi_toy", package = "HaploBreedR")
#' td  <- file.path(toy, "trait_YLD")
#' fit <- try(hbr_advanced_decision(
#'   trait_dirs   = c(YLD = td),
#'   plink_prefix = file.path(toy, "geno_toy"),
#'   output_dir   = tempdir(),
#'   cv_repeats   = 1L, cv_folds = 2L, strict = FALSE,
#'   run_set_selection = FALSE,
#'   save_plot    = FALSE,
#'   verbose      = FALSE))
#' @export
hbr_advanced_decision <- function(
    # shared
  pipeline_result = NULL, pipeline_results = NULL, trait_dirs = NULL,
  output_dir = NULL, phenotypes = NULL, phenotype_id_col = "ID",
  gate_params = NULL, min_concordance = 0.6,
  plink_prefix = NULL, weights = NULL,
  top_parents = 25L, n_cores = 1L, seed = 1L, verbose = TRUE,
  # fit-side
  cv_repeats = 10L, cv_folds = 5L, top_frac = 0.10, min_count = 5L,
  strata = NULL, n_strata = 3L, strict = TRUE,
  final_model = c("auto", "gblup", "two_kernel"),
  boot = 2000L,
  # cross-design side
  cM_per_Mb = NULL, genetic_map = NULL, bin_cM = 1,
  top_crosses = 20L, relatedness_max = 0.30,
  mating = c("DH", "RIL"), sel_frac = 0.10, h_sel = 1,
  cross_rank = c("uc", "blocks", "balanced"),
  haploblock_method = c("recomb", "ld", "fourgamete", "gabriel",
                        "recomb+gabriel", "recomb+fourgamete"),
  haploblock_rmax = 0.01, haploblock_r2 = 0.5,
  haploblock_ld_max_cM = 3, haploblock_max_span_cM = 5,
  haploblock_flank_cM = 1, haploblock_max_snps = 400L,
  haploblock_min_gamete = 0.01,
  run_set_selection = TRUE, set_size = 8L,
  n_candidates = 60L, block_cM = 5) {
  
  fit <- hbr_fit_models(
    pipeline_result = pipeline_result,
    pipeline_results = pipeline_results,
    trait_dirs = trait_dirs,
    output_dir = output_dir,
    phenotypes = phenotypes,
    phenotype_id_col = phenotype_id_col,
    gate_params = gate_params,
    min_concordance = min_concordance,
    plink_prefix = plink_prefix,
    weights = weights,
    cv_repeats = cv_repeats, cv_folds = cv_folds,
    top_frac = top_frac, min_count = min_count,
    strata = strata, n_strata = n_strata, strict = strict,
    final_model = match.arg(final_model),
    top_parents = top_parents,
    n_cores = n_cores, boot = boot, seed = seed, verbose = verbose)
  
  hbr_cross_design(
    fit = fit,
    output_dir = output_dir,
    cM_per_Mb = cM_per_Mb, genetic_map = genetic_map, bin_cM = bin_cM,
    top_crosses = top_crosses, relatedness_max = relatedness_max,
    mating = match.arg(mating), sel_frac = sel_frac, h_sel = h_sel,
    cross_rank = match.arg(cross_rank),
    haploblock_method = match.arg(haploblock_method),
    haploblock_rmax = haploblock_rmax, haploblock_r2 = haploblock_r2,
    haploblock_ld_max_cM = haploblock_ld_max_cM,
    haploblock_max_span_cM = haploblock_max_span_cM,
    haploblock_flank_cM = haploblock_flank_cM,
    haploblock_max_snps = haploblock_max_snps,
    haploblock_min_gamete = haploblock_min_gamete,
    run_set_selection = run_set_selection, set_size = set_size,
    n_candidates = n_candidates, block_cM = block_cM,
    verbose = verbose)
}


# ============================================================================
# Cross-design output writer
# ============================================================================

#' Write the cross-design outputs
#'
#' @param result An `hgsi_advanced` object.
#' @param dir Output directory.
#' @param verbose Print progress.
#' @noRd
.hgsi_write_cross_outputs <- function(result, dir, verbose = TRUE) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  hdir <- file.path(dir, "hgsi")
  dir.create(hdir, recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    tabs <- result[c("parents", "crosses", "blocks",
                     "cv_summary", "cv_compare", "concordance", "fidelity",
                     "provenance")]
    for (nm in names(tabs))
      if (!is.null(tabs[[nm]]) && NROW(tabs[[nm]]))
        data.table::fwrite(tabs[[nm]],
                           file.path(hdir, paste0("hgsi_", nm, ".csv")))
    
    ss <- result$selected_set
    if (!is.null(ss))
      data.table::fwrite(
        data.frame(id = ss$ids, set_value = ss$value,
                   truncation_value = ss$trunc_value,
                   coverage = ss$coverage,
                   mean_relatedness = ss$mean_relatedness,
                   feasible = ss$feasible),
        file.path(hdir, "hgsi_selected_set.csv"))
    
    bs <- result$block_set
    if (!is.null(bs))
      data.table::fwrite(
        data.frame(id = bs$ids, index = bs$index,
                   unique_blocks = bs$unique_blocks,
                   set_blocks_covered = bs$n_blocks_covered,
                   n_blocks = bs$n_blocks, coverage = bs$coverage,
                   truncation_coverage = bs$truncation_coverage,
                   mean_relatedness = bs$mean_relatedness,
                   feasible = bs$feasible,
                   uncovered_blocks = paste(bs$uncovered, collapse = ";")),
        file.path(hdir, "hgsi_block_set.csv"))
    
    # GEBV table
    tbl <- data.table::data.table(id = result$ids)
    for (tr in result$traits) {
      data.table::set(tbl, j = paste0(tr, "_predicted"),
                      value = result$fits[[tr]]$pred)
      data.table::set(tbl, j = paste0(tr, "_reliability"),
                      value = result$fits[[tr]]$reliability)
      data.table::set(tbl, j = paste0(tr, "_observed"),
                      value = result$pheno_mat[, tr])
    }
    data.table::set(tbl, j = "index", value = result$index)
    data.table::fwrite(tbl, file.path(hdir, "hgsi_gebv_predicted.csv"))
    
    # Markdown report
    md <- c(
      "# HaploBreedR advanced decision report", "",
      sprintf("Evaluation label per trait: %s",
              paste(sprintf("%s = %s", names(result$evaluation_label),
                            result$evaluation_label), collapse = "; ")),
      "",
      sprintf("Genetic map source: %s",
              if (is.na(result$map_source))
                "none (crosses ranked by mid-parent index only)"
              else result$map_source),
      sprintf("Mating system: %s | selection fraction: %s",
              result$settings$mating,
              format(result$settings$sel_frac %||% NA)),
      sprintf("Haplotype blocks: method = %s",
              result$settings$haploblocks$method),
      "")
    if (length(result$warnings))
      md <- c(md, "## Warnings", paste0("- ", result$warnings), "")
    md <- c(md,
            "## Models used for ranking",
            paste0("- ", names(result$model_used), ": ", result$model_used), "",
            "## Cross-validation (mean over folds and repeats)",
            "M0 = GBLUP; M2 = features chosen on ALL data (optimistic); ",
            "M3 = features re-discovered inside each training fold.", "",
            "```",
            utils::capture.output(print(result$cv_summary, row.names = FALSE,
                                        digits = 3)),
            "```", "",
            "## Comparisons", "",
            "```",
            utils::capture.output(print(result$cv_compare, row.names = FALSE,
                                        digits = 3)),
            "```", "",
            "## Gate replay concordance", "",
            "```",
            utils::capture.output(print(result$concordance, row.names = FALSE,
                                        digits = 3)),
            "```", "",
            "## Top 10 parents", "",
            "```",
            utils::capture.output(print(utils::head(result$parents, 10))),
            "```", "")
    if (NROW(result$crosses))
      md <- c(md, "## Top 10 crosses", "",
              "```",
              utils::capture.output(print(utils::head(result$crosses, 10))),
              "```", "")
    if (NROW(result$blocks)) {
      hb <- result$settings$haploblocks
      md <- c(md, "## Haplotype blocks", "",
              sprintf("Grouped by %s (rmax = %s, r2 >= %s within %s cM, ",
                      hb$method, hb$rmax, hb$r2, hb$ld_max_cM),
              sprintf("max span %s cM).", hb$max_span_cM), "",
              "```",
              utils::capture.output(print(utils::head(result$blocks, 20),
                                          row.names = FALSE)),
              "```", "")
    }
    writeLines(md, file.path(hdir, "advanced_report.md"))
    
    saveRDS(result, file.path(hdir, "advanced_result.rds"))
    if (length(result$warnings))
      writeLines(result$warnings, file.path(hdir, "hgsi_warnings.txt"))
    if (verbose) message("  Cross-design outputs written to ", dir)
  }, error = function(e)
    warning("could not write cross-design outputs: ",
            conditionMessage(e), call. = FALSE))
  invisible(dir)
}

#' Print an hgsi_advanced object
#'
#' Console summary of a completed [hbr_advanced_decision()] or
#' [hbr_cross_design()] run: evaluation labels, traits, models,
#' parent and cross counts, block and set summaries, genetic-map source
#' and any warnings recorded during the run.
#'
#' @param x An `hgsi_advanced` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @examples
#' x <- structure(
#'   list(evaluation_label = c(YLD = "honest_discovery_in_fold"),
#'        traits = "YLD",
#'        model_used = c(YLD = "gblup"),
#'        ids = sprintf("IND%02d", 1:5),
#'        parents = data.frame(id = sprintf("IND%02d", 1:5), index = 5:1),
#'        crosses = data.frame(parent_a = "IND01", parent_b = "IND02",
#'                             uc = 1.2),
#'        map_source = "bundled rice map",
#'        warnings = character(0),
#'        output_dir = tempdir()),
#'   class = c("hgsi_advanced", "list"))
#' print(x)
#' @export
print.hgsi_advanced <- function(x, ...) {
  cat("<hgsi_advanced>\n")
  cat(sprintf("  evaluation   : %s\n",
              paste(sprintf("%s=%s", names(x$evaluation_label),
                            x$evaluation_label), collapse = "; ")))
  cat(sprintf("  traits       : %s\n", paste(x$traits, collapse = ", ")))
  cat(sprintf("  models       : %s\n",
              paste(sprintf("%s=%s", names(x$model_used), x$model_used),
                    collapse = "; ")))
  cat(sprintf("  n individuals: %d\n", length(x$ids)))
  cat(sprintf("  parents      : %d ranked\n", nrow(x$parents)))
  cat(sprintf("  crosses      : %d\n", NROW(x$crosses)))
  if (!is.null(x$blocks))
    cat(sprintf("  blocks       : %d (%s)\n", nrow(x$blocks),
                x$settings$haploblocks$method))
  if (!is.null(x$selected_set))
    cat(sprintf("  selected set : %d parents (SA)\n",
                length(x$selected_set$set)))
  if (!is.null(x$block_set))
    cat(sprintf("  block set    : %d parents, coverage %.2f\n",
                length(x$block_set$ids), x$block_set$coverage))
  if (!is.null(x$map_source) && !is.na(x$map_source))
    cat(sprintf("  map source   : %s\n", x$map_source))
  if (length(x$warnings))
    cat(sprintf("  warnings     : %d (see $warnings)\n", length(x$warnings)))
  cat(sprintf("  output       : %s\n", x$output_dir))
  invisible(x)
}


#' Summarise an hgsi_advanced object
#'
#' Prints the console summary from [print.hgsi_advanced()] followed by the
#' cross-validation summary table, the paired comparisons, the top ten
#' ranked parents and, when available, the top ten crosses.
#'
#' @param object An `hgsi_advanced` object.
#' @param ... Ignored.
#' @return `object`, invisibly.
#' @examples
#' x <- structure(
#'   list(evaluation_label = c(YLD = "honest_discovery_in_fold"),
#'        traits = "YLD",
#'        model_used = c(YLD = "gblup"),
#'        ids = sprintf("IND%02d", 1:5),
#'        parents = data.frame(id = sprintf("IND%02d", 1:5), index = 5:1),
#'        crosses = data.frame(parent_a = "IND01", parent_b = "IND02",
#'                             uc = 1.2),
#'        cv_summary = data.frame(trait = "YLD", method = "M0", r = 0.5),
#'        cv_compare = data.frame(trait = "YLD", method = "M3",
#'                                baseline = "M0", mean_diff = 0.01),
#'        map_source = "bundled rice map",
#'        warnings = character(0),
#'        output_dir = tempdir()),
#'   class = c("hgsi_advanced", "list"))
#' summary(x)
#' @export
summary.hgsi_advanced <- function(object, ...) {
  print(object)
  cat("\nCross-validation summary:\n")
  print(object$cv_summary, row.names = FALSE, digits = 3)
  cat("\nComparisons:\n")
  print(object$cv_compare, row.names = FALSE, digits = 3)
  cat("\nTop 10 parents:\n")
  print(utils::head(object$parents, 10))
  if (NROW(object$crosses)) {
    cat("\nTop 10 crosses:\n")
    print(utils::head(object$crosses, 10))
  }
  invisible(object)
}
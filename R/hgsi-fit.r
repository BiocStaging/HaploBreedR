# ============================================================================
# R/hgsi-fit.R
# Part 1 of the HGSI pipeline: loading, GRM, haplotype codes, gate concordance,
# cross-validation, final models, index.
#
# Returns an `hgsi_fit` object that hbr_cross_design() consumes. The fit is the
# expensive part of the pipeline (CV with strict = TRUE can take hours); the
# cross design is cheap and can be rerun with different block methods.
# ============================================================================

#' Fit genomic prediction models and rank parents
#'
#' Runs every phenotype-dependent step of the HGSI pipeline: loads the
#' pipeline output, builds the genomic relationship matrix, rebuilds the
#' haplotype codes for all tested genes, replays the discovery gate inside
#' each CV fold (M0 / M2 / M3), fits the final model on all phenotyped
#' individuals, and computes the multi-trait selection index. Returns an
#' `hgsi_fit` object that [hbr_cross_design()] consumes. Use this function
#' on its own when the same fit will be reused for several cross-design
#' configurations; otherwise call [hbr_advanced_decision()] which runs both
#' parts end to end.
#'
#' @param pipeline_result A single [haplobreed_analysis()] result.
#' @param pipeline_results Named list of results; takes precedence over
#'   `pipeline_result`.
#' @param trait_dirs Named character vector of pipeline output directories
#'   for traits without an in-memory result.
#' @param output_dir Directory for `hgsi_fit.rds` and the CV CSVs. `NULL`
#'   or `NA` writes nothing.
#' @param phenotypes Optional wide data.frame (ID + one column per trait,
#'   plus PC columns if the pipeline used PCs).
#' @param phenotype_id_col ID column of `phenotypes`.
#' @param gate_params Optional gate override list, per trait or shared.
#' @param min_concordance Minimum Jaccard overlap for the replay to count
#'   as faithful.
#' @param plink_prefix PLINK prefix; default `pipeline_summary$PLINK_Prefix`.
#' @param weights Named non-negative emphasis per genetic SD.
#' @param cv_repeats,cv_folds Cross-validation settings.
#' @param top_frac Fraction selected for realised gain.
#' @param min_count Minimum training carriers per haplotype indicator.
#' @param strata,n_strata Fold strata, or number of genotype clusters.
#' @param strict Recompute GRM allele frequencies inside each training fold.
#' @param final_model `"auto"`, `"gblup"` or `"two_kernel"`.
#' @param top_parents Number of parents reported.
#' @param n_cores Workers for the GRM.
#' @param boot Bootstrap resamples for the CV comparison.
#' @param seed Random seed.
#' @param verbose Print progress.
#'
#' @return An object of class `hgsi_fit`. Feed it to [hbr_cross_design()].
#'
#' @examples
#' toy <- system.file("extdata", "test_data", package = "HaploBreedR")
#' fit <- try(hbr_fit_models(
#'   trait_dirs   = c(trait1 = toy),
#'   plink_prefix = file.path(toy, "test"),
#'   output_dir   = tempdir(),
#'   cv_repeats   = 1L, cv_folds = 2L, strict = FALSE,
#'   verbose      = FALSE))
#'
#' @export
hbr_fit_models <- function(
    pipeline_result  = NULL,
    pipeline_results = NULL,
    trait_dirs       = NULL,
    output_dir       = NULL,

    phenotypes       = NULL,
    phenotype_id_col = "ID",
    gate_params      = NULL,
    min_concordance  = 0.6,

    plink_prefix     = NULL,
    weights          = NULL,

    cv_repeats       = 10L,
    cv_folds         = 5L,
    top_frac         = 0.10,
    min_count        = 5L,
    strata           = NULL,
    n_strata         = 3L,
    strict           = TRUE,

    final_model      = c("auto", "gblup", "two_kernel"),
    top_parents      = 25L,

    n_cores          = 1L,
    boot             = 2000L,
    seed             = 1L,
    verbose          = TRUE) {

  final_model <- match.arg(final_model)

  warn <- character()
  addw <- function(msg) {
    warn <<- c(warn, msg)
    if (verbose) message("  WARNING: ", msg)
  }
  t_all <- Sys.time()

  # ------------------------------------------------------------------------
  # 1. Normalise input to a named list of pipeline results
  # ------------------------------------------------------------------------
  if (is.null(pipeline_results) && !is.null(pipeline_result)) {
    nm <- as.character(pipeline_result$pipeline_summary$Trait[1])[1]
    if (is.na(nm) || !nzchar(nm)) nm <- "trait"
    pipeline_results <- stats::setNames(list(pipeline_result), nm)
  }
  if (is.null(pipeline_results) && is.null(trait_dirs))
    stop("Supply `pipeline_result`, `pipeline_results`, or `trait_dirs`.",
         call. = FALSE)
  if (!is.null(trait_dirs) && is.null(names(trait_dirs)))
    names(trait_dirs) <- basename(trait_dirs)
  if (!is.null(pipeline_results) && is.null(names(pipeline_results)))
    stop("`pipeline_results` must be a named list.", call. = FALSE)
  traits <- unique(c(names(pipeline_results), names(trait_dirs)))

  # ------------------------------------------------------------------------
  # 2. Load each trait
  # ------------------------------------------------------------------------
  if (verbose) message("[1/6] Loading HaploBreedR output")
  per_trait_gp <- length(names(gate_params)) > 0L &&
    all(names(gate_params) %in% traits)
  gate_for <- function(tr) if (per_trait_gp) gate_params[[tr]] else gate_params
  if (!per_trait_gp && length(traits) > 1L &&
      !is.null(gate_params$trait_threshold))
    addw(paste0("gate_params$trait_threshold is in trait units but is ",
                "applied to all traits; pass gate_params as a list named ",
                "by trait instead."))
  loaded <- lapply(traits, function(tr)
    .hgsi_load_trait(trait_name = tr,
                     pipeline_result = pipeline_results[[tr]],
                     pipeline_dir = if (tr %in% names(trait_dirs))
                       trait_dirs[[tr]] else NULL,
                     phenotypes = phenotypes,
                     phenotype_id_col = phenotype_id_col,
                     gate_params_user = gate_for(tr),
                     verbose = verbose))
  names(loaded) <- traits
  for (tr in traits) for (msg in loaded[[tr]]$notes) addw(msg)

  directions <- vapply(loaded, function(l) l$direction, character(1))
  if (is.null(weights))
    weights <- stats::setNames(rep(1, length(traits)), traits)
  if (is.null(names(weights)) && length(weights) == length(traits))
    names(weights) <- traits
  weights <- weights[traits]
  if (anyNA(weights) || any(weights < 0))
    stop("`weights` must be non-negative and named by trait (",
         paste(traits, collapse = ", "),
         "); the sign comes from the pipeline direction.", call. = FALSE)

  if (is.null(plink_prefix)) plink_prefix <- loaded[[1]]$plink_prefix
  if (is.null(plink_prefix) || is.na(plink_prefix))
    stop("Supply `plink_prefix`.", call. = FALSE)
  if (length(unique(vapply(loaded, function(l) l$plink_prefix,
                           character(1)))) > 1L)
    addw("traits were run on different PLINK prefixes; the first one is ",
         "used for all.")
  if (is.null(output_dir)) output_dir <- loaded[[1]]$out_dir

  # ------------------------------------------------------------------------
  # 3. GRM (streamed via BEDMatrix)
  # ------------------------------------------------------------------------
  if (verbose) message("[2/6] Genomic relationship matrix (streamed)")
  fam <- data.table::fread(paste0(plink_prefix, ".fam"),
                           header = FALSE, data.table = FALSE)
  ids <- hbr_clean_id(fam$V2)
  grm <- .hgsi_build_grm(plink_prefix, n_cores = n_cores, verbose = verbose)
  if (length(ids) != grm$n) stop("the .fam and .bed sizes disagree")
  n <- length(ids)
  G_all <- .scale_diag1(grm$A / grm$d)
  X_all <- matrix(1, n, 1L)

  # ------------------------------------------------------------------------
  # 4. Phenotype matrix and PCs
  # ------------------------------------------------------------------------
  pheno_mat <- matrix(vapply(traits, function(tr) {
    ph <- loaded[[tr]]$phenotype; y <- rep(NA_real_, n)
    j <- match(ids, ph$ID); y[!is.na(j)] <- ph$y[j[!is.na(j)]]; y
  }, numeric(n)), n, length(traits), dimnames = list(ids, traits))
  pcs_l <- lapply(loaded, function(l) {
    if (is.null(l$pcs)) return(NULL)
    as.matrix(l$pcs[match(ids, l$pcs$ID), l$pc_cols, drop = FALSE])
  })
  for (tr in traits) {
    nobs <- sum(!is.na(pheno_mat[, tr]))
    if (verbose)
      message(sprintf("    %s: %d / %d individuals phenotyped", tr, nobs, n))
    if (nobs < 40L)
      stop("too few phenotyped individuals for trait ", tr, call. = FALSE)
    if (nobs < n)
      addw(sprintf("%s: %d individuals lack a phenotype; they are ",
                   tr, n - nobs),
           "predicted, not used for training.")
  }

  strata_v <- if (is.null(strata)) {
    .hgsi_auto_strata(G_all, n_strata, seed)
  } else {
    sv <- if (!is.null(names(strata))) strata[ids] else strata
    if (length(sv) != n)
      stop("`strata` must have one value per individual", call. = FALSE)
    sv
  }

  # ------------------------------------------------------------------------
  # 5. Haplotype codes and gate concordance
  # ------------------------------------------------------------------------
  if (verbose)
    message("[3/6] Haplotype codes for ALL tested genes, fidelity and gate concordance")
  codes_l <- list(); conc_rows <- list(); fid_rows <- list()
  label <- character()
  for (tr in traits) {
    ld <- loaded[[tr]]; y <- pheno_mat[, tr]
    codes <- .hgsi_build_codes(plink_prefix, ld$gene_snps, ids, grm$snps,
                                maf_min = ld$gate$maf_min, verbose = verbose)
    codes <- .hgsi_relabel_codes(codes, ids, ld)
    reference <- stats::setNames(colnames(codes) %in% ld$superior,
                                 colnames(codes))
    conc <- .hgsi_gate_concordance(codes, y, ld$gate, reference, pcs_l[[tr]])
    fid <- .hgsi_check_fidelity(codes, ids, ld, conc$best_hap)
    fid_ok <- !nrow(fid) || (min(fid$agreement) >= 0.95 &&
      all(is.na(fid$superior_overlap) | fid$superior_overlap >= 0.9))
    conc_ok <- is.na(conc$jaccard) || conc$jaccard >= min_concordance
    label[tr] <- .hgsi_eval_label(ld$is_gwas,
                                   ncol(codes) >= 0.8 * length(ld$tested),
                                   conc_ok, fid_ok)
    if (nrow(fid)) fid_rows[[tr]] <- cbind(trait = tr, fid)
    conc_rows[[tr]] <- data.frame(
      trait = tr, direction = ld$direction, agreement = conc$agreement,
      kappa = conc$kappa, jaccard = conc$jaccard, n_genes = conc$n_genes,
      n_superior_pipeline = sum(reference), stringsAsFactors = FALSE)
    codes_l[[tr]] <- codes
  }

  # ------------------------------------------------------------------------
  # 6. Cross-validation
  # ------------------------------------------------------------------------
  if (verbose)
    message("[4/6] Cross-validation with discovery replay (M0 / M2 / M3)")
  cv_l <- list(); cmp_l <- list()
  for (tr in traits) {
    t0 <- Sys.time()
    cv <- .hgsi_cv(pheno_mat[, tr], ids, grm, plink_prefix, codes_l[[tr]],
                   loaded[[tr]]$gene_snps, loaded[[tr]]$gate,
                   pcs = pcs_l[[tr]], strata = strata_v,
                   repeats = cv_repeats, k = cv_folds, top_frac = top_frac,
                   min_count = min_count, seed = seed, strict = strict,
                   verbose = FALSE)
    cv_l[[tr]] <- cv
    cmp <- .hgsi_compare_all(cv, boot = boot, top_frac = top_frac, seed = seed)
    if (!is.null(cmp)) cmp_l[[tr]] <- cbind(trait = tr, cmp)
    .hbr_log_time(t0, paste0("CV ", tr), verbose)
  }
  cv_compare <- do.call(rbind, cmp_l)

  # ------------------------------------------------------------------------
  # 7. Final models
  # ------------------------------------------------------------------------
  if (verbose) message("[5/6] Final models (all phenotyped, predicting everyone)")
  fits <- list(); model_used <- character()
  for (tr in traits) {
    y <- pheno_mat[, tr]; obs <- which(!is.na(y)); codes <- codes_l[[tr]]
    cmp <- cmp_l[[tr]]
    row <- if (is.null(cmp)) NULL else
      cmp[cmp$method == "M3" & cmp$baseline == "M0", , drop = FALSE]
    use_tk <- switch(final_model,
      gblup      = FALSE,
      two_kernel = TRUE,
      auto = !is.null(row) && nrow(row) == 1L && is.finite(row$p) &&
        row$mean_diff > 0 && row$p < 0.05 && row$boot_lo > 0)
    feat_genes <- intersect(loaded[[tr]]$superior, colnames(codes))
    fit <- NULL
    if (use_tk && length(feat_genes)) {
      des <- .hgsi_hap_design(codes, feat_genes, obs, min_count)
      if (!is.null(des$X)) {
        cols <- sort(unique(stats::na.omit(match(
          unlist(loaded[[tr]]$gene_snps[feat_genes]), grm$snps))))
        feat <- .hgsi_feature_A(plink_prefix, grm$p, grm$keep, cols)
        Gr <- .scale_diag1(.hgsi_remainder(grm, feat))
        fit <- .hgsi_fit_all(y, X_all, Gr, des$X)
        fit$model <- "two_kernel"
        fit$kernel_scale <- mean(diag(grm$A - feat$A))
        fit$exclude <- feat$cols; fit$hap_info <- des$info
      }
    }
    if (is.null(fit)) {
      fit <- .hgsi_fit_all(y, X_all, G_all)
      fit$model <- "gblup"; fit$kernel_scale <- mean(diag(grm$A))
      fit$exclude <- integer(0); fit$hap_info <- NULL
    }
    fits[[tr]] <- fit; model_used[tr] <- fit$model
    if (verbose)
      message(sprintf("    %s: %s (w=%.2f)", tr, fit$model,
                      ifelse(is.na(fit$w), 0, fit$w)))
  }

  # ------------------------------------------------------------------------
  # 8. Index
  # ------------------------------------------------------------------------
  if (verbose)
    message("[6/6] Index (BLUP GEBVs; weights per genetic SD, signed by pipeline direction)")
  gebv_mat <- matrix(vapply(traits, function(tr) fits[[tr]]$g, numeric(n)),
                     n, length(traits), dimnames = list(ids, traits))
  rel_mat  <- matrix(vapply(traits, function(tr) fits[[tr]]$reliability,
                            numeric(n)),
                     n, length(traits), dimnames = list(ids, traits))
  sd_g  <- vapply(traits, function(tr) sqrt(fits[[tr]]$s2g), numeric(1))
  w_eff <- .hgsi_index_weights(weights, directions, sd_g)
  idx   <- as.numeric(gebv_mat %*% w_eff)
  ord   <- order(idx, decreasing = TRUE)

  parents_all <- data.table::data.table(id = ids[ord], index = idx[ord],
                                        rank = seq_len(n))
  extra <- do.call(cbind, lapply(traits, function(tr) {
    m <- data.frame(gebv_mat[ord, tr], rel_mat[ord, tr])
    names(m) <- paste0(tr, c("_gebv", "_reliability")); m
  }))
  parents_all <- cbind(parents_all, extra)
  parents <- parents_all[seq_len(min(top_parents, n))]

  # ------------------------------------------------------------------------
  # 9. Cross-trait conflicts (cheap, belongs with the fit)
  # ------------------------------------------------------------------------
  conflicts <- if (length(traits) > 1L) {
    .hgsi_conflicts(lapply(stats::setNames(traits, traits), function(tr) {
      nd <- cv_l[[tr]]$naive_discovery; g <- loaded[[tr]]$superior
      data.frame(gene = g, hap = nd$best_hap[match(g, nd$gene)],
                 stringsAsFactors = FALSE)
    }))
  } else data.table::data.table()

  # ------------------------------------------------------------------------
  # 10. Assemble the hgsi_fit object
  # ------------------------------------------------------------------------
  result <- list(
    traits = traits, directions = directions,
    ids = ids, pheno_mat = pheno_mat,
    gebv_mat = gebv_mat, rel_mat = rel_mat,
    sd_g = sd_g, w_eff = w_eff, index = idx,
    grm = grm, G_all = G_all, X_all = X_all,
    codes_l = codes_l, loaded = loaded,
    pcs_l = pcs_l, strata_v = strata_v,
    fits = fits, model_used = model_used,
    cv = cv_l, cv_summary = .hgsi_cv_summary(cv_l), cv_compare = cv_compare,
    concordance = do.call(rbind, conc_rows),
    fidelity = if (length(fid_rows)) do.call(rbind, fid_rows) else NULL,
    evaluation_label = label,
    parents = parents, parents_all = parents_all,
    conflicts = conflicts,
    weights = weights, plink_prefix = plink_prefix,
    output_dir = output_dir,
    provenance = do.call(rbind, lapply(traits, function(tr)
      data.frame(trait = tr, direction = directions[[tr]],
                 pcs = if (length(loaded[[tr]]$pc_cols))
                   paste(loaded[[tr]]$pc_cols, collapse = ",") else "none",
                 phenotypes = loaded[[tr]]$provenance$phenotypes,
                 tested_genes = loaded[[tr]]$provenance$tested_genes,
                 functional_snps = loaded[[tr]]$provenance$functional_snps,
                 stringsAsFactors = FALSE))),
    warnings = warn,
    settings = list(
      final_model = final_model, strict = strict, seed = seed,
      cv = list(repeats = cv_repeats, folds = cv_folds,
                top_frac = top_frac, min_count = min_count, boot = boot),
      strata = list(n_strata = n_strata),
      gate = lapply(loaded, function(l) l$gate)))
  class(result) <- c("hgsi_fit", "list")

  # ------------------------------------------------------------------------
  # 11. Write outputs (fit side)
  # ------------------------------------------------------------------------
  if (!is.null(output_dir) && !is.na(output_dir)) {
    .hgsi_write_fit_outputs(result, output_dir, verbose)
  }

  .hbr_log_time(t_all, "fit models", verbose)
  invisible(result)
}


#' Write the fit-side outputs
#'
#' @param fit An `hgsi_fit` object.
#' @param dir Output directory.
#' @param verbose Print progress.
#' @noRd
.hgsi_write_fit_outputs <- function(fit, dir, verbose = TRUE) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  hdir <- file.path(dir, "hgsi")
  dir.create(hdir, recursive = TRUE, showWarnings = FALSE)
  tryCatch({
    tabs <- fit[c("parents", "cv_summary", "cv_compare",
                  "concordance", "fidelity", "provenance")]
    for (nm in names(tabs))
      if (!is.null(tabs[[nm]]) && NROW(tabs[[nm]]))
        data.table::fwrite(tabs[[nm]],
                           file.path(hdir, paste0("hgsi_", nm, ".csv")))

    tbl <- data.table::data.table(id = fit$ids)
    for (tr in fit$traits) {
      data.table::set(tbl, j = paste0(tr, "_predicted"),
                      value = fit$fits[[tr]]$pred)
      data.table::set(tbl, j = paste0(tr, "_reliability"),
                      value = fit$fits[[tr]]$reliability)
      data.table::set(tbl, j = paste0(tr, "_observed"),
                      value = fit$pheno_mat[, tr])
    }
    data.table::set(tbl, j = "index", value = fit$index)
    data.table::fwrite(tbl, file.path(hdir, "hgsi_gebv_predicted.csv"))

    saveRDS(fit, file.path(dir, "hgsi_fit.rds"))
    if (verbose) message("  Fit outputs written to ", dir)
  }, error = function(e)
    warning("could not write fit outputs: ", conditionMessage(e),
            call. = FALSE))
  invisible(dir)
}


#' Print an hgsi_fit object
#'
#' @param x An `hgsi_fit` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' #' Print an hgsi_fit object
#'
#' @param x An `hgsi_fit` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @examples
#' fit <- structure(
#'   list(traits = "YLD",
#'        evaluation_label = c(YLD = "honest_discovery_in_fold"),
#'        model_used = c(YLD = "gblup"),
#'        ids = sprintf("IND%02d", 1:5),
#'        parents = data.frame(id = sprintf("IND%02d", 1:5), index = 5:1),
#'        warnings = character(0),
#'        output_dir = tempdir()),
#'   class = c("hgsi_fit", "list"))
#' print(fit)
#' @export
print.hgsi_fit <- function(x, ...) {
  cat("<hgsi_fit>\n")
  cat(sprintf("  traits       : %s\n", paste(x$traits, collapse = ", ")))
  cat(sprintf("  evaluation   : %s\n",
              paste(sprintf("%s=%s", names(x$evaluation_label),
                            x$evaluation_label), collapse = "; ")))
  cat(sprintf("  models       : %s\n",
              paste(sprintf("%s=%s", names(x$model_used), x$model_used),
                    collapse = "; ")))
  cat(sprintf("  n individuals: %d\n", length(x$ids)))
  cat(sprintf("  parents      : %d ranked\n", nrow(x$parents)))
  if (length(x$warnings))
    cat(sprintf("  warnings     : %d (see $warnings)\n", length(x$warnings)))
  cat(sprintf("  fit file     : %s/hgsi_fit.rds\n", x$output_dir))
  invisible(x)
}

#' Summarise an hgsi_fit object
#'
#' @param object An `hgsi_fit` object.
#' @param ... Ignored.
#' @return `object`, invisibly.
#' #' Summarise an hgsi_fit object
#'
#' @param object An `hgsi_fit` object.
#' @param ... Ignored.
#' @return `object`, invisibly.
#' @examples
#' fit <- structure(
#'   list(traits = "YLD",
#'        evaluation_label = c(YLD = "honest_discovery_in_fold"),
#'        model_used = c(YLD = "gblup"),
#'        ids = sprintf("IND%02d", 1:5),
#'        parents = data.frame(id = sprintf("IND%02d", 1:5), index = 5:1),
#'        cv_summary = data.frame(trait = "YLD", method = "M0", r = 0.5),
#'        cv_compare = data.frame(trait = "YLD", method = "M3",
#'                                baseline = "M0", mean_diff = 0.01),
#'        warnings = character(0),
#'        output_dir = tempdir()),
#'   class = c("hgsi_fit", "list"))
#' summary(fit)
#' @export
summary.hgsi_fit <- function(object, ...) {
  print(object)
  cat("\nCross-validation summary:\n")
  print(object$cv_summary, row.names = FALSE, digits = 3)
  cat("\nComparisons:\n")
  print(object$cv_compare, row.names = FALSE, digits = 3)
  cat("\nTop 10 parents:\n")
  print(utils::head(object$parents, 10))
  invisible(object)
}
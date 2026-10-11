# ============================================================================
# R/association.r
# Haplotype-trait association testing (Welch t-test, ANOVA, permutation,
# bootstrap) and the four-gate superior-haplotype call.
#
# ============================================================================


# ----------------------------------------------------------------------------
# Permutation test
# ----------------------------------------------------------------------------

#' Permutation test for haplotype-group differences
#'
#' Permutes the haplotype labels `n_perm` times and recomputes the
#' one-way ANOVA F-statistic on each permutation. The reported p-value
#' uses the standard (b + 1) / (m + 1) correction, where b is the number
#' of permutation statistics at least as large as the observed one and
#' m is the number of valid permutations.
#'
#' @param y Numeric response vector.
#' @param g Grouping factor with at least two levels.
#' @param n_perm Number of label permutations.
#'
#' @return A list with `p_value`, `observed_stat`, `statistic` (always
#'   `"F"`), and `n_valid`.
#'
#' @noRd
hbr_permutation_test <- function(y, g, n_perm = 1000) {
  
  n_perm <- as.integer(n_perm)
  if (is.na(n_perm) || n_perm < 1L)
    stop("`n_perm` must be a positive integer.", call. = FALSE)
  
  fstat <- function(yy, gg) tryCatch(
    summary(stats::aov(yy ~ gg))[[1]]$`F value`[1],
    error = function(e) NA_real_)
  
  obs <- fstat(y, g)
  if (is.na(obs))
    return(list(p_value = NA_real_, observed_stat = NA_real_,
                statistic = "F", n_valid = 0L))
  
  perm_stats <- vapply(seq_len(n_perm),
                       function(i) fstat(y, sample(g)),
                       numeric(1))
  
  n_valid <- sum(!is.na(perm_stats))
  if (n_valid < 1L)
    return(list(p_value = NA_real_, observed_stat = obs,
                statistic = "F", n_valid = 0L))
  
  n_exceed <- sum(perm_stats >= obs, na.rm = TRUE)
  p_value  <- (n_exceed + 1) / (n_valid + 1)
  
  list(p_value       = p_value,
       observed_stat = obs,
       statistic     = "F",
       n_valid       = n_valid)
}


# ----------------------------------------------------------------------------
# Bootstrap CI
# ----------------------------------------------------------------------------

#' Bootstrap confidence interval for a sample mean
#'
#' Computes a percentile bootstrap confidence interval for the mean of
#' a numeric vector using [boot::boot()]. If the vector has fewer than
#' 5 non-missing values, or the bootstrap fails, returns `NA` bounds and
#' a normal-approximation standard error.
#'
#' @param x Numeric vector.
#' @param n_boot Number of bootstrap resamples.
#'
#' @return A list with `mean`, `median`, `lower`, `upper`, and `se`.
#'
#' @noRd
hbr_bootstrap_ci <- function(x, n_boot = 1000) {
  x <- x[!is.na(x)]
  if (length(x) < 5)
    return(list(mean = mean(x), median = stats::median(x),
                lower = NA, upper = NA, se = NA))
  bm <- function(d, i) mean(d[i])
  tryCatch({
    br <- boot::boot(x, bm, R = n_boot)
    ci <- boot::boot.ci(br, type = "perc")
    list(mean   = mean(x),
         median = stats::median(x),
         lower  = ci$percent[4],
         upper  = ci$percent[5],
         se     = stats::sd(br$t))
  }, error = function(e)
    list(mean = mean(x), median = stats::median(x),
         lower = NA, upper = NA,
         se = stats::sd(x) / sqrt(length(x))))
}


# ----------------------------------------------------------------------------
# Association test
# ----------------------------------------------------------------------------

#' Test haplotype-trait association
#'
#' For two haplotypes, Welch's t-test. For three or more, one-way ANOVA
#' with Duncan (default) and/or Tukey post-hoc. Optional principal-
#' component correction regresses the trait on the PCs first and tests
#' the residuals. Permutation and bootstrap resampling are always run.
#'
#' @param hap_data Output of [hbr_define_haplotypes()].
#' @param pheno_data Phenotype data frame with an ID column and the trait.
#' @param trait_col Trait column name.
#' @param pc_cols Principal-component columns. `NULL` or empty skips PC
#'   correction.
#' @param posthoc `"duncan"`, `"tukey"`, or both.
#' @param n_perm Permutation iterations.
#' @param n_boot Bootstrap iterations.
#'
#' @return A list with `test`, `p_value`, `permutation_p`, `groups`,
#'   `data`, `omega_squared`, `cv_percent`, plus the fitted model and
#'   per-gate diagnostics. Returns `NULL` if fewer than 10 usable
#'   observations remain.
#'
#' @noRd
hbr_test_association <- function(
    hap_data,
    pheno_data,
    trait_col,
    pc_cols  = c("PC1", "PC2", "PC3", "PC4", "PC5"),
    posthoc  = c("duncan", "tukey"),
    n_perm   = 1000,
    n_boot   = 1000) {
  
  posthoc <- tolower(trimws(as.character(posthoc)))
  posthoc <- unique(posthoc[!is.na(posthoc) & nzchar(posthoc)])
  if (length(posthoc) == 0L || !all(posthoc %in% c("duncan", "tukey")))
    stop("posthoc must be 'duncan', 'tukey', or both.", call. = FALSE)
  run_duncan <- "duncan" %in% posthoc
  run_tukey  <- "tukey"  %in% posthoc
  
  dat <- merge(pheno_data,
               hap_data$haplotypes[, .(ID, Haplotype)],
               by = "ID")
  dat <- dat[!is.na(get(trait_col))]
  if (length(intersect(pheno_data$ID, hap_data$haplotypes$ID)) < 20 ||
      nrow(dat) < 10)
    return(NULL)
  
  n_hap <- data.table::uniqueN(dat$Haplotype)
  if (n_hap < 2L) return(NULL)
  
  avail_pcs <- intersect(pc_cols, names(dat))
  use_pcs   <- length(avail_pcs) > 0
  test_col  <- trait_col
  
  if (use_pcs) {
    pcm <- tryCatch(
      stats::lm(stats::as.formula(paste(trait_col, "~",
                                        paste(avail_pcs, collapse = "+"))),
                data = dat),
      error = function(e) NULL)
    dat <- data.table::copy(dat)
    dat[["__r__"]] <- if (!is.null(pcm)) stats::residuals(pcm)
    else dat[[trait_col]]
    test_col <- "__r__"
  }
  
  perm <- tryCatch(
    hbr_permutation_test(y = dat[[test_col]], g = dat$Haplotype,
                         n_perm = n_perm),
    error = function(e)
      list(p_value = NA_real_, observed_stat = NA_real_,
           statistic = "F", n_valid = 0L))
  
  ci_all <- function() {
    lapply(stats::setNames(nm = unique(dat$Haplotype)), function(h) {
      hbr_bootstrap_ci(dat[Haplotype == h, get(trait_col)], n_boot = n_boot)
    })
  }
  
  build_groups <- function(lr, ci_l) {
    groups <- merge(
      lr,
      dat[, .(Mean   = mean(get(trait_col), na.rm = TRUE),
              Median = stats::median(get(trait_col), na.rm = TRUE),
              SE     = stats::sd(get(trait_col), na.rm = TRUE) / sqrt(.N)),
          by = Haplotype],
      by = "Haplotype")
    groups$CI_Lower <- vapply(groups$Haplotype,
                              function(h) ci_l[[as.character(h)]]$lower,
                              numeric(1))
    groups$CI_Upper <- vapply(groups$Haplotype,
                              function(h) ci_l[[as.character(h)]]$upper,
                              numeric(1))
    groups <- merge(groups, dat[, .(N = .N), by = Haplotype],
                    by = "Haplotype")
    groups <- groups[order(-Mean)]
    groups$NT_Sequence <- hap_data$haplotype_map$NT_Sequence[
      match(groups$Haplotype, hap_data$haplotype_map$Haplotype)]
    groups$AA_Sequence <- hap_data$haplotype_map$AA_Sequence[
      match(groups$Haplotype, hap_data$haplotype_map$Haplotype)]
    groups
  }
  
  base_return <- function(tname, p_val, groups, omega = NA, eff = NA, ...) {
    list(test                   = tname,
         p_value                = p_val,
         permutation_p          = perm$p_value,
         permutation_statistic  = perm$observed_stat,
         permutation_n_valid    = perm$n_valid,
         groups                 = groups,
         data                   = dat,
         omega_squared          = omega,
         effect_size            = eff,
         assumptions            = NULL,
         pc_corrected           = use_pcs,
         available_pcs          = avail_pcs,
         test_col               = test_col,
         posthoc                = posthoc,
         ...)
  }
  
  # ------------------------------------------------------------------------
  # Two haplotypes: t-test
  # ------------------------------------------------------------------------
  if (n_hap == 2L) {
    grps <- sort(unique(as.character(dat$Haplotype)))
    x1 <- dat[Haplotype == grps[1], get(test_col)]
    x2 <- dat[Haplotype == grps[2], get(test_col)]
    
    tt <- tryCatch(stats::t.test(x1, x2, var.equal = TRUE),
                   error = function(e) NULL)
    if (is.null(tt)) return(NULL)
    
    p_val <- tt$p.value
    means <- c(mean(x1, na.rm = TRUE), mean(x2, na.rm = TRUE))
    names(means) <- grps
    ord <- order(-means)
    letters_vec <- c("a", "a")
    if (is.finite(p_val) && p_val < 0.05) letters_vec[ord[2]] <- "b"
    
    lr <- data.table::data.table(
      Haplotype          = grps,
      Duncan_Group       = letters_vec,
      Tukey_Group        = letters_vec,
      Significance_Group = letters_vec,
      Letter             = letters_vec)
    
    groups <- build_groups(lr, ci_all())
    tname  <- paste0("t-test", if (use_pcs) "(PC-adj)" else "")
    
    if (!use_pcs) {
      s1 <- stats::sd(x1, na.rm = TRUE); n1 <- length(x1)
      s2 <- stats::sd(x2, na.rm = TRUE); n2 <- length(x2)
      pooled_sd  <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
      grand_mean <- mean(c(x1, x2), na.rm = TRUE)
      cv_pc <- if (is.finite(grand_mean) && grand_mean != 0)
        pooled_sd / grand_mean * 100 else NA_real_
    } else cv_pc <- NA_real_
    
    return(base_return(
      tname, p_val, groups,
      omega             = NA_real_,
      eff               = NA_real_,
      ttest             = tt,
      cv_percent        = cv_pc,
      ms_within         = NA_real_,
      anova_significant = is.finite(p_val) && p_val < 0.05))
  }
  
  # ------------------------------------------------------------------------
  # Three or more: one-way ANOVA + post-hoc
  # ------------------------------------------------------------------------
  aov_formula <- stats::as.formula(paste0("`", test_col, "` ~ Haplotype"))
  aov_m <- tryCatch(stats::aov(aov_formula, data = dat),
                    error = function(e) NULL)
  if (is.null(aov_m)) return(NULL)
  
  ss      <- summary(aov_m)[[1]]
  anova_p <- as.numeric(ss$`Pr(>F)`[1])
  ss_b    <- as.numeric(ss$`Sum Sq`[1])
  ss_t    <- sum(as.numeric(ss$`Sum Sq`), na.rm = TRUE)
  df_b    <- as.numeric(ss$Df[1])
  ms_w    <- as.numeric(ss$`Mean Sq`[2])
  omg     <- max(0, (ss_b - df_b * ms_w) / (ss_t + ms_w))
  
  cv_pc <- if (!use_pcs)
    sqrt(ms_w) / mean(dat[[trait_col]], na.rm = TRUE) * 100
  else NA_real_
  
  gh_lvls <- unique(as.character(dat$Haplotype))
  lr <- data.table::data.table(
    Haplotype          = gh_lvls,
    Duncan_Group       = "a",
    Tukey_Group        = "a",
    Significance_Group = "a")
  
  if (is.finite(anova_p) && anova_p < 0.05) {
    
    if (run_duncan && requireNamespace("agricolae", quietly = TRUE)) {
      duncan_result <- tryCatch({
        dr <- agricolae::duncan.test(aov_m, "Haplotype", alpha = 0.05)
        mn <- data.table::as.data.table(dr$groups, keep.rownames = TRUE)
        data.table::setnames(mn,
                             c("Haplotype", "Mean_Duncan", "Duncan_Group"))
        mn[, .(Haplotype, Duncan_Group = as.character(Duncan_Group))]
      }, error = function(e) NULL)
      
      if (!is.null(duncan_result)) {
        lr <- merge(lr, duncan_result, by = "Haplotype",
                    all.x = TRUE, suffixes = c("", ".new"))
        lr[, Duncan_Group := ifelse(is.na(Duncan_Group.new),
                                    Duncan_Group, Duncan_Group.new)]
        lr[, Duncan_Group.new := NULL]
      }
    }
    
    if (run_tukey &&
        requireNamespace("emmeans", quietly = TRUE) &&
        requireNamespace("multcompView", quietly = TRUE)) {
      tukey_result <- tryCatch({
        emm <- emmeans::emmeans(aov_m, ~Haplotype)
        pw  <- as.data.frame(emmeans::contrast(emm, method = "pairwise",
                                               adjust = "tukey"))
        cm  <- matrix(1, nrow = length(gh_lvls), ncol = length(gh_lvls),
                      dimnames = list(gh_lvls, gh_lvls))
        for (k in seq_len(nrow(pw))) {
          pr <- trimws(strsplit(as.character(pw$contrast[k]),
                                " - ", fixed = TRUE)[[1]])
          if (length(pr) == 2L && all(pr %in% gh_lvls)) {
            cm[pr[1], pr[2]] <- pw$p.value[k]
            cm[pr[2], pr[1]] <- pw$p.value[k]
          }
        }
        cm[is.na(cm)] <- 1
        cld_r <- multcompView::multcompLetters(cm, threshold = 0.05)
        data.table::data.table(
          Haplotype   = names(cld_r$Letters),
          Tukey_Group = as.character(cld_r$Letters))
      }, error = function(e) NULL)
      
      if (!is.null(tukey_result)) {
        lr <- merge(lr, tukey_result, by = "Haplotype",
                    all.x = TRUE, suffixes = c("", ".new"))
        lr[, Tukey_Group := ifelse(is.na(Tukey_Group.new),
                                   Tukey_Group, Tukey_Group.new)]
        lr[, Tukey_Group.new := NULL]
      }
    }
  }
  
  if (run_duncan)     lr[, Significance_Group := Duncan_Group]
  else if (run_tukey) lr[, Significance_Group := Tukey_Group]
  lr[, Letter := Significance_Group]
  
  groups <- build_groups(lr, ci_all())
  tname  <- paste0("ANOVA+", paste(posthoc, collapse = "+"),
                   if (use_pcs) "(PC-adj)" else "")
  
  base_return(
    tname, anova_p, groups,
    omega             = omg,
    eff               = omg,
    anova             = aov_m,
    cv_percent        = cv_pc,
    ms_within         = ms_w,
    anova_significant = is.finite(anova_p) && anova_p < 0.05)
}


# ----------------------------------------------------------------------------
# Superior-haplotype gate
# ----------------------------------------------------------------------------

#' Identify the superior haplotype at a gene
#'
#' Applies four gates to the association result of a single gene:
#' (A) nominal ANOVA/Welch p < `alpha`, (B) permutation p < `alpha`,
#' (C) Cohen's d >= `min_cohens_d`, and (D) an optional trait threshold.
#' The haplotype with the most extreme mean in the requested `direction`
#' is chosen as the candidate; ties are broken by group size.
#'
#' @param test_res Output of [hbr_test_association()].
#' @param direction `"higher"` or `"lower"`.
#' @param trait_col Trait column name.
#' @param min_cohens_d Minimum effect size for Gate C.
#' @param alpha Significance threshold for Gates A and B.
#' @param trait_threshold Optional trait threshold for Gate D.
#' @param n_panel Number of individuals in the phenotyped panel.
#' @param small_sample_correction If `TRUE` (default) the effect size is
#'   Hedges' g (Cohen's d times the small-sample factor J = 1 - 3/(4df - 1)),
#'   which removes the upward bias of d in small haplotype groups. The output
#'   columns keep the names `cohens_d` / `Cohens_d`. Set `FALSE` to reproduce
#'   runs made before this correction.
#'
#' @return A list with the superior haplotype and per-gate diagnostics:
#'   `haplotype`, `nt_sequence`, `aa_sequence`, `mean`, `median`, `se`,
#'   `ci_lower`, `ci_upper`, `n`, `freq_analysed_pct`, `freq_panel_pct`,
#'   `letter`, `second_haplotype`, `improvement_pct`, `cohens_d`,
#'   `diff_mean`, `diff_ci_lower`, `diff_ci_upper`, `is_superior`,
#'   `gate_a_sig`, `gate_b_perm_sig`, `gate_c_effect`,
#'   `gate_d_threshold`, `superior_label`.
#'
#' @noRd
hbr_identify_superior <- function(test_res, direction, trait_col,
                                  min_cohens_d = 0.3, alpha = 0.05,
                                  trait_threshold = NULL, n_panel = NULL,
                                  small_sample_correction = TRUE) {
  
  grp <- data.table::as.data.table(test_res$groups)
  
  extreme_mean <- if (direction == "higher") max(grp$Mean) else min(grp$Mean)
  extreme_grp  <- grp[grp$Mean == extreme_mean, ]
  
  if (nrow(extreme_grp) > 1) {
    best <- extreme_grp[extreme_grp$N == max(extreme_grp$N), ]
    if (nrow(best) > 1) best <- best[1, ]
  } else {
    best <- extreme_grp
  }
  
  o <- order(grp$Mean, decreasing = (direction == "higher"))
  if (nrow(grp) >= 2) {
    sec <- grp[o[2], ]
    sec_hap <- sec$Haplotype
    impr <- if (abs(sec$Mean) > .Machine$double.eps)
      ((best$Mean - sec$Mean) / abs(sec$Mean)) * 100 else NA
    dat <- test_res$data
    eff_col <- if (!is.null(test_res$test_col)) test_res$test_col else trait_col
    g1 <- dat[[eff_col]][dat$Haplotype == best$Haplotype]
    g2 <- dat[[eff_col]][dat$Haplotype == sec$Haplotype]
    psd <- sqrt(((length(g1) - 1) * stats::var(g1, na.rm = TRUE) +
                   (length(g2) - 1) * stats::var(g2, na.rm = TRUE)) /
                  (length(g1) + length(g2) - 2))
    cohens_d <- if (!is.na(psd) && psd > 0)
      abs(mean(g1, na.rm = TRUE) - mean(g2, na.rm = TRUE)) / psd else NA
    # Hedges' (1981) small-sample correction: Cohen's d is biased upward when
    # groups are small (e.g. 3 individuals per haplotype). J = 1 - 3 / (4 df - 1).
    if (isTRUE(small_sample_correction) && !is.na(cohens_d)) {
      n12 <- length(g1) + length(g2)
      if (n12 > 2L) cohens_d <- cohens_d * (1 - 3 / (4 * (n12 - 2) - 1))
    }
    dm <- best$Mean - sec$Mean
    dl <- dm - 1.96 * sqrt(best$SE^2 + sec$SE^2)
    du <- dm + 1.96 * sqrt(best$SE^2 + sec$SE^2)
  } else {
    sec_hap  <- NA_character_
    impr     <- NA_real_
    cohens_d <- NA_real_
    dm <- dl <- du <- NA_real_
  }
  
  gate_a <- isTRUE(test_res$p_value < alpha)
  gate_b <- isTRUE(test_res$permutation_p < alpha)
  gate_c <- !is.na(cohens_d) && cohens_d >= min_cohens_d
  
  gate_d <- TRUE
  if (!is.null(trait_threshold)) {
    gate_d <- if (direction == "higher")
      isTRUE(best$Mean >= trait_threshold)
    else isTRUE(best$Mean <= trait_threshold)
  }
  
  is_sup <- gate_a && gate_b && gate_c && gate_d
  
  fail <- c(
    if (!gate_a) sprintf("p>=%.2g", alpha),
    if (!gate_b) sprintf("perm-p>=%.2g", alpha),
    if (!gate_c) sprintf("d=%.2f < %.1f", cohens_d, min_cohens_d),
    if (!gate_d) sprintf("mean=%.3g does not cross trait_threshold=%.3g",
                         best$Mean, trait_threshold))
  
  freq_analysed_pct <- round(best$N / sum(grp$N) * 100, 2)
  freq_panel_pct <- if (!is.null(n_panel) && n_panel > 0)
    round(best$N / n_panel * 100, 2) else NA_real_
  
  list(
    haplotype = best$Haplotype,
    nt_sequence = best$NT_Sequence,
    aa_sequence = if ("AA_Sequence" %in% names(best)) best$AA_Sequence
    else NA,
    mean   = best$Mean,
    median = if ("Median" %in% names(best)) best$Median else NA,
    se     = if ("SE" %in% names(best)) best$SE else NA,
    ci_lower = if ("CI_Lower" %in% names(best)) best$CI_Lower else NA,
    ci_upper = if ("CI_Upper" %in% names(best)) best$CI_Upper else NA,
    n = best$N,
    freq_analysed_pct = freq_analysed_pct,
    freq_panel_pct    = freq_panel_pct,
    letter            = if ("Letter" %in% names(best))
      best[["Letter"]] else NA_character_,
    second_haplotype  = sec_hap,
    is_superior       = is_sup,
    gate_a_sig        = gate_a,
    gate_b_perm_sig   = gate_b,
    gate_c_effect     = gate_c,
    gate_d_threshold  = gate_d,
    superior_label    = if (is_sup) best$Haplotype else
      paste0("No superior [", paste(fail, collapse = "; "), "]"),
    improvement_pct   = impr,
    cohens_d          = cohens_d,
    diff_mean         = dm,
    diff_ci_lower     = dl,
    diff_ci_upper     = du)
}

#' Re-apply the superior-haplotype rule with different settings
#'
#' The pipeline calls a haplotype superior when four gates pass: a
#' nominal significance threshold on the trait association (A), a
#' permutation-based significance threshold (B), an effect-size
#' minimum on the best-versus-runner-up comparison (C), and — when
#' requested — an absolute trait value the best haplotype must reach
#' (D). `hbr_regate()` recomputes which genes pass these gates under
#' new thresholds, using the statistics already stored in the run, and
#' reports how the call changes.
#'
#' @details
#' Two calling conventions are supported.
#'
#' **Data.frame input.** When `x` is a data.frame carrying at least a
#' `P_value` column (and, for the effect-size gate, `Cohens_d` or
#' `Hedges_g`), the function returns a copy of the table with
#' `Adjusted_P`, `Significant_FDR`, `Significant_Nominal` and
#' `Is_Superior` recomputed. This is the legacy behaviour and is fully
#' backward compatible.
#'
#' **Run input.** When `x` is an `hbr_run` (from
#' [hbr_read_results()]) or anything that function accepts, all four
#' gates are recomputed — including the permutation gate and, when the
#' pipeline stored a trait threshold, the trait gate. Arguments left
#' at `NULL` inherit the run's stored values, so the default call
#' reproduces the pipeline exactly. The return is a classed object with
#' the per-gene calls, the effective settings, the genes gained and
#' lost versus the stored call, and a Jaccard overlap statistic.
#'
#' The PCs the pipeline used and the Hedges correction cannot be
#' changed here, because both alter the stored statistics; re-run the
#' pipeline for that.
#'
#' @param x A data.frame (legacy) or an `hbr_run` (new).
#' @param alpha Significance threshold for the nominal and permutation
#'   gates, and for Benjamini-Hochberg FDR when `require_fdr = TRUE`.
#'   `NULL` inherits the run's stored value.
#' @param min_cohens_d Minimum effect size for gate C. `NULL` inherits.
#' @param require_fdr If `TRUE`, a gene is superior only if it also
#'   passes BH-FDR across all tested genes. `NULL` inherits.
#' @param trait_threshold Absolute trait value gate D. `NA` removes the
#'   threshold. `NULL` inherits.
#' @param direction Kept for backward compatibility with the legacy
#'   signature; ignored. The run carries its own direction.
#' @param ... Ignored.
#'
#' @return When `x` is a data.frame: a `data.table` with the recomputed
#'   columns. When `x` is an `hbr_run`: an object of class
#'   `hbr_regate` with elements `genes`, `settings`, `overlap` and
#'   `check`.
#'
#' @seealso [hbr_read_results()] to build the run input,
#'   [haplobreed_analysis()] for the pipeline that produced it.
#'
#' @examples
#' # Legacy: data.frame input
#' summary_table <- data.frame(
#'   P_value  = c(0.001, 0.02, 0.5),
#'   Cohens_d = c(0.8, 0.4, 0.1))
#' hbr_regate(summary_table, alpha = 0.05, min_cohens_d = 0.3)
#'
#' # New: hbr_run input re-applies all four gates
#' toy <- system.file("extdata", "hgsi_toy", "trait_YLD",
#'                    package = "HaploBreedR")
#' if (nzchar(toy) && dir.exists(toy)) {
#'   run <- hbr_read_results(toy)
#'   rg  <- hbr_regate(run, alpha = 0.01, min_cohens_d = 0.5)
#'   rg$overlap$gained
#' }
#'
#' @export
hbr_regate <- function(x, alpha = NULL, min_cohens_d = NULL,
                       require_fdr = NULL, trait_threshold = NULL,
                       direction = NULL, ...) {
  
  # ---- Legacy path: data.frame input -------------------------------------
  if (is.data.frame(x)) {
    a  <- if (is.null(alpha)) 0.05 else alpha
    md <- if (is.null(min_cohens_d)) 0.5 else min_cohens_d
    fd <- if (is.null(require_fdr)) TRUE else require_fdr
    return(.hbr_regate_df(x, a, md, fd))
  }
  
  # ---- New path: hbr_run input -------------------------------------------
  run <- if (inherits(x, "hbr_run")) x else hbr_read_results(x)
  if (is.list(run) && !inherits(run, "hbr_run"))
    stop("supply one trait (a single hbr_run)", call. = FALSE)
  
  tb <- run$table
  g  <- run$gate
  need <- c("P_value", "Permutation_P", "Cohens_d", "Superior_Mean")
  miss <- setdiff(need, names(tb))
  if (length(miss))
    stop("the results table lacks: ", paste(miss, collapse = ", "),
         call. = FALSE)
  
  gate_calls <- function(a, md, thr, fdr) {
    p  <- as.numeric(tb$P_value)
    pp <- as.numeric(tb$Permutation_P)
    d  <- as.numeric(tb$Cohens_d)
    m  <- as.numeric(tb$Superior_Mean)
    ga <- !is.na(p) & p < a
    gb <- !is.na(pp) & pp < a
    gc <- !is.na(d) & d >= md
    gd <- if (is.finite(thr))
      !is.na(m) & (if (run$direction == "higher") m >= thr else m <= thr)
    else rep(TRUE, length(p))
    nominal <- ga & gb & gc & gd
    adj <- hbr_apply_fdr_correction(p, "BH")
    list(nominal = nominal,
         adj     = adj,
         superior = if (isTRUE(fdr))
           nominal & !is.na(adj) & adj < a else nominal)
  }
  
  stored <- list(alpha = g$alpha, min_d = g$min_d,
                 thr = g$trait_threshold, fdr = g$require_fdr)
  used <- list(
    alpha = if (is.null(alpha)) stored$alpha else alpha,
    min_d = if (is.null(min_cohens_d)) stored$min_d else min_cohens_d,
    thr   = if (is.null(trait_threshold)) stored$thr else trait_threshold,
    fdr   = if (is.null(require_fdr)) stored$fdr else require_fdr)
  
  new  <- gate_calls(used$alpha, used$min_d, used$thr, used$fdr)
  old  <- gate_calls(stored$alpha, stored$min_d, stored$thr, stored$fdr)
  orig <- if ("Is_Superior" %in% names(tb))
    as.logical(tb$Is_Superior) %in% TRUE else tb$is_superior
  
  genes <- data.frame(
    gene = tb$gene,
    Chr  = if ("Chr" %in% names(tb)) tb$Chr else NA_character_,
    P_value = as.numeric(tb$P_value),
    Adjusted_P = new$adj,
    Cohens_d = as.numeric(tb$Cohens_d),
    Superior_Haplotype = if ("Superior_Haplotype" %in% names(tb))
      tb$Superior_Haplotype else NA_character_,
    superior_stored = orig,
    superior_new    = new$superior,
    stringsAsFactors = FALSE)
  
  inter <- sum(orig & new$superior)
  uni   <- sum(orig | new$superior)
  
  structure(list(
    genes = genes,
    settings = list(stored = stored, used = used,
                    direction = run$direction),
    overlap = list(
      n_stored = sum(orig),
      n_new    = sum(new$superior),
      gained   = genes$gene[!orig & new$superior],
      lost     = genes$gene[orig & !new$superior],
      jaccard  = if (uni) inter / uni else NA_real_),
    check = list(
      n_disagree = sum(old$superior != orig),
      agrees     = all(old$superior == orig))),
    class = "hbr_regate")
}

#' Recompute the superior-haplotype call on a summary table
#'
#' Legacy worker: returns a copy of a summary data.frame with the
#' significance and superior flags recomputed. Used by the data.frame
#' branch of [hbr_regate()].
#'
#' @param dt A data.frame or data.table with a `P_value` column.
#' @param alpha Significance threshold.
#' @param min_cohens_d Minimum effect size.
#' @param require_fdr If `TRUE`, use BH-FDR-adjusted p-values.
#' @return A data.table with the recomputed columns.
#' @noRd
.hbr_regate_df <- function(dt, alpha, min_cohens_d, require_fdr) {
  dt <- data.table::as.data.table(data.table::copy(dt))
  if (!"P_value" %in% names(dt))
    stop("summary_table must contain a 'P_value' column", call. = FALSE)
  
  dt[, Adjusted_P := stats::p.adjust(P_value, method = "BH")]
  dt[, Significant_FDR := !is.na(Adjusted_P) & Adjusted_P < alpha]
  dt[, Significant_Nominal := !is.na(P_value) & P_value < alpha]
  
  d_col  <- if ("Cohens_d" %in% names(dt)) "Cohens_d" else "Hedges_g"
  eff_ok <- if (d_col %in% names(dt)) dt[[d_col]] >= min_cohens_d else TRUE
  
  sig_check <- if (require_fdr)
    dt$Significant_FDR else dt$Significant_Nominal
  dt[, Is_Superior := sig_check & eff_ok]
  dt
}

#' Print an `hbr_regate` result
#'
#' @param x An object of class `hbr_regate`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @examples
#' toy <- system.file("extdata", "hgsi_toy", "trait_YLD",
#'                    package = "HaploBreedR")
#' if (nzchar(toy) && dir.exists(toy)) {
#'   run <- hbr_read_results(toy)
#'   print(hbr_regate(run, alpha = 0.01))
#' }
#' @export
print.hbr_regate <- function(x, ...) {
  u <- x$settings$used
  o <- x$overlap
  cat(sprintf(
    "Re-gated (direction %s): alpha %s | min effect %s | ",
    x$settings$direction, .hbr_fmt(u$alpha), .hbr_fmt(u$min_d)))
  cat(sprintf("trait threshold %s | FDR %s\n",
              if (is.finite(u$thr)) u$thr else "none",
              .hbr_fmt(u$fdr)))
  cat(sprintf(
    "  stored call: %d superior | new call: %d superior | Jaccard %s\n",
    o$n_stored, o$n_new, .hbr_fmt(round(o$jaccard, 3))))
  if (length(o$gained))
    cat("  gained:", paste(utils::head(o$gained, 10), collapse = ", "),
        if (length(o$gained) > 10) "..." else "", "\n")
  if (length(o$lost))
    cat("  lost:  ", paste(utils::head(o$lost, 10), collapse = ", "),
        if (length(o$lost) > 10) "..." else "", "\n")
  if (!x$check$agrees)
    cat(sprintf(
      "  note: the stored settings do not reproduce the stored call ",
      "for %d gene(s); the table may have been edited.\n",
      x$check$n_disagree))
  invisible(x)
}

# R/hgsi-engine.R
# Two-kernel genomic prediction engine, discovery replay and cross-validation.
# Everything that uses phenotypes is estimated on training individuals only.

#' Open a BEDMatrix handle (cached, silent)
#'
#' Opens `plink_prefix` once per R session and returns the cached handle
#' on later calls. The two info lines that `BEDMatrix` prints while
#' reading `.fam` and `.bim` are captured internally so they never reach
#' the console.
#'
#' @param plink_prefix PLINK file prefix.
#' @return A `BEDMatrix` object.
#' @noRd
.hgsi_open_bed <- function(plink_prefix) {
  suppressMessages(
    BEDMatrix::BEDMatrix(plink_prefix, simple_names = TRUE))
}
#' Scale a kernel to mean diagonal 1
#'
#' @param K Square matrix.
#' @return `K / mean(diag(K))`.
#' @noRd
.scale_diag1 <- function(K) {
  d <- mean(diag(K))
  if (!is.finite(d) || d <= 0) stop("kernel has non-positive mean diagonal")
  K / d
}

# ---------- streamed GRM (BEDMatrix; VanRaden 2008 method 1) -----------------

#' Genomic relationship matrix, streamed
#'
#' VanRaden (2008) method 1, computed in SNP blocks from a PLINK file.
#' @param plink_prefix PLINK prefix.
#' @param maf_min Minimum minor allele frequency.
#' @param block_cols SNPs per block.
#' @param n_cores Workers for block computation.
#' @param verbose Print progress.
#' @param train_idx Rows used for allele frequencies (default all; training rows in strict CV).
#' @return List: `A` (unscaled cross-product), `d`, `p`, `keep`, `n`, `m`, `maf_min`, `snps`; `G = A / d`.
#' @noRd
.hgsi_build_grm <- function(plink_prefix, maf_min = 0.05, block_cols = 20000L,
                            n_cores = 1L, verbose = TRUE, train_idx = NULL) {
  if (!requireNamespace("BEDMatrix", quietly = TRUE))
    stop("BEDMatrix required. install.packages('BEDMatrix')", call. = FALSE)
  G <- suppressMessages(
    BEDMatrix::BEDMatrix(plink_prefix, simple_names = TRUE))
  n <- nrow(G); m <- ncol(G); snps <- colnames(G)
  if (is.null(snps)) stop("BEDMatrix returned no SNP names", call. = FALSE)
  if (is.null(train_idx)) train_idx <- seq_len(n)
  if (verbose) message(sprintf("  GRM: n=%d, m=%d (streaming)", n, m))

  blk <- function(Gb, cols) {
    x <- Gb[, cols, drop = FALSE]
    storage.mode(x) <- "double"
    xt <- x[train_idx, , drop = FALSE]
    nobs <- colSums(!is.na(xt))
    pb <- colSums(xt, na.rm = TRUE) / (2 * pmax(nobs, 1L))
    ok <- nobs > 0L & pmin(pb, 1 - pb) >= maf_min
    A_blk <- NULL
    if (any(ok)) {
      z <- sweep(x[, ok, drop = FALSE], 2, 2 * pb[ok])
      z[is.na(z)] <- 0
      A_blk <- tcrossprod(z)
    }
    list(A = A_blk, p = pb, keep = ok, cols = cols)
  }

  A <- matrix(0, n, n); p <- numeric(m); keep <- logical(m)
  blocks <- split(seq_len(m), ceiling(seq_len(m) / block_cols))
  if (n_cores > 1L && requireNamespace("future.apply", quietly = TRUE)) {
    future::plan(future::multisession, workers = n_cores)
    on.exit(future::plan(future::sequential), add = TRUE)
    res <- future.apply::future_lapply(blocks, function(cols) {
      blk(.hgsi_open_bed(plink_prefix), cols)
    }, future.seed = TRUE)
  } else {
    res <- vector("list", length(blocks))
    for (b in seq_along(blocks)) {
      res[[b]] <- blk(G, blocks[[b]])
      if (verbose && b %% 5L == 0L)
        message(sprintf("    block %d / %d", b, length(blocks)))
    }
  }
  for (r in res) {
    if (!is.null(r$A)) A <- A + r$A
    p[r$cols] <- r$p; keep[r$cols] <- r$keep
  }
  d <- 2 * sum(p[keep] * (1 - p[keep]))
  list(A = A, d = d, p = p, keep = keep, n = n, m = m, maf_min = maf_min,
       snps = snps)
}

#' Cross-product of a SNP subset
#'
#' @param plink_prefix PLINK prefix.
#' @param p Allele frequencies from [.hgsi_build_grm()].
#' @param keep SNPs retained by the MAF filter.
#' @param cols Integer column indices of the subset.
#' @param block_cols SNPs per block.
#' @return List: `A`, `d`, `cols` (retained columns).
#' @noRd
.hgsi_feature_A <- function(plink_prefix, p, keep, cols, block_cols = 20000L) {
  G <- suppressMessages(
    BEDMatrix::BEDMatrix(plink_prefix, simple_names = TRUE))
  n <- nrow(G)
  cols <- sort(unique(cols[!is.na(cols)]))
  cols <- cols[keep[cols]]
  A <- matrix(0, n, n)
  if (!length(cols)) return(list(A = A, d = 0, cols = integer(0)))
  for (s in seq(1L, length(cols), by = block_cols)) {
    cc <- cols[s:min(length(cols), s + block_cols - 1L)]
    x <- G[, cc, drop = FALSE]
    storage.mode(x) <- "double"
    z <- sweep(x, 2, 2 * p[cc]); z[is.na(z)] <- 0
    A <- A + tcrossprod(z)
  }
  list(A = A, d = 2 * sum(p[cols] * (1 - p[cols])), cols = cols)
}

#' Remainder kernel after removing feature SNPs
#'
#' Exact because the cross-product is a sum over SNPs.
#' @param grm Output of [.hgsi_build_grm()].
#' @param feat Output of [.hgsi_feature_A()] or `NULL`.
#' @return `(A - A_feat) / (d - d_feat)`.
#' @noRd
.hgsi_remainder <- function(grm, feat = NULL) {
  if (is.null(feat) || feat$d == 0) return(grm$A / grm$d)
  (grm$A - feat$A) / (grm$d - feat$d)
}

# ---------- spectral two-kernel REML (Kang et al. 2008) ----------------------

#' REML profile over the variance ratio
#'
#' Rotated data `eta` with eigenvalues `lam` (Kang et al. 2008); the genetic variance is profiled out.
#' @param eta Rotated residual vector.
#' @param lam Eigenvalues (non-negative).
#' @return List: `delta`, `s2g`, `s2e`, `negll`.
#' @noRd
.hgsi_reml_profile <- function(eta, lam) {
  nq <- length(eta)
  f <- function(logd) {
    dl <- exp(logd)
    s <- sum(eta^2 / (lam + dl)) / nq
    0.5 * (nq * log(2 * pi * s) + sum(log(lam + dl)) + nq)
  }
  o <- stats::optimize(f, c(-10, 10))
  dl <- exp(o$minimum)
  s2g <- sum(eta^2 / (lam + dl)) / nq
  list(delta = dl, s2g = s2g, s2e = dl * s2g, negll = o$objective)
}

#' Two-kernel REML
#'
#' Estimates the weight `w` of the feature kernel by grid search plus refinement; each `w` is fitted by spectral REML.
#' @param y Phenotypes (training).
#' @param X Fixed-effect design.
#' @param G_R Remainder kernel (training block).
#' @param K_H Feature kernel (training block) or `NULL` for GBLUP.
#' @param w_grid Grid for `w`.
#' @param refine Refine around the best grid value.
#' @return List: `w`, `s2g`, `s2e`, `s2R`, `s2H`, `negll`.
#' @noRd
.hgsi_reml2k <- function(y, X, G_R, K_H = NULL,
                         w_grid = seq(0, 1, by = 0.1), refine = TRUE) {
  n <- length(y)
  qrX <- qr(X); q <- qrX$rank
  if (q >= n) stop("no residual degrees of freedom")
  Qp <- qr.Q(qrX, complete = TRUE)[, seq.int(q + 1L, n), drop = FALSE]
  yq <- crossprod(Qp, y)
  one <- function(w) {
    K <- if (is.null(K_H)) G_R else (1 - w) * G_R + w * K_H
    e <- eigen(crossprod(Qp, K %*% Qp), symmetric = TRUE)
    eta <- as.numeric(crossprod(e$vectors, yq))
    r <- .hgsi_reml_profile(eta, pmax(e$values, 0))
    c(w = w, negll = r$negll, s2g = r$s2g, s2e = r$s2e)
  }
  grid <- if (is.null(K_H)) 0 else w_grid
  prof <- do.call(rbind, lapply(grid, one))
  best <- prof[which.min(prof[, "negll"]), ]
  if (refine && !is.null(K_H)) {
    lo <- max(0, best[["w"]] - 0.05); hi <- min(1, best[["w"]] + 0.05)
    if (hi > lo) {
      o <- stats::optimize(function(w) one(w)[["negll"]], c(lo, hi))
      cand <- one(o$minimum)
      if (cand[["negll"]] < best[["negll"]]) best <- cand
    }
  }
  w <- best[["w"]]
  list(w = w, s2g = best[["s2g"]], s2e = best[["s2e"]],
       s2R = (1 - w) * best[["s2g"]], s2H = w * best[["s2g"]],
       negll = best[["negll"]])
}

#' Fit on training individuals and predict validation individuals
#'
#' All phenotype-dependent quantities (REML, fixed effects, BLUP) use training individuals only; validation individuals are predicted through the cross-kernel. Reliability is 1 - PEV / (s2g * K_ii), with PEV including the fixed-effect correction.
#' @param y Phenotypes, `NA` allowed.
#' @param X Fixed-effect design (all individuals).
#' @param G_R Remainder kernel, mean diagonal 1 (all individuals).
#' @param X_H Centred haplotype-indicator design or `NULL`.
#' @param train,val Row indices.
#' @param w_grid Grid for the feature-kernel weight.
#' @return List: `pred`, `g`, `reliability`, `w`, `s2R`, `s2H`, `s2e`, `s2g`, `hap_effects`, `vir` (V^-1 r), `train`, `val`.
#' @noRd
.hgsi_fit_predict <- function(y, X, G_R, X_H = NULL, train, val,
                              w_grid = seq(0, 1, by = 0.1)) {
  train <- intersect(train, which(!is.na(y)))
  if (length(train) < 20L) stop("too few training observations")
  val <- intersect(val, seq_along(y))
  if (!length(val)) stop("no validation individuals")
  yt <- y[train]; Xt <- X[train, , drop = FALSE]; Xv <- X[val, , drop = FALSE]

  K_H_all <- NULL; cH <- NA_real_
  if (!is.null(X_H)) {
    cH <- mean(rowSums(X_H[train, , drop = FALSE]^2))      # training rows only
    if (!is.finite(cH) || cH <= 0) X_H <- NULL
    else K_H_all <- tcrossprod(X_H) / cH
  }
  gr_tr <- G_R[train, train, drop = FALSE]
  kh_tr <- if (is.null(K_H_all)) NULL else K_H_all[train, train, drop = FALSE]
  rm_ <- .hgsi_reml2k(yt, Xt, gr_tr, kh_tr, w_grid = w_grid)
  w <- rm_$w

  Kw_all <- if (is.null(K_H_all)) G_R else (1 - w) * G_R + w * K_H_all
  V <- rm_$s2g * Kw_all[train, train, drop = FALSE]
  diag(V) <- diag(V) + rm_$s2e
  Vi <- chol2inv(chol(V))
  XtVi <- crossprod(Xt, Vi)
  XtViXt <- XtVi %*% Xt
  beta <- solve(XtViXt, XtVi %*% yt)
  r <- yt - Xt %*% beta
  Vir <- Vi %*% r

  K_vt <- Kw_all[val, train, drop = FALSE]
  g_val <- as.numeric(rm_$s2g * (K_vt %*% Vir))
  pred_val <- as.numeric(Xv %*% beta + g_val)

  # PEV with the fixed-effect correction (P matrix); reliability = 1 - PEV/s2g*K_ii
  P <- Vi - Vi %*% Xt %*% solve(XtViXt) %*% XtVi
  kd_val <- diag(Kw_all)[val]
  pev_val <- rm_$s2g * kd_val - rm_$s2g^2 * rowSums((K_vt %*% P) * K_vt)
  rel_val <- pmin(pmax(1 - pev_val / (rm_$s2g * kd_val), 0), 1)

  hap_eff <- if (is.null(X_H)) NULL else
    as.numeric((rm_$s2H / cH) * crossprod(X_H[train, , drop = FALSE], Vir))

  list(pred = pred_val, g = g_val, reliability = as.numeric(rel_val),
       w = w, s2R = rm_$s2R, s2H = rm_$s2H, s2e = rm_$s2e, s2g = rm_$s2g,
       hap_effects = hap_eff, vir = as.numeric(Vir), train = train, val = val)
}

#' Fit on all phenotyped individuals and predict everyone
#'
#' @param y,X,G_R,X_H As in [.hgsi_fit_predict()].
#' @return As [.hgsi_fit_predict()].
#' @noRd
.hgsi_fit_all <- function(y, X, G_R, X_H = NULL) {
  .hgsi_fit_predict(y, X, G_R, X_H, train = which(!is.na(y)),
                    val = seq_along(y))
}

# ---------- strata (unsupervised, genotype-only) -----------------------------

#' Genotype clusters for stratified folds
#'
#' k-means on the leading kernel principal components (unsupervised, genotype only).
#' @param G Relationship matrix.
#' @param n_strata Number of clusters (1 = no strata).
#' @param seed Random seed.
#' @return Integer cluster per individual.
#' @noRd
.hgsi_auto_strata <- function(G, n_strata = 3L, seed = 1L) {
  n <- nrow(G)
  if (n_strata <= 1L) return(rep(1L, n))
  H <- diag(n) - matrix(1 / n, n, n)
  e <- eigen(H %*% G %*% H, symmetric = TRUE)
  k <- min(3L, sum(e$values > 1e-8))
  pcs <- e$vectors[, seq_len(k), drop = FALSE] %*%
    diag(sqrt(e$values[seq_len(k)]), k)
  km <- .with_seed(seed, stats::kmeans(pcs, centers = n_strata, nstart = 10L))
  as.integer(km$cluster)
}

# ---------- design matrix, discovery gate, folds, metrics --------------------

#' Centred haplotype-indicator design
#'
#' One column per haplotype with at least `min_count` carriers in `train`; centred on training frequencies. Codes containing `"N"` count as missing (mean).
#' @param codes Matrix from [.hgsi_build_codes()].
#' @param genes Genes to include.
#' @param train Training rows.
#' @param min_count Minimum carriers.
#' @return List: `X` (or `NULL`) and `info` (`gene`, `hap`, `freq`).
#' @noRd
.hgsi_hap_design <- function(codes, genes, train, min_count = 5L) {
  cols <- list(); info <- list()
  for (g in genes) {
    if (!g %in% colnames(codes)) next
    v <- codes[, g]
    v[!is.na(v) & grepl("N", v, fixed = TRUE)] <- NA
    tb <- table(v[train])
    keep <- names(tb)[tb >= min_count]
    if (length(keep) < 2L) next
    for (h in keep) {
      x <- as.numeric(v == h); mu <- mean(x[train], na.rm = TRUE)
      x <- x - mu; x[is.na(x)] <- 0
      cols[[length(cols) + 1L]] <- x
      info[[length(info) + 1L]] <- data.frame(gene = g, hap = h, freq = mu,
                                              stringsAsFactors = FALSE)
    }
  }
  if (!length(cols)) return(list(X = NULL, info = NULL))
  list(X = do.call(cbind, cols), info = do.call(rbind, info))
}

#' Replay of the pipeline's superior-haplotype gate
#'
#' Groups by haplotype (size and frequency filters), then per gene: one-way F-test (equals the pooled t-test for two groups) on the trait, residualised on the pipeline's PCs if it used any; best and runner-up haplotypes by raw group mean in the pipeline's direction; Hedges' g (or Cohen's d) of best vs runner-up on the residuals; optional trait threshold; BH across genes when FDR was required. The pipeline's permutation gate is not replayed. Uses `train` individuals only.
#' @param codes Matrix from [.hgsi_build_codes()].
#' @param y Phenotypes.
#' @param train Row indices.
#' @param gate Gate list (`direction`, `alpha`, `min_d`, `min_ind`, `min_hap_freq`, `require_fdr`, `trait_threshold`, `hedges`).
#' @param pcs Matrix of the pipeline's PCs (rows = individuals) or `NULL`.
#' @return Data.frame: `gene`, `n_groups`, `p_raw`, `p_adj`, `d`, `best_hap`, `is_superior`.
#' @noRd
.hgsi_gate_fold <- function(codes, y, train, gate, pcs = NULL) {
  train <- train[!is.na(y[train])]
  empty <- data.frame(gene = character(), n_groups = integer(),
                      p_raw = numeric(), p_adj = numeric(), d = numeric(),
                      best_hap = character(), is_superior = logical(),
                      stringsAsFactors = FALSE)
  if (length(train) < 20L) return(empty)
  sgn <- if (gate$direction == "higher") 1 else -1
  genes <- colnames(codes)
  out <- vector("list", length(genes))
  for (i in seq_along(genes)) {
    v <- codes[train, i]
    ok <- !is.na(v)
    tb <- table(v[ok])
    gk <- names(tb)[tb >= gate$min_ind & tb / sum(tb) >= gate$min_hap_freq]
    if (length(gk) < 2L) next
    sel <- which(ok & v %in% gk)
    tr_i <- train[sel]; vv <- v[sel]; yy <- y[tr_i]
    if (!is.null(pcs)) {
      cc <- stats::complete.cases(pcs[tr_i, , drop = FALSE])
      tr_i <- tr_i[cc]; vv <- vv[cc]; yy <- yy[cc]
    }
    if (length(yy) < 10L || length(unique(vv)) < 2L) next
    rr <- if (is.null(pcs)) yy else
      stats::lm.fit(cbind(1, pcs[tr_i, , drop = FALSE]), yy)$residuals

    gmean <- tapply(yy, vv, mean); gn <- tapply(yy, vv, length)
    ord <- order(-sgn * gmean, -gn)
    best <- names(gmean)[ord[1]]; second <- names(gmean)[ord[2]]

    k <- length(gmean); N <- length(rr)
    ss_w <- sum(tapply(rr, vv, function(z) sum((z - mean(z))^2)))
    ss_t <- sum((rr - mean(rr))^2)
    if (N - k < 1 || ss_w <= 0) next
    Fv <- ((ss_t - ss_w) / (k - 1)) / (ss_w / (N - k))
    pv <- stats::pf(Fv, k - 1, N - k, lower.tail = FALSE)

    a <- rr[vv == best]; b <- rr[vv == second]
    sp <- sqrt(((length(a) - 1) * stats::var(a) + (length(b) - 1) * stats::var(b)) /
                 (length(a) + length(b) - 2))
    d <- if (is.finite(sp) && sp > 0) abs(mean(a) - mean(b)) / sp else NA_real_
    if (isTRUE(gate$hedges) && is.finite(d))      # Hedges' g, as in the pipeline
      d <- d * (1 - 3 / (4 * (length(a) + length(b) - 2) - 1))
    thr_ok <- if (is.finite(gate$trait_threshold))
      (if (sgn > 0) gmean[[best]] >= gate$trait_threshold
       else gmean[[best]] <= gate$trait_threshold) else TRUE
    out[[i]] <- data.frame(gene = genes[i], n_groups = k, p_raw = pv, d = d,
                           best_hap = best, thr_ok = thr_ok,
                           stringsAsFactors = FALSE)
  }
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) return(empty)
  res <- do.call(rbind, out)
  res$p_adj <- if (isTRUE(gate$require_fdr))
    hbr_apply_fdr_correction(res$p_raw, "BH") else res$p_raw
  res$is_superior <- !is.na(res$p_adj) & res$p_adj < gate$alpha &
    !is.na(res$d) & res$d >= gate$min_d & res$thr_ok
  res[, c("gene", "n_groups", "p_raw", "p_adj", "d", "best_hap", "is_superior")]
}

#' Agreement of replayed and pipeline superior genes
#'
#' @param codes Matrix from [.hgsi_build_codes()].
#' @param y Phenotypes.
#' @param gate Gate list.
#' @param reference Named logical: gene -> superior in the pipeline.
#' @param pcs PCs matrix or `NULL`.
#' @return List: `agreement`, `kappa`, `jaccard`, `n_genes`, `table`, `best_hap`.
#' @noRd
.hgsi_gate_concordance <- function(codes, y, gate, reference, pcs = NULL) {
  res <- .hgsi_gate_fold(codes, y, which(!is.na(y)), gate, pcs)
  genes <- colnames(codes)
  rep_flag <- stats::setNames(genes %in% res$gene[res$is_superior], genes)
  ref <- as.logical(reference[genes])
  tab <- table(replay = factor(rep_flag, c(FALSE, TRUE)),
               reference = factor(ref, c(FALSE, TRUE)))
  po <- sum(diag(tab)) / sum(tab)
  pe <- sum(rowSums(tab) * colSums(tab)) / sum(tab)^2
  either <- rep_flag | ref
  jacc <- if (any(either)) sum(rep_flag & ref) / sum(either) else NA_real_
  list(agreement = po,
       kappa = if (isTRUE(all.equal(pe, 1))) NA_real_ else (po - pe) / (1 - pe),
       jaccard = jacc, n_genes = length(genes), table = tab,
       best_hap = stats::setNames(res$best_hap, res$gene))
}

#' Stratified random fold assignment
#'
#' @param strata Stratum per individual.
#' @param k Number of folds.
#' @param seed Random seed.
#' @return Integer fold per individual.
#' @noRd
.hgsi_make_folds <- function(strata, k, seed = 1L) {
  .with_seed(seed, {
    fold <- integer(length(strata))
    for (s in unique(strata)) {
      idx <- which(strata == s)
      idx <- idx[sample.int(length(idx))]
      fold[idx] <- rep_len(sample.int(k), length(idx))
    }
    fold
  })
}

#' Prediction metrics for one fold
#'
#' `pred` must include the intercept so the RMSE is valid. `gain` is the selection differential of the top fraction in the pipeline's direction.
#' @param pred,obs Predictions and observations.
#' @param top_frac Fraction selected.
#' @param sgn `+1` for higher, `-1` for lower.
#' @return Named numeric: `r`, `slope` (obs on pred), `rmse`, `gain`.
#' @noRd
.hgsi_metrics <- function(pred, obs, top_frac, sgn) {
  r <- if (stats::sd(pred) > 0 && stats::sd(obs) > 0) stats::cor(pred, obs) else NA_real_
  slope <- if (stats::var(pred) > 0) stats::cov(pred, obs) / stats::var(pred) else NA_real_
  rmse <- sqrt(mean((obs - pred)^2, na.rm = TRUE))
  nsel <- max(1L, ceiling(top_frac * length(obs)))
  o <- order(sgn * pred, decreasing = TRUE)[seq_len(nsel)]
  c(r = r, slope = slope, rmse = rmse, gain = sgn * (mean(obs[o]) - mean(obs)))
}

# ---------- cross-validation with discovery replay ---------------------------

#' Cross-validation with discovery replay
#'
#' Fits M0 (GBLUP), M2 (features discovered on all data; optimistic) and M3 (features re-discovered in each training fold; honest). A fold with no discovered gene falls back to M0; folds are never dropped. With `strict = TRUE` the GRM allele frequencies are recomputed per training set.
#' @param y Phenotypes.
#' @param ids Individual IDs.
#' @param grm Output of [.hgsi_build_grm()].
#' @param plink_prefix PLINK prefix.
#' @param codes Matrix from [.hgsi_build_codes()].
#' @param gene_snps Named list gene -> SNP IDs.
#' @param gate Gate list.
#' @param pcs PCs matrix or `NULL`.
#' @param strata Fold strata.
#' @param repeats,k Repeats and folds.
#' @param top_frac Fraction selected for `gain`.
#' @param min_count Minimum carriers per haplotype indicator.
#' @param seed Random seed.
#' @param strict Recompute the GRM within each training set.
#' @param verbose Print progress.
#' @return List: `metrics` (per fold and method), `preds`, `naive_discovery`.
#' @noRd
.hgsi_cv <- function(y, ids, grm, plink_prefix, codes, gene_snps, gate,
                     pcs = NULL, strata = NULL, repeats = 10L, k = 5L,
                     top_frac = 0.1, min_count = 5L, seed = 1L,
                     strict = FALSE, verbose = TRUE) {
  n <- length(y)
  if (is.null(strata)) strata <- rep(1L, n)
  sgn <- if (gate$direction == "higher") 1 else -1
  obs <- which(!is.na(y))
  G_all <- .scale_diag1(grm$A / grm$d)
  X <- matrix(1, n, 1L)

  snp_cols <- function(genes) {
    cc <- match(unlist(gene_snps[genes], use.names = FALSE), grm$snps)
    sort(unique(cc[!is.na(cc)]))
  }
  cache <- new.env(parent = emptyenv())
  remainder_for <- function(genes, grm_use, use_cache) {
    key <- paste(sort(genes), collapse = "|")
    if (use_cache && !is.null(cache[[key]])) return(cache[[key]])
    feat <- .hgsi_feature_A(plink_prefix, grm_use$p, grm_use$keep, snp_cols(genes))
    Gr <- .scale_diag1(.hgsi_remainder(grm_use, feat))
    if (use_cache) {
      if (length(ls(cache)) > 100L) rm(list = ls(cache), envir = cache)
      assign(key, Gr, envir = cache)
    }
    Gr
  }
  # two-kernel fit for a gene set; falls back to M0 when no usable features
  fit_set <- function(genes, G_f, grm_f, tr, va, pred0, use_cache) {
    if (!length(genes)) return(list(pred = pred0, w = NA_real_, ng = 0L, fb = TRUE))
    des <- .hgsi_hap_design(codes, genes, tr, min_count)
    if (is.null(des$X)) return(list(pred = pred0, w = NA_real_, ng = 0L, fb = TRUE))
    f <- .hgsi_fit_predict(y, X, remainder_for(genes, grm_f, use_cache),
                           des$X, tr, va)
    list(pred = f$pred, w = f$w, ng = length(unique(des$info$gene)), fb = FALSE)
  }

  # naive discovery on ALL phenotyped individuals (the optimistic comparator)
  naive_tab <- .hgsi_gate_fold(codes, y, obs, gate, pcs)
  naive_genes <- naive_tab$gene[naive_tab$is_superior]

  metrics <- list(); preds <- list()
  for (rep_i in seq_len(repeats)) {
    fold <- .hgsi_make_folds(strata, k, seed = seed + rep_i)
    for (f in seq_len(k)) {
      va <- which(fold == f & !is.na(y)); tr <- which(fold != f & !is.na(y))
      if (length(va) < 5L || length(tr) < 20L) next
      if (isTRUE(strict)) {
        grm_f <- .hgsi_build_grm(plink_prefix, maf_min = grm$maf_min,
                                 train_idx = tr, verbose = FALSE)
        G_f <- .scale_diag1(grm_f$A / grm_f$d); use_cache <- FALSE
      } else {
        grm_f <- grm; G_f <- G_all; use_cache <- TRUE
      }
      pred0 <- .hgsi_fit_predict(y, X, G_f, NULL, tr, va)$pred
      m2 <- fit_set(naive_genes, G_f, grm_f, tr, va, pred0, use_cache)
      fold_tab <- .hgsi_gate_fold(codes, y, tr, gate, pcs)
      m3 <- fit_set(fold_tab$gene[fold_tab$is_superior], G_f, grm_f, tr, va,
                    pred0, use_cache)

      for (mm in c("M0", "M2", "M3")) {
        o <- switch(mm, M0 = list(pred = pred0, w = NA_real_, ng = 0L, fb = FALSE),
                    M2 = m2, M3 = m3)
        mt <- .hgsi_metrics(o$pred, y[va], top_frac, sgn)
        metrics[[length(metrics) + 1L]] <- data.frame(
          rep = rep_i, fold = f, method = mm, ntrain = length(tr), ntest = length(va),
          r = mt[["r"]], slope = mt[["slope"]], rmse = mt[["rmse"]], gain = mt[["gain"]],
          w_hat = o$w, n_genes = o$ng, fallback = o$fb, stringsAsFactors = FALSE)
        preds[[length(preds) + 1L]] <- data.frame(
          rep = rep_i, fold = f, method = mm, id = ids[va], y = y[va],
          pred = o$pred, stringsAsFactors = FALSE)
      }
    }
    if (verbose) message("    CV repeat ", rep_i, " / ", repeats)
  }
  if (!length(metrics)) stop("no CV fold could be evaluated (too few individuals?)")
  list(metrics = do.call(rbind, metrics), preds = do.call(rbind, preds),
       naive_discovery = naive_tab)
}

#' Summarise CV metrics by trait and method
#'
#' @param cv_l Named list of [.hgsi_cv()] results.
#' @return Data.frame of mean `r`, `slope`, `rmse`, `gain`, genes used and fallback fraction.
#' @noRd
.hgsi_cv_summary <- function(cv_l) {
  out <- do.call(rbind, lapply(names(cv_l), function(tr) {
    m <- cv_l[[tr]]$metrics
    do.call(rbind, lapply(split(m, m$method), function(d)
      data.frame(trait = tr, method = d$method[1], n_folds = nrow(d),
                 r = mean(d$r, na.rm = TRUE), r_sd = stats::sd(d$r, na.rm = TRUE),
                 slope = mean(d$slope, na.rm = TRUE), rmse = mean(d$rmse, na.rm = TRUE),
                 gain = mean(d$gain, na.rm = TRUE), mean_genes = mean(d$n_genes),
                 fallback_frac = mean(d$fallback), stringsAsFactors = FALSE)))
  }))
  rownames(out) <- NULL
  out
}

# ---------- comparison statistics --------------------------------------------

#' Paired comparison of two methods in CV
#'
#' Corrected resampled t-test (Nadeau and Bengio 2003; variance factor 1/J + n_test/n_train, J-1 df) plus a paired bootstrap over individuals on repeat-averaged predictions. Neither is conclusive alone.
#' @param cv Output of [.hgsi_cv()].
#' @param baseline,method Methods compared (`method - baseline`).
#' @param metric Column of `cv$metrics`.
#' @param boot Bootstrap resamples.
#' @param seed Random seed.
#' @return One-row data.frame, or `NULL` if fewer than 3 folds.
#' @noRd
.hgsi_compare <- function(cv, baseline = "M0", method = "M3",
                          metric = "r", boot = 2000L, seed = 1L,
                          top_frac = 0.1) {
  m <- cv$metrics; pr <- cv$preds
  a <- m[m$method == method, c("rep", "fold", metric, "ntrain", "ntest")]
  b <- m[m$method == baseline, c("rep", "fold", metric)]
  j <- merge(a, b, by = c("rep", "fold"), suffixes = c(".a", ".b"))
  dv <- j[[paste0(metric, ".a")]] - j[[paste0(metric, ".b")]]
  ok <- is.finite(dv); j <- j[ok, , drop = FALSE]; dv <- dv[ok]; J <- length(dv)
  if (J < 3L) return(NULL)
  ratio <- mean(j$ntest) / mean(j$ntrain)
  se <- sqrt((1 / J + ratio) * stats::var(dv))
  tval <- if (is.finite(se) && se > 0) mean(dv) / se else NA_real_
  pval <- if (is.finite(tval)) 2 * stats::pt(-abs(tval), df = J - 1L) else NA_real_
  pa <- pr[pr$method == method, ]; pb <- pr[pr$method == baseline, ]
  ida <- tapply(pa$pred, pa$id, mean); idb <- tapply(pb$pred, pb$id, mean)
  yy <- tapply(pa$y, pa$id, mean)
  common <- intersect(names(ida), names(idb))
  stat_fn <- switch(metric,
                    r     = function(pr, y) stats::cor(pr, y),
                    rmse  = function(pr, y) sqrt(mean((y - pr)^2)),
                    gain  = function(pr, y) {
                      sel <- order(-pr)[seq_len(max(1L, round(top_frac * length(pr))))]
                      mean(y[sel]) - mean(y)
                    },
                    slope = function(pr, y) stats::coef(stats::lm(y ~ pr))[2],
                    function(pr, y) stats::cor(pr, y)
  )
  bd <- .with_seed(seed, replicate(boot, {
    s <- sample(common, replace = TRUE)
    stat_fn(ida[s], yy[s]) - stat_fn(idb[s], yy[s])
  }))
  data.frame(method = method, baseline = baseline, metric = metric,
             mean_diff = mean(dv), se = se, t = tval, df = J - 1L, p = pval,
             boot_lo = unname(stats::quantile(bd, 0.025, na.rm = TRUE)),
             boot_hi = unname(stats::quantile(bd, 0.975, na.rm = TRUE)),
             stringsAsFactors = FALSE)
}

#' M2 vs M0, M3 vs M0 and M2 vs M3
#'
#' A positive M2 vs M3 difference is the optimism of discovery on all data.
#' @param cv Output of [.hgsi_cv()].
#' @param boot Bootstrap resamples.
#' @param seed Random seed.
#' @return Data.frame with one row per comparison.
#' @noRd
.hgsi_compare_all <- function(cv, boot = 2000L, seed = 1L, top_frac = 0.1) {
  rbind(.hgsi_compare(cv, "M0", "M2", boot = boot, seed = seed, top_frac = top_frac),
        .hgsi_compare(cv, "M0", "M3", boot = boot, seed = seed, top_frac = top_frac),
        .hgsi_compare(cv, "M3", "M2", boot = boot, seed = seed, top_frac = top_frac))
}
# ---------- index ---------------------------------------------------------------

#' Selection-index weights
#'
#' `sign * weight / genetic SD` per trait; the sign is `-1` for traits the pipeline ran with direction `"lower"`, so a larger index is always better.
#' @param weights Non-negative emphasis per genetic SD.
#' @param directions `"higher"` or `"lower"` per trait.
#' @param sd_g Genetic SD per trait (square root of the REML genetic variance).
#' @return Numeric vector of weights.
#' @noRd
.hgsi_index_weights <- function(weights, directions, sd_g) {
  as.numeric(weights) * ifelse(as.character(directions) == "lower", -1, 1) /
    as.numeric(sd_g)
}


#' Streamed Genomic Relationship Matrix
#'
#' Computes the VanRaden (2008) Method 1 genomic relationship matrix in 
#' streaming blocks directly from binary PLINK files.
#'
#' @param plink_prefix Character. Path prefix to PLINK binary files.
#' @param maf_min Numeric. Minimum minor allele frequency (default 0.05).
#' @param block_cols Integer. SNPs evaluated per memory chunk (default 20,000).
#' @param n_cores Integer. Worker threads for parallel block evaluation.
#' @param verbose Logical. Print progress.
#' @param train_idx Optional integer vector of individuals for allele frequency estimation.
#'
#' @return List containing `A`, `d`, `p`, `keep`, `n`, `m`, and `snps`.
#' @examples
#' toy <- system.file("extdata", "toy", package = "HaploBreedR")
#' if (requireNamespace("BEDMatrix", quietly = TRUE)) {
#'   grm <- hbr_build_grm(file.path(toy, "toy"), verbose = FALSE)
#'   length(grm$snps)
#' }
#' @export
hbr_build_grm <- function(plink_prefix, maf_min = 0.05, block_cols = 20000L,
                          n_cores = 1L, verbose = TRUE, train_idx = NULL) {
  .hgsi_build_grm(plink_prefix, maf_min = maf_min, block_cols = block_cols,
                  n_cores = n_cores, verbose = verbose, train_idx = train_idx)
}

#' Spectral Two-Kernel REML
#'
#' Estimates variance components for a two-kernel mixed model using Kang et al. (2008)
#' spectral profiling over a convex kernel mixture weight.
#'
#' @param y Numeric vector of phenotypes.
#' @param X Matrix of fixed-effect covariates.
#' @param G_R Square relationship matrix (background polygenic kernel).
#' @param K_H Square relationship matrix (feature haplotype kernel) or NULL.
#' @param w_grid Numeric grid for kernel mixing weight w in [0, 1].
#' @param refine Logical. Perform golden section search refinement around best grid point.
#'
#' @return List containing `w`, `s2g`, `s2e`, `s2R`, `s2H`, and `negll`.
#' @examples
#' set.seed(1)
#' n <- 30
#' X <- cbind(1, rnorm(n))
#' Z <- matrix(rnorm(n * 20), n, 20)
#' G <- tcrossprod(Z) / 20
#' y <- rnorm(n)
#' hbr_reml_two_kernel(y, X, G)
#' @export
hbr_reml_two_kernel <- function(y, X, G_R, K_H = NULL,
                                w_grid = seq(0, 1, by = 0.1), refine = TRUE) {
  .hgsi_reml2k(y, X, G_R, K_H = K_H, w_grid = w_grid, refine = refine)
}

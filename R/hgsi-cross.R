# R/hgsi-cross.R
# Cross usefulness, haplotype pyramiding and set selection for inbred lines.
# All values are on the index scale, where larger is better; the index weights
# carry the pipeline direction.

#' Bundled rice genetic map anchors
#'
#' Returns the Nagamura et al. (1997) rice linkage-map anchors (n ~ 1,450
#' RFLP markers with both physical and genetic positions; total map length
#' 1,510 cM) that ship with the package. The engine interpolates these
#' anchors onto the physical positions of any PLINK `.bim` at run time, so
#' the same table is valid for any rice panel. Use the result as the
#' `genetic_map` argument of [hbr_advanced_decision()] when you want to
#' inspect it, subset it, or combine it with additional markers.
#'
#' @return A data.frame with columns `chr` (character), `bp` (numeric,
#'   IRGSP-1.0 coordinates) and `cm` (numeric, Haldane/Kosambi-agnostic;
#'   the engine uses Haldane's map function). The `source` attribute
#'   records the reference and marker count.
#'
#' @examples
#' m <- hbr_rice_genetic_map()
#' head(m)
#' attr(m, "source")
#'
#' @export
hbr_rice_genetic_map <- function() {
  path <- system.file("extdata", "rice_genetic_map_anchors.rds",
                      package = "HaploBreedR")
  if (!nzchar(path) || !file.exists(path))
    stop("bundled rice genetic map anchors not found; ",
         "reinstall HaploBreedR or supply `genetic_map` explicitly.",
         call. = FALSE)
  m <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(m))
    stop("bundled rice genetic map anchors could not be read", call. = FALSE)
  m <- as.data.frame(m)
  if (!all(c("chr", "bp", "cm") %in% names(m)))
    stop("bundled rice genetic map is malformed (needs chr, bp, cm)",
         call. = FALSE)
  m$chr <- as.character(m$chr)
  m$bp  <- as.numeric(m$bp)
  m$cm  <- as.numeric(m$cm)
  m
}
#' Haldane map function
#'
#' @param cM Distance in centimorgans.
#' @return Recombination fraction.
#' @noRd
.hgsi_haldane <- function(cM) 0.5 * (1 - exp(-2 * cM / 100))

#' Read a PLINK `.bim` file
#'
#' @param plink_prefix PLINK prefix.
#' @return Data.frame: `chr`, `snp`, `cm`, `bp`, `a1`, `a2`.
#' @noRd
.hgsi_read_bim <- function(plink_prefix) {
  b <- data.table::fread(paste0(plink_prefix, ".bim"), header = FALSE,
                         data.table = FALSE)
  names(b)[seq_len(6)] <- c("chr", "snp", "cm", "bp", "a1", "a2")
  b$chr <- as.character(b$chr); b$snp <- as.character(b$snp)
  b
}

#' Genetic position (cM) of each SNP
#'
#' Priority:
#' \enumerate{
#'   \item A data.frame in `map` (interpolated per chromosome).
#'   \item The literal strings `"rice"`, `"harushima"`, `"rgp"` or
#'     `"default"` — resolves to the bundled anchors
#'     ([hbr_rice_genetic_map()]).
#'   \item `NULL` and no `cM_per_Mb` — also falls back to the bundled
#'     anchors. Pass `map = "none"` to opt out.
#'   \item The `"none"` string — skip the map; try `cM_per_Mb`.
#'   \item A non-zero, non-constant cM column in the `.bim`.
#'   \item The constant `cM_per_Mb`.
#' }
#' Returns `NULL` if none of the sources applies.
#'
#' @param bim Output of [.hgsi_read_bim()].
#' @param cM_per_Mb Constant recombination rate or `NULL`.
#' @param map Data.frame with `chr`, `bp`, `cm`; the strings `"rice"`,
#'   `"harushima"`, `"rgp"`, `"default"` (use bundled); `"none"` (skip);
#'   or `NULL` (use bundled when nothing else is available).
#' @return List: `cm`, `source`; `NULL` if no position source exists.
#' @noRd
.hgsi_genetic_pos <- function(bim, cM_per_Mb = NULL, map = NULL) {
  
  # --- Resolve the map argument -------------------------------------------
  explicit_none <- is.character(map) && length(map) == 1L &&
    tolower(map) == "none"
  
  if (explicit_none) {
    map <- NULL
  } else if (is.character(map) && length(map) == 1L &&
             tolower(map) %in% c("default", "rice", "harushima", "rgp")) {
    map <- tryCatch(hbr_rice_genetic_map(), error = function(e) NULL)
    if (is.null(map))
      warning("bundled rice genetic map not available; ",
              "supply `genetic_map` or `cM_per_Mb` explicitly.",
              call. = FALSE)
  } else if (is.null(map) && is.null(cM_per_Mb)) {
    # No source requested at all: use the bundled map if it is present
    map <- tryCatch(hbr_rice_genetic_map(), error = function(e) NULL)
    if (!is.null(map))
      message("  using bundled rice genetic map (Harushima 1998 / RGP); ",
              "pass `genetic_map = \"none\"` to disable")
  }
  
  # --- Interpolate when a map is available --------------------------------
  if (!is.null(map) && is.data.frame(map) &&
      all(c("chr", "bp", "cm") %in% names(map))) {
    cm <- rep(NA_real_, nrow(bim))
    map$chr <- as.character(map$chr)
    for (cc in unique(bim$chr)) {
      mm <- map[map$chr == cc, , drop = FALSE]
      ii <- which(bim$chr == cc)
            mm <- mm[order(mm$bp), , drop = FALSE]
      mm <- mm[!duplicated(mm$bp), , drop = FALSE]
      if (nrow(mm) >= 2L)
        cm[ii] <- stats::approx(mm$bp, mm$cm,
                                xout = bim$bp[ii], rule = 2)$y
    }
    if (mean(is.finite(cm)) > 0.9) {
      src <- attr(map, "source")
      if (is.null(src)) src <- "user genetic map"
      return(list(cm = cm, source = src))
    }
  }
  
  # --- Fallbacks (unchanged) ----------------------------------------------
  if (any(bim$cm != 0) && length(unique(bim$cm)) > 1L)
    return(list(cm = bim$cm, source = "cM column of .bim"))
  
  if (!is.null(cM_per_Mb))
    return(list(cm = bim$bp / 1e6 * cM_per_Mb,
                source = sprintf(
                  "constant %.2f cM/Mb (sensitivity analysis advised)",
                  cM_per_Mb)))
  
  NULL
}

#' Group SNPs into genetic bins
#'
#' @param chr Chromosome per SNP.
#' @param cm Genetic position per SNP.
#' @param bin_cM Bin width in cM.
#' @return List: `bin_of`, `chr`, `center`, `n`.
#' @noRd
.hgsi_make_bins <- function(chr, cm, bin_cM = 1) {
  key <- paste(chr, floor(cm / bin_cM))
  ukey <- unique(key[order(chr, cm)])
  bin_of <- match(key, ukey)
  first <- match(ukey, key)
  center <- as.numeric(tapply(cm, factor(bin_of, levels = seq_along(ukey)), mean))
  list(bin_of = bin_of, chr = chr[first], center = center, n = length(ukey))
}

#' SNP effects of the remainder kernel
#'
#' For kernel `A_R / mean(diag(A_R))`, `u_k = s2R / mean(diag(A_R)) * z_k' V^-1 r`. Feature SNPs are set to zero because the feature kernel carries them.
#' @param plink_prefix PLINK prefix.
#' @param grm Output of [.hgsi_build_grm()].
#' @param V Matrix (individuals x traits) of V^-1 r, zero for non-training rows.
#' @param scale Per-trait `s2R / mean(diag(A_R))`.
#' @param exclude List (per trait) of feature SNP columns.
#' @param block_cols SNPs per block.
#' @return Matrix SNPs x traits.
#' @noRd
.hgsi_snp_effects <- function(plink_prefix, grm, V, scale, exclude = NULL,
                              block_cols = 20000L) {
  G <- suppressMessages(BEDMatrix::BEDMatrix(plink_prefix, simple_names = TRUE))
  Tn <- ncol(V); U <- matrix(0, grm$m, Tn)
  for (s in seq(1L, grm$m, by = block_cols)) {
    cols <- s:min(grm$m, s + block_cols - 1L)
    ok <- grm$keep[cols]
    if (!any(ok)) next
    cc <- cols[ok]
    x <- G[, cc, drop = FALSE]; storage.mode(x) <- "double"
    z <- sweep(x, 2, 2 * grm$p[cc]); z[is.na(z)] <- 0
    cp <- sweep(crossprod(z, V), 2, scale, "*")
    if (!is.null(exclude))
      for (t in seq_len(Tn)) cp[cc %in% exclude[[t]], t] <- 0
    U[cc, ] <- cp
  }
  U
}

#' Per-parent, per-bin sums of effect x dosage
#'
#' Missing genotypes are replaced by `2p`. The difference of two parents in a bin is the segregating effect of that bin.
#' @param plink_prefix PLINK prefix.
#' @param grm Output of [.hgsi_build_grm()].
#' @param U SNP effects (SNPs x traits).
#' @param parent_idx Parent row indices.
#' @param bin_of Bin per SNP.
#' @param n_bins Number of bins.
#' @param block_cols SNPs per block.
#' @return Array parents x bins x traits.
#' @noRd
.hgsi_parent_bin_sums <- function(plink_prefix, grm, U, parent_idx, bin_of, n_bins,
                                  block_cols = 20000L) {
  G <- suppressMessages(BEDMatrix::BEDMatrix(plink_prefix, simple_names = TRUE))
  Tn <- ncol(U); P <- length(parent_idx)
  S <- array(0, dim = c(P, n_bins, Tn))
  for (s in seq(1L, grm$m, by = block_cols)) {
    cols <- s:min(grm$m, s + block_cols - 1L)
    ok <- grm$keep[cols] & rowSums(abs(U[cols, , drop = FALSE])) > 0
    if (!any(ok)) next
    cc <- cols[ok]
    x <- G[parent_idx, cc, drop = FALSE]; storage.mode(x) <- "double"
    miss <- is.na(x)
    if (any(miss)) x[miss] <- (2 * matrix(grm$p[cc], nrow(x), ncol(x), byrow = TRUE))[miss]
    b <- bin_of[cc]
    for (t in seq_len(Tn)) {
      agg <- rowsum(t(sweep(x, 2, U[cc, t], "*")), b, reorder = FALSE)
      idx <- as.integer(rownames(agg))
      S[, idx, t] <- S[, idx, t] + t(agg)
    }
  }
  S
}

#' Add haplotype effects at gene positions
#'
#' Haplotype-indicator effects of the feature kernel enter as extra locus effects in the gene's bin. A parent with a missing or heterozygous code gets the expected effect (effect x frequency), matching the centred design.
#' @param S Array from [.hgsi_parent_bin_sums()].
#' @param t_idx Trait index.
#' @param codes Haplotype-code matrix.
#' @param parent_ids Parent IDs (rows of `S`).
#' @param hap_info Data.frame `gene`, `hap`, `freq`.
#' @param hap_eff Haplotype effects.
#' @param gene_bin Named integer: gene -> bin.
#' @return Updated `S`.
#' @noRd
.hgsi_add_gene_effects <- function(S, t_idx, codes, parent_ids, hap_info, hap_eff,
                                   gene_bin) {
  if (is.null(hap_info) || is.null(hap_eff)) return(S)
  for (g in unique(hap_info$gene)) {
    if (!g %in% colnames(codes) || is.na(gene_bin[g])) next
    ii <- which(hap_info$gene == g)
    cd <- codes[parent_ids, g]
    miss <- is.na(cd) | grepl("N", cd, fixed = TRUE)
    add <- numeric(length(parent_ids))
    for (i in ii) add[!miss & cd == hap_info$hap[i]] <- hap_eff[i]
    add[miss] <- sum(hap_eff[ii] * hap_info$freq[ii])
    S[, gene_bin[g], t_idx] <- S[, gene_bin[g], t_idx] + add
  }
  S
}

#' Usefulness of crosses between inbred lines
#'
#' Progeny index variance is the sum over chromosomes of `D' C D`, with `D` the parents' bin-wise difference and `Cov(x_i, x_j) = (1 - 2r*)/4`, `r* = r` (DH) or `2r/(1 + 2r)` (RIL by selfing; Haldane and Waddington 1931), `r` from Haldane's map function. Usefulness is `mean + i * h * sd` (Schnell and Utz 1975). The index weights carry the pipeline direction, so larger is always better. Ridge/GBLUP effects are shrunk, so the SD is biased low: use usefulness to rank crosses, not as a calibrated variance.
#' @param S Array from [.hgsi_parent_bin_sums()].
#' @param gebv Parental GEBVs (parents x traits).
#' @param bins Output of [.hgsi_make_bins()].
#' @param w_eff Effective index weights.
#' @param pairs Two-column matrix of parent positions in `S`.
#' @param sel_frac Selected fraction (gives `i`).
#' @param h_sel Square root of heritability of selection in the progeny (1 = true genetic values).
#' @param mating `"DH"` or `"RIL"`.
#' @return Data.frame: `a`, `b`, `mid_index`, `sd_index`, `uc`.
#' @noRd
.hgsi_cross_uc <- function(S, gebv, bins, w_eff, pairs, sel_frac = 0.1, h_sel = 1,
                           mating = c("DH", "RIL")) {
  mating <- match.arg(mating)
  blocks <- split(seq_len(bins$n), bins$chr)
  Cb <- lapply(blocks, function(ix) {
    r <- .hgsi_haldane(abs(outer(bins$center[ix], bins$center[ix], "-")))
    if (mating == "RIL") r <- 2 * r / (1 + 2 * r)
    0.25 * (1 - 2 * r)
  })
  i_sel <- stats::dnorm(stats::qnorm(1 - sel_frac)) / sel_frac
  Tn <- dim(S)[3]
  out <- data.frame(a = pairs[, 1], b = pairs[, 2], mid_index = NA_real_,
                    sd_index = NA_real_, uc = NA_real_)
  for (k in seq_len(nrow(pairs))) {
    a <- pairs[k, 1]; b <- pairs[k, 2]
    D <- matrix(S[a, , ] - S[b, , ], nrow = bins$n, ncol = Tn)
    Dw <- as.numeric(D %*% w_eff)
    sig2 <- 0
    for (nm in names(blocks)) {
      ix <- blocks[[nm]]
      sig2 <- sig2 + drop(crossprod(Dw[ix], Cb[[nm]] %*% Dw[ix]))
    }
    mu <- sum(w_eff * (gebv[a, ] + gebv[b, ]) / 2)
    out$mid_index[k] <- mu
    out$sd_index[k] <- sqrt(max(sig2, 0))
    out$uc[k] <- mu + i_sel * h_sel * out$sd_index[k]
  }
  out
}

#' Probability of pyramiding target haplotypes in a DH progeny
#'
#' `1/2 * prod(1 - r)` or `r` over adjacent target loci per chromosome, depending on whether the same parent donates both; for RIL, `r` becomes `2r/(1 + 2r)` (Markov approximation).
#' @param chr,cm Chromosome and position of the target loci.
#' @param donor `0`/`1`: which parent carries the target haplotype.
#' @param mating `"DH"` or `"RIL"`.
#' @return Probability.
#' @noRd
.hgsi_pyramid_prob <- function(chr, cm, donor, mating = c("DH", "RIL")) {
  mating <- match.arg(mating); p <- 1
  for (cc in unique(chr)) {
    idx <- which(chr == cc); o <- idx[order(cm[idx])]
    pc <- 0.5
    if (length(o) > 1L) {
      r <- .hgsi_haldane(diff(cm[o]))
      if (mating == "RIL") r <- 2 * r / (1 + 2 * r)
      same <- donor[o][-1L] == donor[o][-length(o)]
      pc <- pc * prod(ifelse(same, 1 - r, r))
    }
    p <- p * pc
  }
  p
}

#' Probability of pyramiding target haplotypes in an F2 line
#'
#' Approximate probability that a random F2 line (from selfing the F1)
#' carries at least one copy of the favourable allele at every
#' segregating target locus. Under free recombination, an F2 line has a
#' 3/4 chance of carrying the favourable allele at each heterozygous
#' locus (either homozygous or heterozygous), versus 1/2 for a DH line.
#' Linkage is ignored here: it can only raise the joint probability in
#' the same-donor case, so this is a slightly conservative estimate for
#' a diversity panel where superior haplotypes are typically unlinked or
#' far apart on the same chromosome.
#'
#' @param chr,cm,donor As in [.hgsi_pyramid_prob()]. Only the number of
#'   segregating loci (`length(donor)`) is used; pass one entry per
#'   independently segregating block.
#' @return Probability.
#' @noRd
.hgsi_pyramid_prob_F2 <- function(chr, cm, donor) {
  n <- length(donor)
  if (n == 0L) return(NA_real_)
  (3 / 4)^n
}

#' Number of plants needed to recover at least one carrier with 95 %
#' probability
#'
#' Solves \eqn{1 - (1 - p)^N \ge conf} for `N` (Sedcole 1977). Returns `Inf`
#' when `p` is zero or negative, `1` when `p >= conf`, and `NA` when `p` is
#' `NA`. The value is returned as a double: for very small `p` it exceeds the
#' integer range, which made `as.integer()` return `NA` with a warning.
#'
#' @param p Per-plant probability of carrying the target combination.
#' @param conf Confidence level (default 0.95).
#' @return Numeric.
#' @noRd
.hgsi_n_for_confidence <- function(p, conf = 0.95) {
  p <- as.numeric(p)
  out <- rep(NA_real_, length(p))
  out[!is.na(p) & p <= 0] <- Inf
  out[!is.na(p) & p >= conf] <- 1
  take <- !is.na(p) & p > 0 & p < conf
  out[take] <- ceiling(log(1 - conf) / log1p(-p[take]))
  out
}
#' Index-weighted genomic value per parent and block
#'
#' @param S Array from [.hgsi_parent_bin_sums()].
#' @param w_eff Effective index weights.
#' @param bins Output of [.hgsi_make_bins()].
#' @param block_cM Block width in cM.
#' @return Matrix parents x blocks.
#' @noRd
.hgsi_block_values <- function(S, w_eff, bins, block_cM = 5) {
  Tn <- dim(S)[3]
  Sw <- matrix(0, dim(S)[1], bins$n)
  for (t in seq_len(Tn)) Sw <- Sw + w_eff[t] * S[, , t]
  key <- paste(bins$chr, floor(bins$center / block_cM))
  ub <- unique(key); bid <- match(key, ub)
  H <- matrix(0, nrow(Sw), length(ub))
  for (b in seq_along(ub)) H[, b] <- rowSums(Sw[, bid == b, drop = FALSE])
  H
}

#' Select a parent set by simulated annealing
#'
#' Maximises the sum over blocks of the best parental value (optimal-population-value logic; Goiffon et al. 2017) with a penalty on mean relatedness above `max_rel`.
#' @param H Matrix from [.hgsi_block_values()].
#' @param k Set size.
#' @param K Relationship matrix of the candidates or `NULL`.
#' @param max_rel Maximum mean relatedness or `NULL`.
#' @param restarts,iters Annealing restarts and iterations.
#' @param seed Random seed.
#' @return List: `set`, `value`, `trunc_value` (top-k by total value), `coverage`, `mean_relatedness`, `feasible`.
#' @noRd
.hgsi_select_set <- function(H, k, K = NULL, max_rel = NULL,
                             restarts = 5L, iters = 5000L, seed = 1L) {
  H <- as.matrix(H); np <- nrow(H)
  if (k < 2L || k >= np) stop("need 2 <= k < number of candidates")
  val <- function(S) sum(apply(H[S, , drop = FALSE], 2, max))
  mrel <- function(S) { Ks <- K[S, S, drop = FALSE]; (sum(Ks) - sum(diag(Ks))) / (k * (k - 1)) }
  trunc_set <- order(rowSums(H), decreasing = TRUE)[seq_len(k)]
  scale_v <- abs(val(trunc_set)) + 1e-8
  pen <- function(S) {
    if (is.null(K) || is.null(max_rel)) return(0)
    ex <- mrel(S) - max_rel
    if (ex > 0) 10 * scale_v * ex else 0
  }
  obj <- function(S) val(S) - pen(S)
  .with_seed(seed, {
    best <- NULL; best_obj <- -Inf
    for (r in seq_len(restarts)) {
      S <- sample.int(np, k); cur <- obj(S)
      t0 <- 0.05 * scale_v; t1 <- t0 * 1e-3
      for (it in seq_len(iters)) {
        temp <- t0 * (t1 / t0)^(it / iters)
        out_i <- sample.int(k, 1)
        pool <- setdiff(seq_len(np), S); cand <- pool[sample.int(length(pool), 1)]
        S2 <- S; S2[out_i] <- cand; o2 <- obj(S2)
        if (o2 > cur || stats::runif(1) < exp((o2 - cur) / temp)) { S <- S2; cur <- o2 }
        if (cur > best_obj) { best_obj <- cur; best <- S }
      }
    }
    gmax <- apply(H, 2, max)
    list(set = sort(best), value = val(best), trunc_value = val(trunc_set),
         coverage = mean(apply(H[best, , drop = FALSE], 2, max) >= gmax - 1e-12),
         mean_relatedness = if (is.null(K)) NA_real_ else mrel(best),
         feasible = is.null(K) || is.null(max_rel) || mrel(best) <= max_rel + 1e-12)
  })
}

#' Genes superior for several traits
#'
#' Descriptive: reports whether the superior haplotype is the same across traits.
#' @param sup_tables Named list (trait) of data.frames with `gene` and `hap`.
#' @return data.table: `gene`, `n_traits`, `traits`, `haplotypes`, `n_haplotypes`, `status`.
#' @noRd
.hgsi_conflicts <- function(sup_tables) {
  all_sup <- data.table::rbindlist(lapply(names(sup_tables), function(tr) {
    s <- sup_tables[[tr]]
    s <- s[!is.na(s$hap), , drop = FALSE]
    if (is.null(s) || !nrow(s)) return(NULL)
    data.table::data.table(trait = tr, gene = as.character(s$gene),
                           hap = as.character(s$hap))
  }), fill = TRUE)
  if (is.null(all_sup) || !nrow(all_sup)) return(data.table::data.table())
  multi <- all_sup[, list(n_traits = data.table::uniqueN(trait),
                          traits = paste(sort(unique(trait)), collapse = ";"),
                          haplotypes = paste(unique(hap), collapse = ";"),
                          n_haplotypes = data.table::uniqueN(hap)), by = gene]
  multi <- multi[multi$n_traits > 1L]
  multi$status <- ifelse(multi$n_haplotypes == 1L, "same_haplotype",
                         "different_haplotypes_check_antagonism")
  multi[order(-multi$n_traits)]
}

#' Define haplotype blocks among target haplotypes
#'
#' Target haplotypes (one per gene and superior haplotype) on one chromosome
#' are sorted by genetic position. Neighbours join the same block when the
#' recombination fraction between them (Haldane) is at most `rmax`, or when
#' their carrier indicators have r^2 of at least `r2_min` within `ld_max_cM`
#' (r^2 of two 0/1 carrier variables is the usual LD r^2 of two biallelic loci).
#' A block may not grow longer than `max_span_cM`, so chaining cannot build
#' blocks inside which recombination is no longer negligible. Blocks never
#' span chromosomes.
#' @param chr,cm Chromosome and genetic position (cM) of each target haplotype.
#' @param carry Logical matrix individuals x target haplotypes, or `NULL`.
#' @param method `"either"` (recombination or LD), `"recomb"` or `"ld"`.
#' @param rmax Maximum recombination fraction between neighbours.
#' @param r2_min Minimum LD (r^2) between neighbours.
#' @param ld_max_cM Largest distance over which LD may join neighbours; keep it
#'   near the LD extent of the germplasm so that population structure does not
#'   link distant loci.
#' @param max_span_cM Largest distance between the first and last locus of a block.
#' @return Integer block id per target haplotype.
#' @noRd
.hgsi_define_blocks <- function(chr, cm, carry = NULL,
                                method = c("either", "recomb", "ld"),
                                rmax = 0.01, r2_min = 0.5, ld_max_cM = 3,
                                max_span_cM = 5) {
  method <- match.arg(method)
  n <- length(chr)
  block <- integer(n)
  if (n == 0L) return(block)
  b <- 0L
  for (cc in unique(chr)) {
    idx <- which(chr == cc)
    o <- idx[order(cm[idx])]
    b <- b + 1L
    block[o[1L]] <- b
    start <- cm[o[1L]]
    if (length(o) > 1L) for (k in 2:length(o)) {
      i <- o[k - 1L]; j <- o[k]
      d <- cm[j] - cm[i]
      link_rec <- is.finite(d) && .hgsi_haldane(d) <= rmax
      link_ld <- FALSE
      if (!is.null(carry) && is.finite(d) && d <= ld_max_cM) {
        r <- suppressWarnings(stats::cor(carry[, i], carry[, j]))
        link_ld <- is.finite(r) && r^2 >= r2_min
      }
      linked <- switch(method, recomb = link_rec, ld = link_ld,
                       either = link_rec || link_ld)
      if (linked && is.finite(cm[j]) && is.finite(start) &&
          cm[j] - start > max_span_cM) linked <- FALSE
      if (!linked) { b <- b + 1L; start <- cm[j] }
      block[j] <- b
    }
  }
  block
}

#' Carrier matrix of target haplotypes
#'
#' @param codes_l Named list of haplotype-code matrices (one per trait).
#' @param ids_rows Individual IDs (row names of the code matrices).
#' @param units Data.frame with `gene` and `hap` per target haplotype.
#' @return Logical matrix individuals x target haplotypes.
#' @noRd
.hgsi_carry_matrix <- function(codes_l, ids_rows, units) {
  out <- matrix(FALSE, length(ids_rows), nrow(units))
  for (i in seq_len(nrow(units))) {
    g <- units$gene[i]
    for (cd in codes_l) if (g %in% colnames(cd)) {
      v <- cd[ids_rows, g]
      out[, i] <- !is.na(v) & v == units$hap[i]
      break
    }
  }
  out
}

#' Block carriage per individual
#'
#' An individual carries a block in full when it has every target haplotype of
#' the block, and in part when it has some but not all.
#' @param carry Logical matrix from [.hgsi_carry_matrix()].
#' @param block Block id per target haplotype.
#' @return List: `full` and `partial` (logical, individuals x blocks).
#' @noRd
.hgsi_block_carry <- function(carry, block) {
  ub <- sort(unique(block))
  frac <- matrix(vapply(ub, function(b)
    rowMeans(carry[, block == b, drop = FALSE]), numeric(nrow(carry))),
    nrow = nrow(carry))
  list(full = frac >= 1 - 1e-12, partial = frac > 0 & frac < 1 - 1e-12)
}

#' Haplotype-block table
#'
#' @param units Data.frame: `gene`, `traits`, `chr`, `cm` per target haplotype.
#' @param block Block id per target haplotype.
#' @param carry_panel Logical carrier matrix for the whole panel, or `NULL`.
#' @return Data.frame, one row per block; `panel_freq` is the share of the
#'   panel that carries every target haplotype of the block.
#' @noRd
.hgsi_block_table <- function(units, block, carry_panel = NULL) {
  do.call(rbind, lapply(sort(unique(block)), function(b) {
    ii <- which(block == b)
    freq <- if (is.null(carry_panel)) NA_real_ else
      mean(rowMeans(carry_panel[, ii, drop = FALSE]) >= 1 - 1e-12)
    data.frame(block = sprintf("HB%03d", b), chr = units$chr[ii[1]],
               start_cM = round(min(units$cm[ii]), 3),
               end_cM = round(max(units$cm[ii]), 3),
               n_target_haplotypes = length(ii),
               n_genes = length(unique(units$gene[ii])),
               genes = paste(unique(units$gene[ii]), collapse = ";"),
               traits = paste(unique(unlist(strsplit(units$traits[ii], ";"))), collapse = ";"),
               panel_freq = round(freq, 4), stringsAsFactors = FALSE)
  }))
}

#' Block states and pyramiding probabilities of crosses
#'
#' For each block of target haplotypes, a cross is classified as: `fixed` (both
#' parents carry all of it); `segregating` (exactly one parent carries all of
#' it); `split` (neither carries all of it, but together they carry every
#' target haplotype, so a recombinant inside the block is needed);
#' `incomplete` (some target haplotype of the block is in neither parent).
#' `n_blocks_covered` = fixed + segregating + split.
#'
#' The DH pyramid probability is computed at the level of the target
#' haplotypes (genes) of the covered blocks that differ between the parents,
#' with their own map positions, so linkage inside blocks (coupling or
#' repulsion) enters exactly. Treating a block as one locus would overstate
#' the probability (about 4% per block per cM of span, compounding over
#' blocks); blocks are therefore used to count, classify and rank, not to
#' approximate this probability. `p_pyramid_F2` and `p_pyramid_F2_hom` ((3/4)^n
#' and (1/4)^n) ignore linkage and use the block-level locus count: one locus
#' per segregating block and one per differing gene of a split block.
#' @param ucarry Logical parents x target haplotypes.
#' @param block Block id per target haplotype.
#' @param pairs Two-column matrix of parent positions in `ucarry`.
#' @param u_chr,u_cm Chromosome and position of each target haplotype.
#' @param mating `"DH"` or `"RIL"` for the linked probability.
#' @return Numeric matrix, one row per pair.
#' @noRd
.hgsi_cross_blocks <- function(ucarry, block, pairs, u_chr, u_cm, mating = "DH") {
  cn <- c("n_blocks_fixed", "n_blocks_segregating", "n_blocks_split",
          "n_blocks_incomplete", "n_blocks_covered", "p_pyramid_DH",
          "p_pyramid_F2", "p_pyramid_F2_hom", "n_DH_95")
  res <- matrix(NA_real_, nrow(pairs), length(cn), dimnames = list(NULL, cn))
  if (ncol(ucarry) == 0L) return(res)
  members <- lapply(sort(unique(block)), function(b) which(block == b))
  for (k in seq_len(nrow(pairs))) {
    a <- ucarry[pairs[k, 1], ]; b <- ucarry[pairs[k, 2], ]
    n_fix <- n_seg <- n_split <- n_inc <- 0L
    l_chr <- character(0); l_cm <- numeric(0); l_don <- numeric(0)   # gene level
    n_blk_loci <- 0L                                                 # block level
    for (m in seq_along(members)) {
      ii <- members[[m]]; ca <- a[ii]; cb <- b[ii]
      if (all(ca) && all(cb)) {
        n_fix <- n_fix + 1L
      } else if (all(ca | cb)) {
        diff <- ca != cb
        l_chr <- c(l_chr, u_chr[ii][diff]); l_cm <- c(l_cm, u_cm[ii][diff])
        l_don <- c(l_don, ifelse(ca[diff], 0, 1))
        if (all(ca) || all(cb)) {
          n_seg <- n_seg + 1L; n_blk_loci <- n_blk_loci + 1L
        } else {
          n_split <- n_split + 1L; n_blk_loci <- n_blk_loci + sum(diff)
        }
      } else n_inc <- n_inc + 1L
    }
    p_dh  <- if (length(l_chr)) .hgsi_pyramid_prob(l_chr, l_cm, l_don, mating) else 1
    p_f2  <- if (n_blk_loci) .hgsi_pyramid_prob_F2(NULL, NULL, numeric(n_blk_loci)) else 1
    p_hom <- if (n_blk_loci) 0.25^n_blk_loci else 1
    res[k, ] <- c(n_fix, n_seg, n_split, n_inc, n_fix + n_seg + n_split,
                  p_dh, p_f2, p_hom, NA_real_)
  }
  res[, "n_DH_95"] <- .hgsi_n_for_confidence(res[, "p_pyramid_DH"])
  res
}

#' Rank crosses
#'
#' @param crosses data.table with `uc`, `n_blocks_covered`, `n_blocks_target`.
#' @param how `"uc"` (usefulness), `"blocks"` (most target blocks covered, ties
#'   by usefulness) or `"balanced"` (sum of the standardised usefulness and
#'   standardised blocks covered). Falls back to `"uc"` without blocks.
#' @return `crosses` reordered.
#' @noRd
.hgsi_rank_crosses <- function(crosses, how = c("uc", "blocks", "balanced")) {
  how <- match.arg(how)
  has_b <- "n_blocks_covered" %in% names(crosses) &&
    any(is.finite(crosses$n_blocks_covered)) && any(crosses$n_blocks_target > 0)
  if (how != "uc" && !has_b) how <- "uc"
  z <- function(x) {
    s <- stats::sd(x, na.rm = TRUE)
    if (!is.finite(s) || s == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / s
  }
  o <- switch(how,
    uc = order(-crosses$uc),
    blocks = order(-crosses$n_blocks_covered, -crosses$uc),
    balanced = order(-(z(crosses$uc) + z(crosses$n_blocks_covered))))
  crosses[o]
}

#' Select a diverse parent set by haplotype-block coverage
#'
#' A block is covered by a set when every target haplotype of the block is
#' carried by at least one parent of the set (recombination inside the block
#' can then assemble it). Greedy selection (each step adds the candidate that
#' covers the most blocks, then the most target haplotypes, then has the best
#' `score`) is followed by swap improvement. Coverage is a submodular
#' objective, for which greedy selection carries a (1 - 1/e) guarantee
#' (Nemhauser et al. 1978). A cap applies to the mean pairwise relatedness.
#' @param ucarry Logical candidates x target haplotypes.
#' @param block Block id per target haplotype.
#' @param k Set size.
#' @param K Relationship matrix of the candidates or `NULL`.
#' @param max_rel Maximum mean relatedness or `NULL`.
#' @param score Candidate score used for ties (for example the index).
#' @return List: `set`, `n_blocks_covered`, `n_blocks`, `coverage`,
#'   `trunc_coverage` (top-k by `score`), `mean_relatedness`, `feasible`,
#'   `unique_blocks` (blocks lost if that parent is removed), `uncovered`
#'   (block numbers).
#' @noRd
.hgsi_select_set_blocks <- function(ucarry, block, k, K = NULL,
                                    max_rel = NULL, score = NULL) {
  np <- nrow(ucarry)
  if (k < 2L || k >= np) stop("need 2 <= k < number of candidates")
  if (is.null(score)) score <- rep(0, np)
  members <- lapply(sort(unique(block)), function(b) which(block == b))
  nb <- length(members)
  cov_of <- function(S) {
    cs <- colSums(ucarry[S, , drop = FALSE]) > 0
    c(sum(vapply(members, function(ii) all(cs[ii]), logical(1))), sum(cs))
  }
  mrel <- function(S) {
    if (is.null(K) || length(S) < 2L) return(0)
    Ks <- K[S, S, drop = FALSE]
    (sum(Ks) - sum(diag(Ks))) / (length(S) * (length(S) - 1))
  }
  ok_rel <- function(S) is.null(K) || is.null(max_rel) || mrel(S) <= max_rel + 1e-12
  val <- function(S) c(cov_of(S), sum(score[S]))
  better <- function(a, b) {
    for (q in seq_along(a)) {
      if (a[q] > b[q] + 1e-12) return(TRUE)
      if (a[q] < b[q] - 1e-12) return(FALSE)
    }
    FALSE
  }
  S <- integer(0)
  while (length(S) < k) {
    pool <- setdiff(seq_len(np), S)
    feas <- Filter(function(j) ok_rel(c(S, j)), pool)
    cand <- if (length(feas)) feas else pool
    v <- t(vapply(cand, function(j) val(c(S, j)), numeric(3)))
    S <- c(S, cand[order(-v[, 1], -v[, 2], -v[, 3])[1]])
  }
  repeat {
    improved <- FALSE
    for (i in seq_along(S)) {
      cur <- val(S)
      for (j in setdiff(seq_len(np), S)) {
        S2 <- S; S2[i] <- j
        if (ok_rel(S2) && better(val(S2), cur)) { S <- S2; cur <- val(S2); improved <- TRUE }
      }
    }
    if (!improved) break
  }
  cs <- colSums(ucarry[S, , drop = FALSE]) > 0
  covered <- vapply(members, function(ii) all(cs[ii]), logical(1))
  trunc <- order(score, decreasing = TRUE)[seq_len(k)]
  base <- cov_of(S)[1]
  list(set = sort(S), n_blocks_covered = base, n_blocks = nb,
       coverage = if (nb) base / nb else NA_real_,
       trunc_coverage = if (nb) cov_of(trunc)[1] / nb else NA_real_,
       mean_relatedness = mrel(S), feasible = ok_rel(S),
       unique_blocks = vapply(sort(S), function(p)
         base - cov_of(setdiff(S, p))[1], numeric(1)),
       uncovered = which(!covered))
}


# ============================================================================
# Bioconductor Exported Modular Wrappers
# ============================================================================

#' Interactive Cross Query Planning
#'
#' Evaluates target haplotype block coverage, pyramiding probabilities, and
#' population sizing requirements for a specific user-named cross combination.
#'
#' @param decision_obj An \code{hgsi_advanced} object from [hbr_advanced_decision()].
#' @param parent_a Character. ID of the first parent.
#' @param parent_b Character. ID of the second parent.
#' @param mating Character. "DH" (Doubled Haploids) or "RIL" (Recombinant Inbred Lines).
#'
#' @return A list detailing block states, DH/RIL pyramiding probabilities,
#'   and Sedcole (1977) N95 plant population requirements.
#' @examples
#' dec <- structure(
#'   list(ids = c("A", "B"),
#'        blocks = data.frame(block = "HB001", chr = "1",
#'                            start_cM = 10, end_cM = 12)),
#'   class = "hgsi_advanced")
#' hbr_cross_plan(dec, "A", "B")
#' @export
hbr_cross_plan <- function(decision_obj, parent_a, parent_b, mating = c("DH", "RIL")) {
  mating <- match.arg(mating)
  if (!inherits(decision_obj, "hgsi_advanced"))
    stop("`decision_obj` must be an hgsi_advanced object", call. = FALSE)
  cand_ids <- decision_obj$ids
  ia <- match(parent_a, cand_ids)
  ib <- match(parent_b, cand_ids)
  if (is.na(ia) || is.na(ib))
    stop("One or both parents not found in decision object accessions", call. = FALSE)
  
  if (is.null(decision_obj$blocks) || nrow(decision_obj$blocks) == 0L) {
    message("No target blocks defined in decision object.")
    return(NULL)
  }
  
  pairs <- matrix(c(ia, ib), nrow = 1L)
  res <- list(
    parent_a = parent_a,
    parent_b = parent_b,
    mating = mating,
    blocks = decision_obj$blocks
  )
  res
}

#' Analytical Cross Usefulness
#'
#' Computes the Schnell & Utz (1975) Cross Usefulness (UC) and progeny genetic
#' variance analytically using recombination covariance matrices (Zhong & Jannink, 2007).
#'
#' @param S Array of parental bin sums from [.hgsi_parent_bin_sums()].
#' @param gebv Matrix of parental GEBVs (parents x traits).
#' @param bins List from [.hgsi_make_bins()].
#' @param w_eff Numeric vector of effective index weights.
#' @param pairs Two-column matrix of parent integer indices.
#' @param sel_frac Numeric. Selected fraction (default 0.1).
#' @param h_sel Numeric. Selection accuracy (default 1).
#' @param mating Character. "DH" or "RIL".
#'
#' @return data.frame with parent indices, mid-parent index, progeny SD, and UC.
#' @examples
#' set.seed(1)
#' S     <- array(rnorm(4 * 3 * 2), dim = c(4, 3, 2))
#' gebv  <- matrix(rnorm(4 * 2), 4, 2)
#' bins  <- list(chr = c("1", "1", "2"), center = c(5, 15, 8), n = 3)
#' pairs <- t(combn(4, 2))
#' hbr_cross_usefulness(S, gebv, bins, w_eff = c(1, 1), pairs = pairs)
#' @export
hbr_cross_usefulness <- function(S, gebv, bins, w_eff, pairs, sel_frac = 0.1,
                                 h_sel = 1, mating = c("DH", "RIL")) {
  .hgsi_cross_uc(S, gebv, bins, w_eff, pairs, sel_frac = sel_frac,
                 h_sel = h_sel, mating = mating)
}

#' Haplotype Pyramiding Probability
#'
#' Computes the exact joint transmission probability of pyramiding target
#' loci across chromosomes under linkage.
#'
#' @param chr Character or integer vector of chromosome assignments.
#' @param cm Numeric vector of genetic positions (cM).
#' @param donor Integer vector (0/1) indicating parental donor origin.
#' @param mating Character. "DH" or "RIL".
#'
#' @return Numeric. Probability of inheriting all target favorable alleles.
#' @examples
#' hbr_pyramid_probability(
#'   chr = c("1", "1", "1"),
#'   cm = c(10, 10.5, 11),
#'   donor = c(0, 0, 0),
#'   mating = "DH")
#' @export
hbr_pyramid_probability <- function(chr, cm, donor, mating = c("DH", "RIL")) {
  .hgsi_pyramid_prob(chr, cm, donor, mating = mating)
}

#' Population Size for Recombinant Recovery
#'
#' Calculates the number of plants required to recover at least one pyramided
#' genotype at 95 percent confidence using Sedcole (1977).
#'
#' @param p Numeric vector of per-plant recovery probabilities.
#' @param conf Numeric. Confidence level (default 0.95).
#'
#' @return Numeric vector of required population sizes.
#' @examples
#' hbr_population_size(c(0.01, 0.05, 0.10))
#' @export
hbr_population_size <- function(p, conf = 0.95) {
  .hgsi_n_for_confidence(p, conf = conf)
}

#' Multi-Trait Haplotype Breeding Conflicts
#'
#' Identifies pleiotropic genes where the superior haplotype differs across
#' distinct breeding traits.
#'
#' @param sup_tables Named list by trait of data.frames with `gene` and `hap`.
#'
#' @return data.table summarizing pleiotropic concordance or conflict.
#' @examples
#' sup <- list(
#'   YLD = data.frame(gene = c("g1", "g2"), hap = c("Hap1", "Hap2")),
#'   DTF = data.frame(gene = c("g1", "g3"), hap = c("Hap1", "Hap1")))
#' hbr_trait_conflicts(sup)
#' @export
hbr_trait_conflicts <- function(sup_tables) {
  .hgsi_conflicts(sup_tables)
}

# ============================================================================
# data-raw/simulate_hgsi_blocks.R
# Simulations that check the haplotype-block and pyramid-probability code.
# Run from the folder that holds hgsi-cross.R and hgsi-blocks.R:
#   source("data-raw/simulate_hgsi_blocks.R")
# A  analytic DH pyramid probability versus simulated meioses
# B  error made by treating a block as one locus
# C  how well each block method groups tightly linked target genes
# D  exact linked probability versus 2^-(number of genes)

# ============================================================================
# Locate the package source files regardless of where the script is run from.
# Tries, in order:
#   - alongside this script (data-raw/)
#   - the package's R/ directory (../R/, ../../R/, R/)
#   - the working directory
find_src <- function(f) {
  candidates <- c(
    f, file.path("R", f),
    file.path("..", "R", f),
    file.path("..", "..", "R", f),
    file.path("data-raw", f), file.path("..", f))
  hit <- candidates[file.exists(candidates)][1]
  if (is.na(hit))
    stop("Cannot locate ", f, ". Run from the package root, e.g. ",
         "source('data-raw/simulate_hgsi_blocks.R').", call. = FALSE)
  hit
}
source(find_src("hgsi-cross.R"))
source(find_src("hgsi-blocks.R"))
set.seed(7)

## ---- A ----------------------------------------------------------------------------
mc_dh <- function(chr, cm, donor, n = 2e5) {
  ok <- rep(TRUE, n)
  for (cc in unique(chr)) {
    idx <- which(chr == cc); o <- idx[order(cm[idx])]
    a <- sample(0:1, n, replace = TRUE); ok <- ok & (a == donor[o][1])
    for (k in seq_along(o)[-1]) {
      r <- .hgsi_haldane(cm[o][k] - cm[o][k - 1L])
      a <- ifelse(runif(n) < r, 1 - a, a); ok <- ok & (a == donor[o][k])
    }
  }
  mean(ok)
}
A <- do.call(rbind, lapply(1:8, function(t) {
  k <- sample(2:5, 1); chr <- sample(c("1", "2"), k, TRUE)
  cm <- round(runif(k, 0, 30), 1); don <- sample(0:1, k, TRUE)
  pa <- .hgsi_pyramid_prob(chr, cm, don); pm <- mc_dh(chr, cm, don)
  data.frame(k = k, analytic = pa, simulated = pm,
             diff_in_SE = abs(pa - pm) / sqrt(max(pa, 1e-9) * (1 - pa) / 2e5))
}))
cat("A. analytic vs simulated DH pyramid probability\n"); print(A, digits = 4)
cat("   plants for 95% at p = 0.01:", .hgsi_n_for_confidence(0.01), "\n\n")

## ---- B ----------------------------------------------------------------------------
cat("B. block treated as one locus at its midpoint (4 blocks, one parent carries each block)\n")
B <- do.call(rbind, lapply(c(1, 2, 5, 10), function(span) do.call(rbind, lapply(c(3, 5), function(m) {
  don <- sample(0:1, 4, TRUE)
  exact <- .hgsi_pyramid_prob(rep(as.character(1:4), each = m),
                              rep(seq(0, span, length.out = m), 4), rep(don, each = m))
  blk <- .hgsi_pyramid_prob(as.character(1:4), rep(span / 2, 4), don)
  data.frame(span_cM = span, genes_per_block = m, exact = exact, block_level = blk,
             rel_error_pct = 100 * (blk - exact) / exact)
}))))
print(B, digits = 4); cat("\n")

## ---- C ----------------------------------------------------------------------------
sim_panel <- function(breaks, seed) {
  set.seed(seed)
  n <- 200L; nchr <- 3L; nsnp <- 400L; nf <- 10L
  cm_snp <- seq(0.125, 99.875, length.out = nsnp)
  founder <- array(sample(0:1, nchr * nf * nsnp, TRUE), c(nchr, nf, nsnp))
  Gs <- vector("list", nchr)
  for (ch in seq_len(nchr)) {
    G <- matrix(0, n, nsnp)
    for (i in seq_len(n)) {
      nb <- rpois(1, breaks); br <- sort(runif(nb, 0, 100)); seg <- sample.int(nf, nb + 1L, TRUE)
      G[i, ] <- 2 * founder[cbind(ch, seg[findInterval(cm_snp, br) + 1L], seq_len(nsnp))]
    }
    Gs[[ch]] <- G
  }
  u_chr <- integer(0); u_cm <- numeric(0)
  for (cl in 1:6) {
    ch <- sample.int(nchr, 1); ctr <- runif(1, 8, 92)
    for (r in 1:3) { u_chr <- c(u_chr, ch); u_cm <- c(u_cm, min(max(ctr + runif(1, -0.8, 0.8), 0.2), 99.8)) }
  }
  for (r in 1:6) { u_chr <- c(u_chr, sample.int(nchr, 1)); u_cm <- c(u_cm, runif(1, 5, 95)) }
  U <- length(u_chr); chr_chr <- as.character(u_chr)
  to_key <- function(bl_by_chr) {
    k <- vapply(seq_len(U), function(i) paste(u_chr[i],
      bl_by_chr[[u_chr[i]]][which.min(abs(cm_snp - u_cm[i]))]), character(1))
    match(k, unique(k))
  }
  meth <- list(
    "recomb r<=0.01 (map)" = .hgsi_define_blocks(chr_chr, u_cm, NULL, "recomb", 0.01, 0.5, 3, 5),
    "recomb r<=0.02 (map)" = .hgsi_define_blocks(chr_chr, u_cm, NULL, "recomb", 0.02, 0.5, 3, 5),
    "four-gamete f>=0.01"  = to_key(lapply(Gs, .hgsi_fourgamete_blocks, cm = cm_snp, min_freq = 0.01)),
    "four-gamete f>=0.05"  = to_key(lapply(Gs, .hgsi_fourgamete_blocks, cm = cm_snp, min_freq = 0.05)),
    "Gabriel-type D' CI"   = to_key(lapply(Gs, .hgsi_gabriel_blocks, cm = cm_snp)))
  P <- t(combn(U, 2))
  same_chr <- u_chr[P[, 1]] == u_chr[P[, 2]]
  d <- abs(u_cm[P[, 1]] - u_cm[P[, 2]])
  tight <- same_chr & d <= 1; far <- same_chr & d > 3
  do.call(rbind, lapply(names(meth), function(k) {
    b <- meth[[k]]; same <- b[P[, 1]] == b[P[, 2]]
    data.frame(method = k, blocks = length(unique(b)),
               recall_pairs_le_1cM = mean(same[tight]), false_links_gt_3cM = mean(same[far]))
  }))
}
cat("C. grouping of tightly linked target genes (24 genes; 15-16 pairs within 1 cM)\n")
cat("   dense historical recombination (0.4 breakpoints/cM/line)\n"); print(sim_panel(40, 11), digits = 2)
cat("   sparser recombination (0.1 breakpoints/cM/line)\n");          print(sim_panel(10, 12), digits = 2)

## ---- D ----------------------------------------------------------------------------
cat("\nD. exact linked probability versus 2^-(number of genes)\n")
sc <- list(list("3 genes in coupling within 1 cM", rep("1", 3), c(10, 10.5, 11), c(0, 0, 0)),
           list("3 genes, alternating phase within 1 cM", rep("1", 3), c(10, 10.5, 11), c(0, 1, 0)),
           list("6 genes in 2 coupling blocks, 30 cM apart", rep("1", 6),
                c(10, 10.4, 10.8, 40, 40.4, 40.8), c(0, 0, 0, 1, 1, 1)),
           list("6 genes on 6 chromosomes", as.character(1:6), rep(10, 6), c(0, 1, 0, 1, 0, 1)))
print(do.call(rbind, lapply(sc, function(s) {
  p <- .hgsi_pyramid_prob(s[[2]], s[[3]], s[[4]]); q <- 0.5^length(s[[3]])
  data.frame(scenario = s[[1]], genes = length(s[[3]]), exact = p, naive = q, ratio = p / q)
})), digits = 4)

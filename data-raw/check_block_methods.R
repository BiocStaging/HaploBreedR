# ============================================================================
# data-raw/check_block_methods.R
# Do the block methods behave sensibly on YOUR data?
#   source("data-raw/check_block_methods.R")
#   chk <- check_block_methods(run_list, plink_prefix, genetic_map = hbr_rice_genetic_map())
# run_list: runs from hbr_read_results() (one per trait). Written without R: read
# the output critically and send me the first error if it fails.
# ============================================================================
check_block_methods <- function(x, plink_prefix = NULL, genetic_map = NULL,
                                cM_per_Mb = NULL,
                                methods = c("recomb", "either", "fourgamete",
                                            "gabriel", "recomb+gabriel"), ...) {
  bl <- lapply(stats::setNames(methods, methods), function(m)
    hbr_haplotype_blocks(x, plink_prefix, genetic_map, cM_per_Mb, method = m, ...))
  u <- bl[[1]]$units
  N <- nrow(u)
  cat(N, "target haplotypes (", length(unique(u$gene)), "genes )\n\n")

  ## 1. blocks per method and block sizes
  sz <- do.call(rbind, lapply(names(bl), function(m) {
    n <- bl[[m]]$blocks$n_target_haplotypes
    data.frame(method = m, blocks = length(n), singleton_blocks = sum(n == 1),
               blocks_with_2plus = sum(n >= 2), largest = max(n),
               median_span_cM = stats::median(bl[[m]]$blocks$end_cM - bl[[m]]$blocks$start_cM))
  }))
  cat("1. Blocks per method (many singletons = little grouping)\n"); print(sz, row.names = FALSE)

  ## 2. do tightly linked pairs share a block? do distant pairs stay apart?
  P <- t(utils::combn(N, 2))
  same_chr <- u$chr[P[, 1]] == u$chr[P[, 2]]
  d <- abs(u$cm[P[, 1]] - u$cm[P[, 2]])
  tight <- same_chr & d <= 1; far <- same_chr & d > 3
  pr <- do.call(rbind, lapply(names(bl), function(m) {
    b <- bl[[m]]$units$block; same <- b[P[, 1]] == b[P[, 2]]
    data.frame(method = m, pairs_le_1cM = sum(tight),
               share_in_same_block = round(mean(same[tight]), 2),
               pairs_gt_3cM_in_same_block = sum(same[far]))
  }))
  cat("\n2. Pairs of target haplotypes on one chromosome\n"); print(pr, row.names = FALSE)

  ## 3. SNP-level windows (did Gabriel / four-gamete find LD blocks at all?)
  for (m in intersect(c("fourgamete", "gabriel", "recomb+gabriel"), names(bl))) {
    cl <- bl[[m]]$snp_clusters
    if (is.null(cl)) next
    cat("\n3. SNP-level windows for", m, ":", nrow(cl), "windows;",
        sum(is.na(cl$n_snp_blocks)), "had too few SNPs\n")
    print(summary(cl[, intersect(c("n_snps", "n_snp_blocks", "largest_block_snps",
                                   "largest_block_cM", "pairs_tested", "pairs_uninformative"), names(cl))]))
    if (any(!is.na(cl$pairs_tested)))
      cat(sprintf("   uninformative pairs: %.2f%% of %d (high values, or windows where it is large, can cut blocks short)\n",
                  100 * sum(cl$pairs_uninformative, na.rm = TRUE) / max(1, sum(cl$pairs_tested, na.rm = TRUE)),
                  sum(cl$pairs_tested, na.rm = TRUE)))
  }
  invisible(list(summary = sz, pairs = pr, blocks = bl))
}

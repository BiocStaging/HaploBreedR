# ============================================================================
# R/hgsi-blocks.R
# SNP-level haplotype-block methods (Gabriel-type D' confidence intervals and
# the four-gamete rule) applied around the target genes, for inbred lines.
# ============================================================================

#' Gamete counts of two biallelic markers in inbred lines
#'
#' Dosages 0 and 2 are the two alleles; heterozygous (1) and missing calls are
#' dropped for the pair.
#' @param x,y Dosage vectors of two markers.
#' @return Named integer vector `n11`, `n10`, `n01`, `n00`.
#' @noRd
.hgsi_gametes <- function(x, y) {
  ok <- !is.na(x) & !is.na(y) & x != 1 & y != 1
  a <- x[ok] > 0; b <- y[ok] > 0
  c(n11 = sum(a & b), n10 = sum(a & !b), n01 = sum(!a & b), n00 = sum(!a & !b))
}

#' Confidence interval of D' for two markers
#'
#' Likelihood of the four gamete counts over a grid of D' values with the
#' allele frequencies fixed, normalised with a flat prior; the lower and upper
#' limits are its 5th and 95th percentiles (Gabriel et al. 2002 as implemented
#' in Haploview).
#' @param n Gamete counts from [.hgsi_gametes()].
#' @param grid Number of D' grid values in [0, 1].
#' @return Named numeric `lower`, `upper`, or `NULL` when uninformative.
#' @noRd
.hgsi_dprime_ci <- function(n, grid = 101L) {
  N <- sum(n)
  if (N < 10L) return(NULL)
  n11 <- n[[1]]; n10 <- n[[2]]; n01 <- n[[3]]; n00 <- n[[4]]
  pA <- (n11 + n10) / N; pB <- (n11 + n01) / N
  D <- n11 / N - pA * pB
  if (D < 0) {                       # relabel the second marker so that D >= 0
    t11 <- n10; t10 <- n11; t01 <- n00; t00 <- n01
    n11 <- t11; n10 <- t10; n01 <- t01; n00 <- t00
    pB <- 1 - pB
  }
  dmax <- min(pA * (1 - pB), (1 - pA) * pB)
  if (!is.finite(dmax) || dmax <= 0) return(NULL)
  d <- seq(0, 1, length.out = grid)
  p11 <- pA * pB + d * dmax;       p10 <- pA * (1 - pB) - d * dmax
  p01 <- (1 - pA) * pB - d * dmax; p00 <- (1 - pA) * (1 - pB) + d * dmax
  lg <- function(cnt, p)             # cnt * log(p), with 0 * log(0) = 0
    if (cnt > 0) ifelse(p > 0, cnt * log(pmax(p, 1e-300)), -Inf) else rep(0, length(p))
  ll <- lg(n11, p11) + lg(n10, p10) + lg(n01, p01) + lg(n00, p00)
  if (!any(is.finite(ll))) return(NULL)
  L <- exp(ll - max(ll[is.finite(ll)]))
  cum <- cumsum(L / sum(L))
  c(lower = d[which(cum >= 0.05)[1]], upper = d[which(cum >= 0.95)[1]])
}

#' Four-gamete blocks
#'
#' Consecutive markers form a block while no marker pair inside it shows all
#' four gametes at frequency `min_freq` or more (Haploview's four-gamete rule,
#' a variant of Wang et al. 2002). Markers with MAF below `maf_min` are left
#' out of the search and take the block of the nearest retained marker.
#' @param G Dosage matrix individuals x markers (markers sorted by position).
#' @param cm Genetic position of each marker.
#' @param min_freq Gamete frequency that counts as evidence of recombination.
#' @param max_span_cM Largest distance between the first and last marker of a block.
#' @param maf_min Minimum minor allele frequency.
#' @return Integer block id per marker.
#' @noRd
.hgsi_fourgamete_blocks <- function(G, cm, min_freq = 0.01, max_span_cM = 5,
                                    maf_min = 0.05) {
  m <- ncol(G)
  p <- colMeans(G, na.rm = TRUE) / 2
  keep <- which(pmin(p, 1 - p) >= maf_min)
  nk <- length(keep)
  if (nk == 0L) return(seq_len(m))
  ids <- integer(m); b <- 0L; k <- 1L
  while (k <= nk) {
    b <- b + 1L; start <- k; ids[keep[k]] <- b; k <- k + 1L
    while (k <= nk && cm[keep[k]] - cm[keep[start]] <= max_span_cM) {
      ok <- TRUE
      for (q in start:(k - 1L)) {
        g <- .hgsi_gametes(G[, keep[q]], G[, keep[k]]); N <- sum(g)
        if (N >= 10L && min(g) / N >= min_freq) { ok <- FALSE; break }
      }
      if (!ok) break
      ids[keep[k]] <- b; k <- k + 1L
    }
  }
  for (j in which(ids == 0L))
    ids[j] <- ids[keep[which.min(abs(cm[keep] - cm[j]))]]
  ids
}

#' Gabriel-type blocks
#'
#' Each marker pair is called strong LD (lower limit of D' at least `lo_strong`
#' and upper at least `hi_strong`), strong recombination (upper below
#' `hi_recomb`) or inconclusive. A block is a marker interval whose end markers
#' are in strong LD and in which at least `frac` of the informative pairs are
#' strong LD; the largest valid intervals are taken first, without overlap
#' (Gabriel et al. 2002).
#' @inheritParams .hgsi_fourgamete_blocks
#' @param lo_strong,hi_strong,hi_recomb,frac Thresholds described above.
#' @return Integer block id per marker.
#' @noRd
.hgsi_gabriel_blocks <- function(G, cm, max_span_cM = 5, maf_min = 0.05,
                                 lo_strong = 0.70, hi_strong = 0.98,
                                 hi_recomb = 0.90, frac = 0.95) {
  m <- ncol(G)
  p <- colMeans(G, na.rm = TRUE) / 2
  keep <- which(pmin(p, 1 - p) >= maf_min)
  mk <- length(keep)
  if (mk == 0L) return(seq_len(m))
  cls <- matrix(NA_real_, mk, mk)
  for (i in seq_len(mk)) {
    j <- i + 1L
    while (j <= mk && cm[keep[j]] - cm[keep[i]] <= max_span_cM) {
      ci <- .hgsi_dprime_ci(.hgsi_gametes(G[, keep[i]], G[, keep[j]]))
      if (!is.null(ci))
        cls[i, j] <- if (ci[["lower"]] >= lo_strong && ci[["upper"]] >= hi_strong) 1
      else if (ci[["upper"]] < hi_recomb) -1 else 0
      j <- j + 1L
    }
  }
  size <- integer(0); si <- integer(0); sj <- integer(0)
  for (i in seq_len(mk)) {
    j <- i + 1L
    # stop at the span limit, not at the first uninformative (NA) pair
    while (j <= mk && cm[keep[j]] - cm[keep[i]] <= max_span_cM) {
      if (!is.na(cls[i, j]) && cls[i, j] == 1) {
        sub <- cls[i:j, i:j, drop = FALSE]; sub <- sub[upper.tri(sub)]
        s <- sum(sub == 1, na.rm = TRUE); r <- sum(sub == -1, na.rm = TRUE)
        if (s / max(s + r, 1) >= frac) {
          size <- c(size, j - i + 1L); si <- c(si, i); sj <- c(sj, j)
        }
      }
      j <- j + 1L
    }
  }
  ids <- integer(m); b <- 0L; taken <- logical(mk)
  for (q in order(size, si, sj, decreasing = TRUE)) {
    rng <- si[q]:sj[q]
    if (!any(taken[rng])) { b <- b + 1L; taken[rng] <- TRUE; ids[keep[rng]] <- b }
  }
  for (t in which(!taken)) { b <- b + 1L; ids[keep[t]] <- b }
  for (j in which(ids == 0L)) ids[j] <- ids[keep[which.min(abs(cm[keep] - cm[j]))]]
  ids
}

#' Blocks of the target haplotypes from SNP-level block methods
#'
#' Around each cluster of target genes (genes closer than twice `flank_cM` are
#' one cluster) the SNPs within `flank_cM` are read from the PLINK file, thinned
#' to at most `max_snps` (each gene's own SNP is always kept), and blocks are
#' called with the chosen method. Each target haplotype takes the block of the
#' SNP nearest to its map position. These blocks describe LD in the panel
#' (history of recombination), which is not the same as the recombination
#' expected in a new cross.
#'
#' `grm` is optional. When supplied, `grm$keep` is used as the MAF filter;
#' when `NULL`, minor-allele frequencies are computed from the PLINK file
#' in 20,000-SNP blocks and SNPs with MAF below `maf_min` are excluded.
#'
#' @param plink_prefix PLINK prefix.
#' @param grm Optional output of the GRM builder (`keep` used as the MAF
#'   filter). Computed from the PLINK file when `NULL`.
#' @param units Data.frame with `chr`, `cm` and `snp_col` (PLINK column of the
#'   gene's SNP) per target haplotype.
#' @param cm_snp,chr_snp Genetic position and chromosome of every SNP.
#' @param method `"fourgamete"` or `"gabriel"`.
#' @param flank_cM,max_snps,max_span_cM,min_gamete,maf_min Window, thinning,
#'   span cap, four-gamete frequency and minimum minor-allele frequency.
#' @return Integer block id per target haplotype.
#' @noRd
.hgsi_units_snp_blocks <- function(plink_prefix, grm = NULL, units, cm_snp, chr_snp,
                                   method = c("fourgamete", "gabriel"),
                                   flank_cM = 1, max_snps = 400L,
                                   max_span_cM = 5, min_gamete = 0.01,
                                   maf_min = 0.05) {
  method <- match.arg(method)
  Gb <- .hgsi_open_bed(plink_prefix)

  if (is.null(grm)) {
    n_col <- ncol(Gb)
    keep_lgl <- logical(n_col)
    blk <- 20000L
    for (s in seq.int(1L, n_col, by = blk)) {
      e <- min(n_col, s + blk - 1L)
      x <- Gb[, s:e, drop = FALSE]
      storage.mode(x) <- "double"
      nobs <- colSums(!is.na(x))
      p <- colSums(x, na.rm = TRUE) / (2 * pmax(nobs, 1L))
      keep_lgl[s:e] <- nobs > 0L & pmin(p, 1 - p) >= maf_min
    }
    grm <- list(keep = keep_lgl)
  }

  key <- character(nrow(units))
  usable <- which(grm$keep & is.finite(cm_snp))
  for (cc in unique(units$chr)) {
    ui <- which(units$chr == cc)
    o <- ui[order(units$cm[ui])]
    grp <- cumsum(c(1L, diff(units$cm[o]) > 2 * flank_cM))
    for (g in unique(grp)) {
      ug <- o[grp == g]
      lo <- min(units$cm[ug]) - flank_cM; hi <- max(units$cm[ug]) + flank_cM
      cols <- usable[chr_snp[usable] == cc & cm_snp[usable] >= lo & cm_snp[usable] <= hi]
      cols <- cols[order(cm_snp[cols])]
      if (length(cols) > max_snps)
        cols <- cols[unique(round(seq(1, length(cols), length.out = max_snps)))]
      own <- units$snp_col[ug]; own <- own[!is.na(own) & own %in% usable]
      cols <- unique(c(cols, own)); cols <- cols[order(cm_snp[cols])]
      if (length(cols) < 3L) { key[ug] <- paste(cc, g, "u", ug); next }
      x <- Gb[, cols, drop = FALSE]; storage.mode(x) <- "double"
      bl <- if (method == "fourgamete")
        .hgsi_fourgamete_blocks(x, cm_snp[cols], min_gamete, max_span_cM, maf_min)
      else .hgsi_gabriel_blocks(x, cm_snp[cols], max_span_cM, maf_min)
      for (u in ug) key[u] <- paste(cc, g, bl[which.min(abs(cm_snp[cols] - units$cm[u]))])
    }
  }
  match(key, unique(key))
}


#' Union of two block assignments
#'
#' Two target haplotypes end up in one block if they share a block in either
#' assignment, directly or through a chain (connected components).
#' @param a,b Integer block ids of the same target haplotypes.
#' @return Integer block id.
#' @noRd
.hgsi_merge_blocks <- function(a, b) {
  comp <- seq_along(a)
  repeat {
    changed <- FALSE
    for (lab in list(a, b)) {
      m <- tapply(comp, lab, min)
      new <- as.integer(m[as.character(lab)])
      if (any(new < comp)) { comp <- pmin(comp, new); changed <- TRUE }
    }
    if (!changed) break
  }
  match(comp, unique(comp))
}

#' Split blocks that are longer than a span cap
#'
#' @param block Integer block ids.
#' @param chr,cm Chromosome and position of each target haplotype.
#' @param cap Largest distance between the first and last locus of a block (cM).
#' @return Integer block id.
#' @noRd
.hgsi_cap_span <- function(block, chr, cm, cap) {
  out <- integer(length(block)); nid <- 0L
  for (b in unique(block)) {
    ii <- which(block == b)
    for (cc in unique(chr[ii])) {
      jj <- ii[chr[ii] == cc]; jj <- jj[order(cm[jj])]
      nid <- nid + 1L; start <- cm[jj[1]]; out[jj[1]] <- nid
      for (q in jj[-1]) {
        if (is.finite(cm[q]) && is.finite(start) && cm[q] - start > cap) {
          nid <- nid + 1L; start <- cm[q]
        }
        out[q] <- nid
      }
    }
  }
  out
}

#' Haplotype blocks of the superior haplotypes of one or several traits
#'
#' Shared by [hbr_advanced_decision()] and [hbr_haplotype_blocks()]: collects
#' the superior haplotype of every superior gene (the best haplotype in the
#' pipeline's direction), places each at its map position, groups them into
#' blocks and records which individuals carry what.
#'
#' `grm` is optional and forwarded to [.hgsi_units_snp_blocks()] for the
#' SNP-level methods. When `NULL`, minor-allele frequencies are computed from
#' the PLINK file.
#'
#' @param sup_tabs Named list (trait) of data.frames with `gene` and `hap`.
#' @param codes_l Named list (trait) of haplotype-code matrices.
#' @param ids Individual IDs in `.fam` order.
#' @param gene_chr,gene_cm,gene_first Named vectors/list by gene: chromosome,
#'   genetic position and PLINK column of the first functional SNP.
#' @param plink_prefix,cm_snp,chr_snp PLINK prefix; position and chromosome of every SNP.
#' @param grm Optional output of the GRM builder; see [.hgsi_units_snp_blocks()].
#' @param maf_min Minimum minor-allele frequency for the SNP-level methods.
#' @param method,rmax,r2,ld_max_cM,max_span_cM,flank_cM,max_snps,min_gamete
#'   Block settings (see [hbr_haplotype_blocks()]).
#' @return `NULL` when no target haplotype has a map position, else a list:
#'   `units`, `block`, `carry` (individuals x target haplotypes), `full` and
#'   `partial` (individuals x blocks), `blocks` (table) and `snp_clusters`
#'   (SNP-level windows used, or `NULL`).
#' @noRd
.hgsi_target_blocks <- function(sup_tabs, codes_l, ids, gene_chr, gene_cm,
                                gene_first, plink_prefix, cm_snp, chr_snp,
                                method, rmax, r2, ld_max_cM, max_span_cM,
                                flank_cM, max_snps, min_gamete,
                                grm = NULL, maf_min = 0.05) {
  traits <- names(sup_tabs)
  tg <- unique(do.call(rbind, lapply(sup_tabs, function(s)
    s[!is.na(s$hap), , drop = FALSE])))
  tg <- tg[tg$gene %in% names(gene_chr), , drop = FALSE]
  if (!nrow(tg)) return(NULL)
  rownames(tg) <- NULL
  tg$traits <- vapply(tg$gene, function(g)
    paste(traits[vapply(traits, function(tr) {
      s <- sup_tabs[[tr]]; g %in% s$gene[!is.na(s$hap)]
    }, logical(1))], collapse = ";"), character(1))
  units <- data.frame(gene = tg$gene, hap = tg$hap, traits = tg$traits,
                      chr = unname(gene_chr[tg$gene]),
                      cm = unname(gene_cm[tg$gene]),
                      snp_col = unname(vapply(tg$gene, function(g)
                        as.integer(gene_first[[g]]), integer(1))),
                      stringsAsFactors = FALSE)
  carry <- .hgsi_carry_matrix(codes_l, ids, units)
  snp_method <- sub("^recomb\\+", "", method)
  clusters <- NULL
  block <- if (method %in% c("fourgamete", "gabriel")) {
    b <- .hgsi_units_snp_blocks(
      plink_prefix = plink_prefix,
      grm          = grm,
      units        = units,
      cm_snp       = cm_snp,
      chr_snp      = chr_snp,
      method       = method,
      flank_cM     = flank_cM,
      max_snps     = max_snps,
      max_span_cM  = max_span_cM,
      min_gamete   = min_gamete,
      maf_min      = maf_min)
    clusters <- attr(b, "clusters"); as.integer(b)
  } else if (method %in% c("recomb+gabriel", "recomb+fourgamete")) {
    b1 <- .hgsi_define_blocks(units$chr, units$cm, carry, "recomb", rmax, r2,
                              ld_max_cM, max_span_cM)
    b2 <- .hgsi_units_snp_blocks(
      plink_prefix = plink_prefix,
      grm          = grm,
      units        = units,
      cm_snp       = cm_snp,
      chr_snp      = chr_snp,
      method       = snp_method,
      flank_cM     = flank_cM,
      max_snps     = max_snps,
      max_span_cM  = max_span_cM,
      min_gamete   = min_gamete,
      maf_min      = maf_min)
    clusters <- attr(b2, "clusters")
    .hgsi_cap_span(.hgsi_merge_blocks(b1, as.integer(b2)), units$chr, units$cm,
                   max_span_cM)
  } else {
    .hgsi_define_blocks(units$chr, units$cm, carry, method, rmax, r2,
                        ld_max_cM, max_span_cM)
  }
  bc <- .hgsi_block_carry(carry, block)
  list(units = units, block = block, carry = carry, full = bc$full,
       partial = bc$partial, blocks = .hgsi_block_table(units, block, carry),
       snp_clusters = clusters)
}

#' Haplotype blocks of target loci
#'
#' Groups a set of target loci into haplotype blocks by one of seven
#' methods and reports the block assignment. Two calling conventions
#' are supported.
#'
#' **Low-level.** Supply `plink_prefix`, `units`, `cm_snp` and
#' `chr_snp` directly. `units` must carry `chr`, `cm` and (for the
#' SNP-level methods) `snp_col`; `cm_snp` and `chr_snp` give the
#' genetic and physical position of every SNP in the PLINK file. The
#' function returns an integer vector of block IDs per row of `units`.
#' This is the legacy API and is fully backward compatible.
#'
#' **Run-level.** Pass an `hbr_run` (from [hbr_read_results()]) as the
#' first positional argument or as `plink_prefix`. The function reads
#' the superior haplotypes, rebuilds the codes, replays the discovery
#' gate to pick the best haplotype per gene, places each at its
#' genetic map position, and returns a classed `hbr_blocks` object
#' with the block table, the per-line carrier matrix, and (for
#' SNP-level methods) the SNP windows used.
#'
#' @details
#' The seven block methods are:
#' \describe{
#'   \item{`"recomb"`}{Join neighbours whose Haldane recombination
#'     fraction is at most `rmax`. Describes the recombination
#'     expected in a new cross.}
#'   \item{`"ld"`}{Join neighbours whose carrier indicators have
#'     r-squared at least `r2` within `ld_max_cM`. Describes LD in the
#'     panel.}
#'   \item{`"either"`}{Union of `"recomb"` and `"ld"`.}
#'   \item{`"fourgamete"`}{Four-gamete rule on the SNPs around the
#'     genes (Wang et al. 2002; Haploview default).}
#'   \item{`"gabriel"`}{Gabriel-type D-prime confidence intervals
#'     (Gabriel et al. 2002).}
#'   \item{`"recomb+gabriel"`}{Union of `"recomb"` and `"gabriel"`,
#'     then capped at `max_span_cM`.}
#'   \item{`"recomb+fourgamete"`}{Union of `"recomb"` and
#'     `"fourgamete"`, then capped.}
#' }
#' Blocks never span chromosomes. `"recomb"`, `"ld"` and `"either"`
#' operate on the loci directly; the SNP-level methods read a window
#' of `flank_cM` around each cluster of target genes.
#'
#' @param plink_prefix PLINK prefix, or an `hbr_run` object (in which
#'   case the run-level API is used). When the first positional
#'   argument is an `hbr_run`, this argument receives it.
#' @param units Low-level API: data.frame with `chr` and `cm` per
#'   target locus, and `snp_col` when using the SNP-level methods.
#'   Alternatively, an `hbr_run` object (in which case the run-level
#'   API is used).
#' @param cm_snp,chr_snp Low-level API: genetic positions and
#'   chromosomes of every SNP in the PLINK file.
#' @param method One of the seven methods above.
#' @param flank_cM,max_snps,min_gamete SNP-level methods only: window
#'   around the genes, SNPs kept per cluster, and gamete frequency that
#'   counts as recombination.
#' @param maf_min Minimum minor-allele frequency.
#' @param rmax Maximum Haldane recombination fraction between
#'   neighbours for `"recomb"`.
#' @param r2 Minimum carrier-indicator r-squared for `"ld"`.
#' @param ld_max_cM Largest distance over which LD may join neighbours.
#' @param max_span_cM Largest distance between the first and last locus
#'   of a block.
#' @param r2_min,gabriel_grid,carry Reserved for future versions; kept
#'   for backward compatibility with the low-level API.
#' @param genetic_map,cM_per_Mb Genetic positions for the run-level API:
#'   a data.frame with `chr`, `bp`, `cm` (for example
#'   [hbr_rice_genetic_map()]), or a constant rate. A `.bim` with an
#'   informative cM column is used when neither is supplied.
#' @param verbose Print progress.
#'
#' @return For the low-level API, an integer vector of block IDs per
#'   row of `units`. For the run-level API, an object of class
#'   `hbr_blocks` with elements `blocks`, `units`, `lines`, `full`,
#'   `snp_clusters`, `settings` and `map_source`.
#'
#' @seealso [hbr_read_results()] to build the run input,
#'   [hbr_haplotype_block_table()] for the legacy low-level wrapper.
#'
#' @examples
#' # Legacy low-level API — no genotypes needed for the map methods
#' units <- data.frame(chr = c("1", "1", "1", "2", "2"),
#'                     cm  = c(10, 10.5, 50, 5, 5.2))
#' hbr_haplotype_blocks(method = "recomb", units = units)
#'
#' # Run-level API on the bundled toy data
#' toy <- system.file("extdata", "hgsi_toy", package = "HaploBreedR")
#' if (nzchar(toy) && dir.exists(toy) &&
#'     requireNamespace("BEDMatrix", quietly = TRUE)) {
#'   run <- hbr_read_results(file.path(toy, "trait_YLD"))
#'   blk <- hbr_haplotype_blocks(
#'     run, plink_prefix = file.path(toy, "geno_toy"),
#'     method = "recomb")
#'   blk$blocks
#' }
#'
#' @export
hbr_haplotype_blocks <- function(
    plink_prefix  = NULL,
    units         = NULL,
    cm_snp        = NULL,
    chr_snp       = NULL,
    method        = c("recomb", "ld", "either",
                      "fourgamete", "gabriel",
                      "recomb+gabriel", "recomb+fourgamete"),
    flank_cM      = 1,
    max_snps      = 400L,
    max_span_cM   = 5,
    min_gamete    = 0.01,
    maf_min       = 0.05,
    rmax          = 0.01,
    r2            = 0.5,
    ld_max_cM     = 3,
    r2_min        = 0.25,
    gabriel_grid  = 21L,
    carry         = NULL,
    genetic_map   = NULL,
    cM_per_Mb     = NULL,
    verbose       = FALSE) {
  
  method <- match.arg(method)
  
  # ---- Detect run-level call ---------------------------------------------
  # R's argument matching puts a positional run object in the *first*
  # unused formal, which may be `units` when `plink_prefix` is supplied
  # by name. Check both positions.
  run_input <- NULL
  if (inherits(plink_prefix, "hbr_run")) {
    run_input <- plink_prefix
    plink_prefix <- NULL
  } else if (inherits(units, "hbr_run")) {
    run_input <- units
    units <- NULL
  }
  
  if (!is.null(run_input)) {
    return(.hbr_haplotype_blocks_run(
      run_input,
      plink_prefix = plink_prefix,
      genetic_map  = genetic_map,
      cM_per_Mb    = cM_per_Mb,
      method       = method,
      rmax         = rmax,
      r2           = r2,
      ld_max_cM    = ld_max_cM,
      max_span_cM  = max_span_cM,
      flank_cM     = flank_cM,
      max_snps     = max_snps,
      min_gamete   = min_gamete,
      maf_min      = maf_min,
      verbose      = verbose))
  }
  
  # ---- Legacy path -------------------------------------------------------
  .hbr_haplotype_blocks_legacy(
    plink_prefix, units, cm_snp, chr_snp, method,
    flank_cM, max_snps, max_span_cM, min_gamete, maf_min,
    rmax, r2, ld_max_cM, r2_min, gabriel_grid, carry)
}


#' Run-level haplotype-block calling
#'
#' Internal worker for the run-level branch of
#' [hbr_haplotype_blocks()]. Reads a pipeline run (or list of runs),
#' rebuilds the superior-haplotype codes, replays the discovery gate,
#' and groups the target loci into blocks. Called only from
#' [hbr_haplotype_blocks()].
#'
#' @param x An `hbr_run`, a named list of them, or anything
#'   [hbr_read_results()] takes.
#' @param plink_prefix PLINK prefix; defaults to the run's stored prefix.
#' @param genetic_map,cM_per_Mb Genetic positions.
#' @param method Block method.
#' @param rmax,r2,ld_max_cM,max_span_cM,flank_cM,max_snps,min_gamete
#'   Block settings.
#' @param maf_min Minimum minor-allele frequency.
#' @param verbose Print progress.
#' @return An object of class `hbr_blocks`.
#' @noRd
.hbr_haplotype_blocks_run <- function(
    x, plink_prefix, genetic_map, cM_per_Mb, method,
    rmax, r2, ld_max_cM, max_span_cM, flank_cM, max_snps,
    min_gamete, maf_min, verbose) {
  
  runs <- .hbr_as_runs(x)
  if (is.null(plink_prefix) || !nzchar(plink_prefix))
    plink_prefix <- runs[[1]]$plink_prefix
  if (is.null(plink_prefix) || !nzchar(plink_prefix) ||
      !file.exists(paste0(plink_prefix, ".fam")))
    stop("`plink_prefix` is not set on the run and not supplied as an ",
         "argument (or does not point to an existing .fam file): ",
         if (is.null(plink_prefix)) "<NULL>" else plink_prefix,
         call. = FALSE)
  
  fam <- data.table::fread(paste0(plink_prefix, ".fam"),
                           header = FALSE, data.table = FALSE)
  ids <- hbr_clean_id(fam$V2)
  n   <- length(ids)
  bim <- .hgsi_read_bim(plink_prefix)
  gpos <- .hgsi_genetic_pos(bim, cM_per_Mb, genetic_map)
  if (is.null(gpos))
    stop("blocks need genetic positions: give `genetic_map`, ",
         "`cM_per_Mb`, or a .bim with a cM column", call. = FALSE)
  
  sup_tabs   <- list()
  codes_l    <- list()
  gene_first <- list()
  
  for (tr in names(runs)) {
    ld  <- runs[[tr]]
    sup <- intersect(ld$superior, names(ld$gene_snps))
    if (!length(sup)) {
      sup_tabs[[tr]] <- data.frame(gene = character(),
                                   hap = character())
      next
    }
    codes <- .hgsi_build_codes(plink_prefix, ld$gene_snps[sup], ids,
                               bim$snp, maf_min = ld$gate$maf_min,
                               verbose = verbose)
    codes <- .hgsi_relabel_codes(codes, ids, ld)
    y <- rep(NA_real_, n)
    j <- match(ids, ld$phenotype$ID)
    y[!is.na(j)] <- ld$phenotype$y[j[!is.na(j)]]
    pcs <- if (is.null(ld$pcs)) NULL else
      as.matrix(ld$pcs[match(ids, ld$pcs$ID), ld$pc_cols, drop = FALSE])
    gt <- .hgsi_gate_fold(codes, y, which(!is.na(y)), ld$gate, pcs)
    sup_tabs[[tr]] <- data.frame(
      gene = sup,
      hap  = gt$best_hap[match(sup, gt$gene)],
      stringsAsFactors = FALSE)
    codes_l[[tr]] <- codes
    for (g in sup) if (is.null(gene_first[[g]])) {
      ii <- stats::na.omit(match(ld$gene_snps[[g]], bim$snp))
      if (length(ii)) gene_first[[g]] <- ii[1]
    }
  }
  
  if (!length(gene_first))
    stop("no superior gene with SNPs in the PLINK file", call. = FALSE)
  
  gene_chr <- vapply(gene_first, function(i) bim$chr[i], character(1))
  gene_cm  <- vapply(gene_first, function(i) gpos$cm[i], numeric(1))
  
  tgt <- .hgsi_target_blocks(
    sup_tabs, codes_l, ids, gene_chr, gene_cm, gene_first,
    plink_prefix, gpos$cm, bim$chr, method, rmax, r2,
    ld_max_cM, max_span_cM, flank_cM, max_snps, min_gamete,
    grm = NULL, maf_min = maf_min)
  if (is.null(tgt))
    stop("no superior haplotype with a map position", call. = FALSE)
  
  full <- tgt$full
  colnames(full) <- tgt$blocks$block
  rownames(full) <- ids
  
  units_out <- tgt$units
  units_out$block   <- sprintf("HB%03d", tgt$block)
  units_out$snp_col <- NULL
  
  structure(list(
    blocks = tgt$blocks,
    units  = units_out,
    lines  = data.frame(
      id = ids,
      n_blocks_full    = rowSums(tgt$full),
      n_blocks_partial = rowSums(tgt$partial),
      n_blocks_target  = ncol(tgt$full),
      stringsAsFactors = FALSE),
    full = full,
    snp_clusters = tgt$snp_clusters,
    settings = list(method = method, rmax = rmax, r2 = r2,
                    ld_max_cM = ld_max_cM, max_span_cM = max_span_cM,
                    flank_cM = flank_cM, max_snps = max_snps,
                    min_gamete = min_gamete),
    map_source = gpos$source),
    class = "hbr_blocks")
}


#' Legacy low-level haplotype-block calling
#'
#' The original body of `hbr_haplotype_blocks()` before it became
#' polymorphic. Called by [hbr_haplotype_blocks()] and by
#' [hbr_haplotype_block_table()].
#'
#' @param plink_prefix PLINK prefix.
#' @param units Data.frame with `chr`, `cm` and optionally `snp_col`.
#' @param cm_snp,chr_snp Genetic position and chromosome per SNP.
#' @param method Block method.
#' @param flank_cM,max_snps,max_span_cM,min_gamete,maf_min Block settings.
#' @param rmax,r2,ld_max_cM Block thresholds.
#' @param r2_min,gabriel_grid,carry Reserved.
#' @return Integer block IDs per row of `units`.
#' @noRd
.hbr_haplotype_blocks_legacy <- function(
    plink_prefix, units, cm_snp, chr_snp, method,
    flank_cM, max_snps, max_span_cM, min_gamete, maf_min,
    rmax, r2, ld_max_cM, r2_min, gabriel_grid, carry) {
  
  if (!is.data.frame(units) && !data.table::is.data.table(units))
    stop("`units` must be a data.frame.", call. = FALSE)
  need <- c("chr", "cm")
  miss <- setdiff(need, names(units))
  if (length(miss))
    stop("`units` must have columns: ", paste(miss, collapse = ", "),
         call. = FALSE)
  units <- as.data.frame(units)
  if (nrow(units) == 0L) return(integer(0))
  
  uses_snp <- method %in% c("fourgamete", "gabriel",
                            "recomb+gabriel", "recomb+fourgamete")
  uses_ld  <- method %in% c("ld", "either")
  
  if (uses_snp) {
    if (is.null(plink_prefix) || !nzchar(plink_prefix))
      stop("Method '", method, "' needs `plink_prefix`.", call. = FALSE)
    if (is.null(cm_snp) || is.null(chr_snp))
      stop("Method '", method, "' needs `cm_snp` and `chr_snp`.",
           call. = FALSE)
    if (!"snp_col" %in% names(units))
      stop("Method '", method, "' needs `units$snp_col`.",
           call. = FALSE)
  }
  if (uses_ld && is.null(carry))
    stop("Method '", method, "' needs a `carry` matrix.",
         call. = FALSE)
  if (uses_ld && ncol(carry) != nrow(units))
    stop("`carry` has ", ncol(carry), " columns but `units` has ",
         nrow(units), " rows.", call. = FALSE)
  
  if (method %in% c("recomb", "ld", "either")) {
    return(.hgsi_define_blocks(
      chr = units$chr, cm = units$cm,
      carry = if (identical(method, "recomb")) NULL else carry,
      method = method,
      rmax = rmax, r2_min = r2, ld_max_cM = ld_max_cM,
      max_span_cM = max_span_cM))
  }
  
  if (method %in% c("fourgamete", "gabriel")) {
    return(.hgsi_units_snp_blocks(
      plink_prefix = plink_prefix, units = units,
      cm_snp = cm_snp, chr_snp = chr_snp, method = method,
      flank_cM = flank_cM, max_snps = max_snps,
      max_span_cM = max_span_cM, min_gamete = min_gamete,
      maf_min = maf_min))
  }
  
  snp_method <- sub("^recomb\\+", "", method)
  b_map <- .hgsi_define_blocks(
    chr = units$chr, cm = units$cm, carry = NULL, method = "recomb",
    rmax = rmax, r2_min = r2, ld_max_cM = ld_max_cM,
    max_span_cM = max_span_cM)
  b_snp <- .hgsi_units_snp_blocks(
    plink_prefix = plink_prefix, units = units,
    cm_snp = cm_snp, chr_snp = chr_snp, method = snp_method,
    flank_cM = flank_cM, max_snps = max_snps,
    max_span_cM = max_span_cM, min_gamete = min_gamete,
    maf_min = maf_min)
  merged <- .hgsi_merge_blocks(b_map, as.integer(b_snp))
  capped <- .hgsi_cap_span(merged, units$chr, units$cm, max_span_cM)
  attr(capped, "clusters") <- attr(b_snp, "clusters")
  capped
}


#' Print an `hbr_blocks` result
#'
#' @param x An object of class `hbr_blocks`.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @examples
#' toy <- system.file("extdata", "hgsi_toy", package = "HaploBreedR")
#' if (nzchar(toy) && dir.exists(toy) &&
#'     requireNamespace("BEDMatrix", quietly = TRUE)) {
#'   run <- hbr_read_results(file.path(toy, "trait_YLD"))
#'   print(hbr_haplotype_blocks(
#'     run, plink_prefix = file.path(toy, "geno_toy"),
#'     method = "recomb"))
#' }
#' @export
print.hbr_blocks <- function(x, ...) {
  cat(sprintf("%d target haplotypes in %d blocks (method %s; map: %s)\n",
              nrow(x$units), nrow(x$blocks), x$settings$method,
              x$map_source))
  print(utils::head(x$blocks[, intersect(
    c("block", "chr", "start_cM", "end_cM", "n_genes", "panel_freq"),
    names(x$blocks))], 15), row.names = FALSE)
  if (nrow(x$blocks) > 15) cat("  ... see $blocks\n")
  sizes <- table(factor(pmin(x$blocks$n_target_haplotypes, 5L),
                        levels = 1:5,
                        labels = c("1", "2", "3", "4", "5+")))
  cat("  target haplotypes per block:",
      paste(names(sizes), sizes, sep = ": ", collapse = " | "), "\n")
  cl <- x$snp_clusters
  if (!is.null(cl))
    cat(sprintf("  SNP-level windows: %d\n", nrow(cl)))
  invisible(x)
}


#' Haplotype block table from a set of target loci
#'
#' Runs [hbr_haplotype_blocks()] and returns the block summary table
#' alongside the per-locus assignment. The table has one row per block with
#' the block id, chromosome, start and end position in cM, the number of
#' target loci in the block, the number of distinct genes (when `units` has
#' a `gene` column) and the concatenated gene and trait names.
#'
#' @inheritParams hbr_haplotype_blocks
#' @param traits Optional character vector (one per row of `units`) of the
#'   trait(s) each locus belongs to. When `units` already carries a `traits`
#'   column, that is used instead.
#'
#' @return A list with `blocks` (the summary table) and `assignment`
#'   (the block id per target locus, in the same order as `units`).
#'
#' @examples
#' units <- data.frame(chr = c("1", "1", "1"),
#'                     cm  = c(10, 10.5, 11),
#'                     gene = c("g1", "g2", "g3"))
#' hbr_haplotype_block_table(method = "recomb", units = units)
#'
#' @export
hbr_haplotype_block_table <- function(
    plink_prefix  = NULL,
    units,
    cm_snp        = NULL,
    chr_snp       = NULL,
    traits        = NULL,
    method        = c("recomb", "ld", "either",
                      "fourgamete", "gabriel",
                      "recomb+gabriel", "recomb+fourgamete"),
    flank_cM      = 1,
    max_snps      = 400L,
    max_span_cM   = 5,
    min_gamete    = 0.01,
    maf_min       = 0.05,
    rmax          = 0.01,
    r2            = 0.5,
    ld_max_cM     = 3,
    r2_min        = 0.25,
    gabriel_grid  = 21L,
    carry         = NULL) {

  method <- match.arg(method)
  block <- hbr_haplotype_blocks(
    plink_prefix = plink_prefix, units = units,
    cm_snp = cm_snp, chr_snp = chr_snp, method = method,
    flank_cM = flank_cM, max_snps = max_snps,
    max_span_cM = max_span_cM, min_gamete = min_gamete,
    maf_min = maf_min, rmax = rmax, r2 = r2,
    ld_max_cM = ld_max_cM, r2_min = r2_min,
    gabriel_grid = gabriel_grid, carry = carry)

  u <- as.data.frame(units)
  if (!"gene" %in% names(u)) u$gene <- sprintf("locus_%d", seq_len(nrow(u)))
  if (!"traits" %in% names(u)) {
    u$traits <- if (!is.null(traits)) as.character(traits) else NA_character_
  }

  blocks <- .hgsi_block_table(u, block, carry_panel = NULL)
  attr(blocks, "method") <- method
  list(blocks = blocks, assignment = as.integer(block))
}
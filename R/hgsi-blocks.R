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

#' Standalone Haplotype Block Calling
#'
#' Partitions a set of target loci into haplotype blocks by any of the
#' methods the pipeline uses: map-based recombination between the loci
#' (`"recomb"`), LD between their carrier indicators (`"ld"`), the union of
#' both (`"either"`), SNP-level methods around each target cluster
#' (`"fourgamete"`, `"gabriel"`), or the union of a map method and an
#' SNP method (`"recomb+gabriel"`, `"recomb+fourgamete"`). Blocks never span
#' chromosomes and never exceed `max_span_cM`.
#'
#' @param plink_prefix Character. PLINK file prefix. Required by the
#'   SNP-level methods (`"fourgamete"`, `"gabriel"`) and their combinations;
#'   ignored otherwise.
#' @param units Data.frame with columns `chr` and `cm`, and optionally
#'   `snp_col` (PLINK column index of the locus's own SNP, needed by
#'   `"fourgamete"` and `"gabriel"`).
#' @param cm_snp Numeric vector of genetic positions per SNP in the PLINK
#'   file. Required by the SNP-level methods.
#' @param chr_snp Character vector of chromosomes per SNP in the PLINK
#'   file. Required by the SNP-level methods.
#' @param method Block method, one of `"recomb"`, `"ld"`, `"either"`,
#'   `"fourgamete"`, `"gabriel"`, `"recomb+gabriel"`, `"recomb+fourgamete"`.
#' @param flank_cM Numeric. Window (+/- cM) around each cluster of target
#'   loci from which SNPs are read, for SNP-level methods.
#' @param max_snps Integer. Maximum SNPs per cluster after thinning.
#' @param max_span_cM Numeric. Largest distance between the first and last
#'   locus of a block (cM). Chaining cannot build blocks longer than this.
#' @param min_gamete Numeric. Four-gamete frequency that counts as evidence
#'   of recombination (Haploview default 0.01).
#' @param maf_min Numeric. Minimum minor allele frequency for a SNP to
#'   enter the block search.
#' @param rmax Numeric. Maximum recombination fraction between
#'   neighbouring loci for `"recomb"` and `"either"` (default 0.01, about
#'   1 cM).
#' @param r2 Numeric. Minimum r^2 of carrier indicators to join neighbours
#'   under `"ld"` and `"either"`.
#' @param ld_max_cM Numeric. Largest distance over which LD may join
#'   neighbours (cM).
#' @param r2_min Numeric. Reserved for a future version; currently not
#'   wired into the SNP-level block callers (the four-gamete rule and the
#'   full Gabriel D' test are applied directly to every pair).
#' @param gabriel_grid Integer. Reserved for a future version; currently
#'   not wired into the SNP-level block callers.
#' @param carry Logical matrix individuals x target loci. Required by
#'   `"ld"` and `"either"`; ignored otherwise.
#'
#' @return Integer vector of block IDs per target locus, with attribute
#'   `"clusters"` giving the SNP windows used by SNP-level methods (or
#'   `NULL`).
#'
#' @examples
#' units <- data.frame(chr = c("1", "1", "1", "2", "2"),
#'                     cm  = c(10, 10.5, 50, 5, 5.2))
#' hbr_haplotype_blocks(method = "recomb", units = units)
#'
#' @export
hbr_haplotype_blocks <- function(
    plink_prefix  = NULL,
    units,
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
    carry         = NULL) {

  method <- match.arg(method)

  # --- Validate units -----------------------------------------------------
  if (!is.data.frame(units) && !data.table::is.data.table(units))
    stop("`units` must be a data.frame.", call. = FALSE)
  need <- c("chr", "cm")
  miss <- setdiff(need, names(units))
  if (length(miss))
    stop("`units` must have columns: ", paste(miss, collapse = ", "),
         call. = FALSE)
  units <- as.data.frame(units)
  if (nrow(units) == 0L) return(integer(0))

  # --- Determine what the chosen method needs -----------------------------
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
      stop("Method '", method, "' needs `units$snp_col` ",
           "(PLINK column of each locus's own SNP).", call. = FALSE)
  }
  if (uses_ld && is.null(carry))
    stop("Method '", method, "' needs a `carry` matrix of individuals ",
         "by target loci (ncol = nrow(units)); pass `carry =` or use a ",
         "method that does not require LD.", call. = FALSE)
  if (uses_ld && ncol(carry) != nrow(units))
    stop("`carry` has ", ncol(carry), " columns but `units` has ",
         nrow(units), " rows.", call. = FALSE)

  # --- Dispatch ------------------------------------------------------------
  if (method %in% c("recomb", "ld", "either")) {
    return(.hgsi_define_blocks(
      chr = units$chr, cm = units$cm,
      carry = if (identical(method, "recomb")) NULL else carry,
      method = method,
      rmax = rmax, r2_min = r2, ld_max_cM = ld_max_cM,
      max_span_cM = max_span_cM))
  }

  if (method %in% c("fourgamete", "gabriel")) {
    out <- .hgsi_units_snp_blocks(
      plink_prefix = plink_prefix,
      units        = units,
      cm_snp       = cm_snp,
      chr_snp      = chr_snp,
      method       = method,
      flank_cM     = flank_cM,
      max_snps     = max_snps,
      max_span_cM  = max_span_cM,
      min_gamete   = min_gamete,
      maf_min      = maf_min)
    return(out)
  }

  # --- Combinations: recomb + (gabriel | fourgamete) ----------------------
  snp_method <- sub("^recomb\\+", "", method)
  b_map <- .hgsi_define_blocks(
    chr = units$chr, cm = units$cm, carry = NULL,
    method = "recomb",
    rmax = rmax, r2_min = r2, ld_max_cM = ld_max_cM,
    max_span_cM = max_span_cM)
  b_snp <- .hgsi_units_snp_blocks(
    plink_prefix = plink_prefix,
    units        = units,
    cm_snp       = cm_snp,
    chr_snp      = chr_snp,
    method       = snp_method,
    flank_cM     = flank_cM,
    max_snps     = max_snps,
    max_span_cM  = max_span_cM,
    min_gamete   = min_gamete,
    maf_min      = maf_min)
  merged <- .hgsi_merge_blocks(b_map, as.integer(b_snp))
  capped <- .hgsi_cap_span(merged, units$chr, units$cm, max_span_cM)
  attr(capped, "clusters") <- attr(b_snp, "clusters")
  capped
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
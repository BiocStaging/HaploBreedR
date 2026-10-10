# ============================================================================
# data-raw/make_hgsi_toy.R
# Generates the toy data used in the hbr_advanced_decision() example.
# Run from the package root:  source("data-raw/make_hgsi_toy.R")
#
# Writes inst/extdata/hgsi_toy/:
#   geno_toy.{bed,bim,fam}  150 inbred lines, 600 SNPs, 3 chromosomes x 100 cM
#   trait_YLD/              the files haplobreed_analysis() writes for one trait
# Simulation: mosaic lines of 6 founders; 8 genes with 3 adjacent functional
# SNPs; 3 genes carry haplotype effects; polygenic background; trait YLD,
# direction "higher", no PCs, 10 lines left unphenotyped. The superior-haplotype
# gate uses the pipeline's rules: one-way F-test, Hedges' g of best vs runner-up,
# BH-FDR; alpha 0.05, minimum g 0.3, at least 3 individuals per haplotype.
# ============================================================================

set.seed(2026)
root <- file.path("inst", "extdata", "hgsi_toy")
tdir <- file.path(root, "trait_YLD")
dir.create(file.path(tdir, "superior_genes"), recursive = TRUE,
           showWarnings = FALSE)

n <- 150L; n_chr <- 3L; per_chr <- 200L; n_found <- 6L
ids    <- sprintf("L%03d", seq_len(n))
cm_chr <- seq(0.25, 99.75, length.out = per_chr)          # 0.5 cM spacing
m      <- n_chr * per_chr
chr_of <- rep(seq_len(n_chr), each = per_chr)
cm     <- rep(cm_chr, n_chr)
bp     <- as.integer(cm * 1e4) + 1000L
snp    <- paste0("snp", seq_len(m))

# ---- genotypes: mosaics of founders (inbred, dosage 0/2) --------------------
founder <- array(sample(0:1, n_chr * n_found * per_chr, replace = TRUE),
                 dim = c(n_chr, n_found, per_chr))
D <- matrix(0L, n, m)
for (i in seq_len(n)) {
  for (ch in seq_len(n_chr)) {
    nb   <- rpois(1, 3)
    br   <- sort(runif(nb, 0, 100))
    seg  <- sample.int(n_found, nb + 1L, replace = TRUE)
    f_at <- seg[findInterval(cm_chr, br) + 1L]
    cols <- (ch - 1L) * per_chr + seq_len(per_chr)
    D[i, cols] <- 2L * founder[cbind(ch, f_at, seq_len(per_chr))]
  }
}

# ---- genes: 3 adjacent functional SNPs each ---------------------------------
gene_pos <- data.frame(chr = c(1, 1, 1, 2, 2, 2, 3, 3),
                       cm  = c(20, 55, 85, 15, 50, 80, 30, 70))
genes <- list()
for (k in seq_len(nrow(gene_pos))) {
  j <- (gene_pos$chr[k] - 1L) * per_chr + which.min(abs(cm_chr - gene_pos$cm[k]))
  genes[[sprintf("GENE%02d", k)]] <- (j - 1L):(j + 1L)
}
hapstr <- lapply(genes, function(cols)
  apply(D[, cols, drop = FALSE], 1, paste, collapse = "-"))

# ---- trait ------------------------------------------------------------------
causal <- c(GENE01 = 3, GENE04 = 2.5, GENE07 = 2)
y <- rep(50, n)
for (g in names(causal)) {
  lev <- sort(unique(hapstr[[g]]))
  eff <- stats::setNames(rnorm(length(lev), 0, causal[[g]]), lev)
  y <- y + unname(eff[hapstr[[g]]])
}
poly <- sample.int(m, 60)
y <- y + as.numeric((D[, poly] - 1) %*% rnorm(60, 0, 0.25)) + rnorm(n, 0, 1.5)
y <- round(y, 3)
pheno <- rep(TRUE, n); pheno[sample.int(n, 10)] <- FALSE

# ---- PLINK files ------------------------------------------------------------
con <- file(file.path(root, "geno_toy.bed"), "wb")
writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
bps <- (n + 3L) %/% 4L
for (j in seq_len(m)) {
  code  <- ifelse(D[, j] == 2L, 0L, 3L)                   # 00 = A1/A1, 11 = A2/A2
  code  <- c(code, integer(bps * 4L - n))
  mat   <- matrix(code, nrow = 4L)                        # one column per byte
  bytes <- mat[1, ] + bitwShiftL(mat[2, ], 2L) +
    bitwShiftL(mat[3, ], 4L) + bitwShiftL(mat[4, ], 6L)
  writeBin(as.raw(bytes), con)
}
close(con)

# A1 = ALT, A2 = REF. The functional-SNP table below declares REF = "A", ALT = "G",
# so the .bim must store A1 = "G" and A2 = "A". Swapping these two characters was
# the source of an inverted-allele bug (dosage 2 was being read as REF).
data.table::fwrite(
  data.frame(chr_of, snp, round(cm, 3), bp, "G", "A"),
  file.path(root, "geno_toy.bim"), sep = "\t", col.names = FALSE, quote = FALSE)
data.table::fwrite(
  data.frame(ids, ids, 0L, 0L, 0L, -9L),
  file.path(root, "geno_toy.fam"), sep = " ", col.names = FALSE, quote = FALSE)

# ---- pipeline-style gate ----------------------------------------------------
alpha <- 0.05; min_g <- 0.30; min_ind <- 3L
ph_ids <- ids[pheno]; yp <- y[pheno]
rows <- list(); groups_of <- list()
for (g in names(genes)) {
  h    <- hapstr[[g]][pheno]
  cnt  <- table(h)
  keep <- names(cnt)[cnt >= min_ind]
  if (length(keep) < 2L) next
  sel <- h %in% keep
  hh <- h[sel]; yy <- yp[sel]; idp <- ph_ids[sel]
  tb  <- table(hh)
  lab <- stats::setNames(paste0("Hap", seq_along(tb)),
                         names(tb)[order(-as.integer(tb), names(tb))])
  gm <- tapply(yy, hh, mean); gn <- tapply(yy, hh, length)
  o  <- order(-gm, -gn)
  best <- names(gm)[o[1]]; second <- names(gm)[o[2]]
  p  <- stats::anova(stats::lm(yy ~ factor(hh)))[["Pr(>F)"]][1]
  a <- yy[hh == best]; b <- yy[hh == second]
  sp <- sqrt(((length(a) - 1) * stats::var(a) + (length(b) - 1) * stats::var(b)) /
               (length(a) + length(b) - 2))
  gg <- abs(mean(a) - mean(b)) / sp *
    (1 - 3 / (4 * (length(a) + length(b) - 2) - 1))          # Hedges' g
  rows[[g]] <- data.frame(
    Gene_ID = g, Chr = chr_of[genes[[g]][1]], N_Functional_SNPs_Final = 3L,
    N_Haplotypes = length(gm), Test = "ANOVA+duncan", P_value = p,
    Superior_Haplotype = unname(lab[best]), Superior_NT_Seq = best,
    Superior_Mean = round(gm[[best]], 3), Cohens_d = round(gg, 3),
    Second_Haplotype = unname(lab[second]), stringsAsFactors = FALSE)
  groups_of[[g]] <- data.frame(ID = idp, NT_Seq = hh, AA_Seq = hh,
                               Haplotype = unname(lab[hh]), YLD = yy,
                               stringsAsFactors = FALSE)
}
tab <- do.call(rbind, rows); rownames(tab) <- NULL
tab$Adjusted_P  <- signif(stats::p.adjust(tab$P_value, "BH"), 4)
tab$P_value     <- signif(tab$P_value, 4)
tab$Permutation_P <- pmax(tab$P_value, 0.001)   # stand-in for the pipeline's permutation P
tab$Is_Superior <- tab$Adjusted_P < alpha & tab$Cohens_d >= min_g
if (!any(tab$Is_Superior))
  stop("no superior gene in the toy data; change set.seed() or the effect sizes")

data.table::fwrite(tab, file.path(tdir, "publication_table_functional.csv"))
for (g in tab$Gene_ID[tab$Is_Superior])
  data.table::fwrite(groups_of[[g]],
                     file.path(tdir, "superior_genes", paste0(g, "_individuals.csv")))

fs <- do.call(rbind, lapply(names(genes), function(g) data.frame(
  gene = g, SNP = snp[genes[[g]]], Chr = chr_of[genes[[g]]], POS = bp[genes[[g]]],
  REF = "A", ALT = "G", Region = "CDS_nonsynonymous",
  Genomic_Order = seq_along(genes[[g]]), stringsAsFactors = FALSE)))
data.table::fwrite(fs, file.path(tdir, "functional_snps_all_tested_genes.csv"))
data.table::fwrite(data.frame(ID = ph_ids, YLD = yp),
                   file.path(tdir, "phenotype_data_used.csv"))
data.table::fwrite(data.frame(
  Mode = "DirectGene", Annotation_Source = "CUSTOM", Trait = "YLD",
  Direction = "higher", PLINK_Prefix = file.path(root, "geno_toy"),
  Genes_Tested = nrow(tab),
  Superior_Haplotypes = sum(tab$Is_Superior), Schema_Version = "2",
  FDR_Required = TRUE, Alpha_Gate = alpha, Min_Cohens_D_Gate = min_g,
  Min_Ind_Per_Hap_Gate = min_ind, Min_Hap_Freq_Gate = 0, MAF_Min_Gate = 0.05,
  PC_Cols_Gate = "", Trait_Threshold_Gate = NA_real_, Require_FDR_Gate = TRUE,
  Output_Dir = tdir, stringsAsFactors = FALSE),
  file.path(tdir, "pipeline_summary.csv"))

print(tab[, c("Gene_ID", "N_Haplotypes", "P_value", "Adjusted_P",
              "Cohens_d", "Is_Superior")])
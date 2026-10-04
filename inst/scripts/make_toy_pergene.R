# ============================================================================
# inst/scripts/make_toy_pergene.R
# Generates per-gene CSVs for the toy dataset so plotting examples run
# without a full pipeline. Run from the package root:
#   Rscript inst/scripts/make_toy_pergene.R
# ============================================================================

out <- file.path("inst", "extdata", "toy")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
gene <- "gene1"

# Two haplotypes, three individuals each
hap_map <- data.frame(
  Haplotype   = c("Hap1", "Hap2"),
  Count       = c(6L, 4L),
  Frequency   = c(0.6, 0.4),
  NT_Sequence = c("A-C-G", "A-T-G"),
  AA_Sequence = c("Met-Ala-Val", "Met-Ser-Val"))
write.csv(hap_map, file.path(out, paste0(gene, "_hap_map.csv")),
          row.names = FALSE)

groups <- data.frame(
  Haplotype          = c("Hap1", "Hap2"),
  Mean               = c(110.5, 92.3),
  SE                 = c(4.2, 5.1),
  CI_Lower           = c(102.1, 82.1),
  CI_Upper           = c(118.9, 102.5),
  Significance_Group = c("a", "b"))
write.csv(groups, file.path(out, paste0(gene, "_groups.csv")),
          row.names = FALSE)

func_snps <- data.frame(
  SNP       = c("snp005", "snp020", "snp045"),
  POS       = c(1800L, 4800L, 9800L),
  REF       = c("A", "C", "G"),
  ALT       = c("T", "T", "A"),
  Chr       = c("1", "1", "1"),
  Strand    = c("+", "+", "+"),
  Region    = c("CDS_nonsynonymous", "5UTR", "3UTR"),
  AA_Change = c("A76S", NA, NA),
  Codon_Num = c(76L, NA, NA))
write.csv(func_snps, file.path(out, paste0(gene, "_functional_snps.csv")),
          row.names = FALSE)

ind <- data.frame(
  ID          = sprintf("IND%02d", 1:10),
  Haplotype   = c(rep("Hap1", 6), rep("Hap2", 4)),
  trait_value = c(108, 112, 105, 115, 110, 113,
                  90, 95, 88, 96))
write.csv(ind, file.path(out, paste0(gene, "_individuals.csv")),
          row.names = FALSE)

# ---- Long-format gene/haplotype summary consumed by hbr_plot_haplotype() ----
# 3 SNPs per haplotype, matching nt_sequence "A-C-G" / "A-T-G"
viz <- data.frame(
  gene_id     = rep(gene, 2L),
  chromosome  = rep("1", 2L),
  gene_start  = rep(1000L, 2L),
  gene_end    = rep(2500L, 2L),
  haplotype   = c("Hap1", "Hap2"),
  n           = c(6L, 4L),
  mean        = c(110.5, 92.3),
  se          = c(4.2, 5.1),
  ci_lower    = c(102.1, 82.1),
  ci_upper    = c(118.9, 102.5),
  nt_sequence = c("A-C-G", "A-T-G"),
  aa_sequence = c("Met-Ala-Val", "Met-Ser-Val"),
  snp_1_id    = c("snp005", "snp005"),
  snp_2_id    = c("snp020", "snp020"),
  snp_3_id    = c("snp045", "snp045"),
  snp_1_pos   = c(1800L, 1800L),
  snp_2_pos   = c(4800L, 4800L),
  snp_3_pos   = c(9800L, 9800L),
  snp_1_allele= c("A",     "A"),
  snp_2_allele= c("C",     "T"),
  snp_3_allele= c("G",     "G"),
  stringsAsFactors = FALSE)

# --- metadata columns the gene-object parser reads (hbr_build_gene_object) ---
viz$superior_haplotype <- "Hap1"
viz$is_superior        <- TRUE
viz$p_value            <- 1e-8
viz$permutation_p      <- 1e-3
viz$cohens_d           <- 1.5
viz$improvement_pct    <- 12.3
write.csv(viz, file.path(out, paste0(gene, "_viz.csv")), row.names = FALSE)

cat("Per-gene toy files written to", out, "\n")
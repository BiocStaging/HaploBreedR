# ============================================================================
# inst/scripts/make_test_data.R
# Generate a synthetic PLINK trio + pipeline outputs for testing the engine.
# Run once from the package root: Rscript inst/scripts/make_test_data.R
# ============================================================================

set.seed(42)
out_dir <- file.path("inst", "extdata", "test_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_ind  <- 100L
n_snp  <- 500L
ids    <- sprintf("IND%03d", seq_len(n_ind))
snp_ids <- sprintf("snp%04d", seq_len(n_snp))

# ---- Genotypes (0/1/2 with NA) ----
geno <- matrix(sample(c(0L, 1L, 2L, NA_integer_), n_ind * n_snp,
                      replace = TRUE, prob = c(0.4, 0.1, 0.4, 0.1)),
               n_ind, n_snp, dimnames = list(ids, snp_ids))

# ---- PLINK .raw (numeric dosage 0/1/2, NA) ----
raw <- data.frame(FID = ids, IID = ids, PAT = 0, MAT = 0, SEX = 0,
                  PHENOTYPE = -9, check.names = FALSE)
raw <- cbind(raw, as.data.frame(geno, check.names = FALSE))
write.table(raw, file.path(out_dir, "test.raw"),
            sep = " ", row.names = FALSE, quote = FALSE, na = "NA")

# ---- PLINK .bim, .fam ----
# A1 = ALT, A2 = REF; make sure the two alleles differ per SNP.
bases <- c("A", "C", "G", "T")
A2 <- sample(bases, n_snp, replace = TRUE)
A1 <- vapply(A2, function(b) sample(setdiff(bases, b), 1), character(1))
bim <- data.frame(CHR = 1L, SNP = snp_ids, CM = 0,
                  BP = seq(1000L, by = 200L, length.out = n_snp),
                  A1 = A1, A2 = A2)
write.table(bim, file.path(out_dir, "test.bim"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

fam <- data.frame(FID = ids, IID = ids, PAT = 0, MAT = 0, SEX = 0, PHENO = -9)
write.table(fam, file.path(out_dir, "test.fam"),
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

# ---- PLINK .bed (SNP-major) ----
write_bed <- function(geno, path) {
  n <- nrow(geno); m <- ncol(geno)
  bps <- (n + 3L) %/% 4L
  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)
  for (j in seq_len(m)) {
    code <- ifelse(is.na(geno[, j]), 1L,
                   ifelse(geno[, j] == 2L, 0L, 3L))
    code <- c(code, integer(bps * 4L - n))
    mat  <- matrix(code, nrow = 4L)
    bytes <- mat[1, ] + bitwShiftL(mat[2, ], 2L) +
             bitwShiftL(mat[3, ], 4L) + bitwShiftL(mat[4, ], 6L)
    writeBin(as.raw(bytes), con)
  }
}
write_bed(geno, file.path(out_dir, "test.bed"))

# ---- Phenotype tables ----
pheno <- data.frame(ID = ids,
                    trait1 = rnorm(n_ind, 100, 15),
                    trait2 = rnorm(n_ind, 50, 8))
write.csv(pheno, file.path(out_dir, "test_pheno.csv"), row.names = FALSE)

# phenotype_data_used.csv — read by .hgsi_load_trait()
data.table::fwrite(pheno[, c("ID", "trait1")],
                   file.path(out_dir, "phenotype_data_used.csv"))

# ---- Functional SNPs per tested gene ----
gene_snps <- data.table::data.table(
  gene = rep(sprintf("gene%03d", 1:20), each = 5),
  SNP  = sprintf("snp%04d", sample(1:n_snp, 100)),
  Chr  = "1",
  POS  = NA_integer_,
  REF  = "A",
  ALT  = "G",
  Region        = "CDS_nonsynonymous",
  Genomic_Order = rep(1:5, 20))
data.table::fwrite(gene_snps,
                   file.path(out_dir, "functional_snps_all_tested_genes.csv"))

# ---- publication_table_functional.csv — read by .hgsi_read_tested() ----
sup_genes <- sprintf("gene%03d", 1:5)
pub <- data.frame(
  Gene_ID            = sprintf("gene%03d", 1:20),
  Chr                = "1",
  N_Haplotypes       = sample(2:4, 20, replace = TRUE),
  Test               = "ANOVA+duncan",
  P_value            = c(runif(5, 1e-4, 0.01), runif(15, 0.05, 0.9)),
  Adjusted_P         = c(runif(5, 1e-3, 0.04), runif(15, 0.1,  0.9)),
  Cohens_d           = c(runif(5, 0.5,  2.0),  runif(15, 0.1,  0.4)),
  Superior_Haplotype = c(rep("Hap1", 5), rep(NA, 15)),
  Is_Superior        = c(rep(TRUE, 5),   rep(FALSE, 15)),
  stringsAsFactors   = FALSE)
data.table::fwrite(pub,
                   file.path(out_dir, "publication_table_functional.csv"))

# ---- pipeline_summary.csv — read by .hgsi_read_summary() ----
data.table::fwrite(data.frame(
  Mode                 = "GWAS",
  Annotation_Source    = "CUSTOM",
  Trait                = "trait1",
  Direction            = "higher",
  PLINK_Prefix         = file.path(out_dir, "test"),
  Output_Dir           = out_dir,
  Genes_Tested         = 20L,
  Superior_Haplotypes  = 5L,
  Schema_Version       = "2",
  FDR_Required         = TRUE,
  Alpha_Gate           = 0.05,
  Min_Cohens_D_Gate    = 0.3,
  Min_Ind_Per_Hap_Gate = 3L,
  Min_Hap_Freq_Gate    = 0,
  MAF_Min_Gate         = 0.05,
  PC_Cols_Gate         = "",
  Trait_Threshold_Gate = NA_real_,
  Require_FDR_Gate     = TRUE,
  stringsAsFactors     = FALSE),
  file.path(out_dir, "pipeline_summary.csv"))

# ---- In-memory RDS (also usable via hbr_fit_models(pipeline_result=...)) ----
x <- list(
  superior_results = stats::setNames(lapply(1:5, function(i) {
    list(test_res = list(data = data.frame(ID = ids, y = pheno$trait1)),
         hap_data = list(haplotypes = data.frame(
           ID = ids,
           Haplotype = sample(c("H1", "H2"), n_ind, replace = TRUE))),
         superior = list(haplotype = "H1"))
  }), sup_genes),
  functional_snps_all = gene_snps,
  phenotype_data      = pheno[, c("ID", "trait1")],
  pipeline_summary    = data.frame(
    Trait               = "trait1",
    Direction           = "higher",
    Mode                = "GWAS",
    Genes_Tested        = 20L,
    Superior_Haplotypes = 5L,
    PLINK_Prefix        = file.path(out_dir, "test"),
    Output_Dir          = out_dir,
    Annotation_Source   = "CUSTOM",
    PC_Cols_Gate        = "",
    stringsAsFactors    = FALSE))
saveRDS(x, file.path(out_dir, "pipeline_result_for_breeding.rds"))

cat("Test data written to", out_dir, "\n")
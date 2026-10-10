# ============================================================================
# inst/scripts/make_test_data.R
# Generate a synthetic PLINK trio + pipeline outputs for testing the engine.
# Run once: Rscript inst/scripts/make_test_data.R
# ============================================================================

set.seed(42)
out_dir <- file.path("inst", "extdata", "test_data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_ind  <- 100L
n_snp  <- 500L
ids    <- sprintf("IND%03d", seq_len(n_ind))
snp_ids <- sprintf("snp%04d", seq_len(n_snp))

# Genotypes (0/1/2 with NA)
geno <- matrix(sample(c(0L, 1L, 2L, NA_integer_), n_ind * n_snp,
                      replace = TRUE, prob = c(0.4, 0.1, 0.4, 0.1)),
               n_ind, n_snp, dimnames = list(ids, snp_ids))

# .raw file
raw <- data.frame(FID = ids, IID = ids, PAT = 0, MAT = 0, SEX = 0,
                  PHENOTYPE = -9, check.names = FALSE)
raw <- cbind(raw, as.data.frame(geno, check.names = FALSE))
write.table(raw, file.path(out_dir, "test.raw"),
            sep = " ", row.names = FALSE, quote = FALSE, na = "NA")

# .bim and .fam (needed by BEDMatrix)
bim <- data.frame(CHR = 1L, SNP = snp_ids, CM = 0,
                  BP = seq(1000L, by = 200L, length.out = n_snp),
                  A1 = sample(c("A","C","G","T"), n_snp, replace = TRUE),
                  A2 = sample(c("A","C","G","T"), n_snp, replace = TRUE))
write.table(bim, file.path(out_dir, "test.bim"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

fam <- data.frame(FID = ids, IID = ids, PAT = 0, MAT = 0, SEX = 0, PHENO = -9)
write.table(fam, file.path(out_dir, "test.fam"),
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Create a minimal .bed for BEDMatrix (write binary PLINK format)
# For testing, we write a dummy .bed with the same genotype data
write_bed <- function(geno, path) {
  n <- nrow(geno); m <- ncol(geno)
  bytes_per_snp <- ceiling(n / 4)
  con <- file(path, "wb")
  writeBin(as.raw(c(0x6c, 0x1b, 0x01)), con)  # magic number + SNP-major
  for (j in seq_len(m)) {
    packed <- raw(bytes_per_snp)
    for (i in seq_len(n)) {
      g <- geno[i, j]
      code <- if (is.na(g)) 0L else if (g == 0L) 3L else if (g == 1L) 2L else 1L
      byte_idx <- ((i - 1L) %/% 4L) + 1L
      bit_idx  <- (i - 1L) %% 4L
      packed[byte_idx] <- as.raw(
        as.integer(packed[byte_idx]) | bitwShiftL(code, bit_idx * 2L))
    }
    writeBin(packed, con)
  }
  close(con)
}
write_bed(geno, file.path(out_dir, "test.bed"))

# Phenotype table
pheno <- data.frame(ID = ids,
                    trait1 = rnorm(n_ind, 100, 15),
                    trait2 = rnorm(n_ind, 50, 8))
write.csv(pheno, file.path(out_dir, "test_pheno.csv"), row.names = FALSE)

# Minimal RDS + CSV for the engine
# (mimics pipeline_result_for_breeding.rds structure)
gene_snps <- data.table::data.table(
  gene = rep(sprintf("gene%03d", 1:20), each = 5),
  SNP  = sprintf("snp%04d", sample(1:n_snp, 100)))
data.table::fwrite(gene_snps,
                   file.path(out_dir, "functional_snps_all_tested_genes.csv"))

x <- list(
  superior_results = lapply(1:5, function(i) {
    list(test_res = list(data = data.frame(ID = ids, y = pheno$trait1)),
         hap_data = list(haplotypes = data.frame(ID = ids,
                                                  Haplotype = sample(c("H1","H2"), n_ind, TRUE))),
         superior = list(haplotype = "H1"))
  }),
  pipeline_summary = data.frame(Trait = "trait1", Direction = "higher",
                                Mode = "GWAS", Genes_Tested = 20L,
                                Superior_Haplotypes = 5L))
saveRDS(x, file.path(out_dir, "pipeline_result_for_breeding.rds"))

cat("Test data written to", out_dir, "\n")
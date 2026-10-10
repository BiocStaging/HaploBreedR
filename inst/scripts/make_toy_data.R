# ============================================================================
# inst/scripts/make_toy_data.R
# Generates the tiny demonstration dataset used by the vignette and tests.
# Run once from the package root: Rscript inst/scripts/make_toy_data.R
# ============================================================================

set.seed(42)
out_dir <- file.path("inst", "extdata", "toy")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_ind  <- 10L    # individuals
n_snp  <- 50L    # SNPs
ids    <- sprintf("IND%02d", seq_len(n_ind))
snp_ids <- sprintf("snp%03d", seq_len(n_snp))

# ---- Genotypes (0/1/2 with NA) ----
geno <- matrix(sample(c(0L, 1L, 2L, NA_integer_),
                      n_ind * n_snp, replace = TRUE,
                      prob = c(0.35, 0.15, 0.35, 0.15)),
               n_ind, n_snp,
               dimnames = list(ids, snp_ids))

# ---- PLINK .raw (numeric dosage 0/1/2, NA) ----
raw <- data.frame(FID = ids, IID = ids,
                  PAT = 0, MAT = 0, SEX = 0, PHENOTYPE = -9,
                  check.names = FALSE)
raw <- cbind(raw, as.data.frame(geno, check.names = FALSE))
write.table(raw, file.path(out_dir, "toy.raw"),
            sep = " ", row.names = FALSE, quote = FALSE, na = "NA")

# ---- PLINK .bim, .fam ----
# PLINK .bim convention: A1 = ALT allele (column 5), A2 = REF allele (column 6).
# The pipeline reads A2 as REF and A1 as ALT. Ensure A1 != A2 so no variant
# is trivially invalid.
bases <- c("A", "C", "G", "T")
A2 <- sample(bases, n_snp, replace = TRUE)
A1 <- vapply(A2, function(b) sample(setdiff(bases, b), 1), character(1))
bim <- data.frame(CHR = 1L, SNP = snp_ids, CM = 0,
                  BP = seq(1000L, by = 200L, length.out = n_snp),
                  A1 = A1, A2 = A2)
write.table(bim, file.path(out_dir, "toy.bim"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

fam <- data.frame(FID = ids, IID = ids, PAT = 0, MAT = 0,
                  SEX = 0, PHENO = -9)
write.table(fam, file.path(out_dir, "toy.fam"),
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

# ---- PLINK .bed (SNP-major; BEDMatrix requires this) ----
# PLINK encoding per call: 00 = homozygous A1, 01 = missing,
# 10 = heterozygous, 11 = homozygous A2. Dosage 2 = two ALT (A1) copies -> code 0.
# Dosage 0 = two REF (A2) copies -> code 3.
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
write_bed(geno, file.path(out_dir, "toy.bed"))

# ---- Phenotype table ----
pheno <- data.frame(ID = ids,
                    trait = round(rnorm(n_ind, mean = 100, sd = 15), 2))
write.csv(pheno, file.path(out_dir, "toy_pheno.csv"), row.names = FALSE)

# ---- GWAS hit table ----
gwas <- data.frame(SNP = snp_ids[c(5, 20, 45)],
                   Chromosome = 1L,
                   Position = bim$BP[c(5, 20, 45)])
write.csv(gwas, file.path(out_dir, "toy_gwas.csv"), row.names = FALSE)

# ---- Tiny GFF3 ----
gff <- c(
  "##gff-version 3",
  paste("Chr1", "toy", "gene",       1000, 2500, ".", "+", ".", "ID=gene1", sep = "\t"),
  paste("Chr1", "toy", "mRNA",       1000, 2500, ".", "+", ".", "ID=gene1.1;Parent=gene1", sep = "\t"),
  paste("Chr1", "toy", "CDS",        1200, 1800, ".", "+", "0", "Parent=gene1.1", sep = "\t"),
  paste("Chr1", "toy", "five_prime_UTR", 1000, 1199, ".", "+", ".", "Parent=gene1.1", sep = "\t"),
  paste("Chr1", "toy", "three_prime_UTR", 1801, 2500, ".", "+", ".", "Parent=gene1.1", sep = "\t"))
writeLines(gff, file.path(out_dir, "toy.gff3"))

# ---- Tiny reference FASTA (a placeholder 3 kb sequence) ----
writeLines(c(">Chr1", strrep("A", 3000)),
           file.path(out_dir, "toy.fa"))

cat("Toy dataset written to", out_dir, "\n")
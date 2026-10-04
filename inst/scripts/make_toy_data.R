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

# ---- PLINK .raw (numeric dosage 0/1/2, NA) ----
geno <- matrix(sample(c(0L, 1L, 2L, NA_integer_),
                      n_ind * n_snp, replace = TRUE,
                      prob = c(0.35, 0.15, 0.35, 0.15)),
               n_ind, n_snp,
               dimnames = list(ids, snp_ids))

raw <- data.frame(FID = ids, IID = ids,
                  PAT = 0, MAT = 0, SEX = 0, PHENOTYPE = -9,
                  check.names = FALSE)
raw <- cbind(raw, as.data.frame(geno, check.names = FALSE))
write.table(raw, file.path(out_dir, "toy.raw"),
            sep = " ", row.names = FALSE, quote = FALSE, na = "NA")

# ---- PLINK .bim, .fam ----
bim <- data.frame(CHR = 1L, SNP = snp_ids, CM = 0,
                  BP = seq(1000L, by = 200L, length.out = n_snp),
                  A1 = sample(c("A", "C", "G", "T"), n_snp, replace = TRUE),
                  A2 = sample(c("A", "C", "G", "T"), n_snp, replace = TRUE))
write.table(bim, file.path(out_dir, "toy.bim"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

fam <- data.frame(FID = ids, IID = ids, PAT = 0, MAT = 0,
                  SEX = 0, PHENO = -9)
write.table(fam, file.path(out_dir, "toy.fam"),
            sep = " ", row.names = FALSE, col.names = FALSE, quote = FALSE)

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
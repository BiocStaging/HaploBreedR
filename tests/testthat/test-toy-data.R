test_that("toy dataset has every expected file", {
  toy <- system.file("extdata", "toy", package = "HaploBreedR")
  expect_true(nzchar(toy))
  expect_true(dir.exists(toy))

  files <- list.files(toy)
  expected <- c("toy.raw", "toy.bim", "toy.fam",
                "toy_pheno.csv", "toy_gwas.csv",
                "toy.gff3", "toy.fa")
  for (f in expected) expect_true(f %in% files, info = f)
})

test_that("toy files are small enough for the package tarball", {
  toy <- tools::file_path_sans_ext(
    system.file("extdata", "toy", package = "HaploBreedR"))
  toy <- system.file("extdata", "toy", package = "HaploBreedR")
  sizes <- file.size(list.files(toy, full.names = TRUE))
  expect_true(all(sizes < 100 * 1024))
  expect_true(sum(sizes) < 1024 * 1024)
})
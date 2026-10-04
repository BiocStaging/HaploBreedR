test_that("hbr_palette returns expected structure", {
  expect_length(hbr_palette("haplobreed"), 8)
  expect_length(hbr_palette("okabe_ito", n = 3), 3)
  expect_error(hbr_palette("does-not-exist"))
})

test_that("hbr_theme_names and hbr_palette_names are stable", {
  expect_true("haplobreed" %in% hbr_palette_names())
  expect_true("minimal"    %in% hbr_theme_names())
  expect_s3_class(hbr_theme("minimal"), "theme")
})

test_that("hbr_theme returns all ten themes without error", {
  for (nm in hbr_theme_names()) {
    expect_true(inherits(hbr_theme(nm), "theme"),
                label = paste0("hbr_theme('", nm, "') returns a theme"))
  }
})

test_that("hbr_crop_lookup returns the expected shape", {
  crop <- hbr_crop_lookup()
  expect_s3_class(crop, "data.table")
  expect_true(all(c("crop", "species", "assembly") %in% names(crop)))
  expect_gt(nrow(crop), 5)
})
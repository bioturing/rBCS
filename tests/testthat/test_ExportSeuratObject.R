TEMPDIR <- "../../tmp"
if (!dir.exists(TEMPDIR)) {
  dir.create(TEMPDIR)
}

test_that("general", {
  object <- readRDS('../seurat4_object.rds')
  result.file <- file.path(TEMPDIR, "test_ExportSeuratObject.bcs")
  expect_true(ExportSeuratObject(object, result.file, overwrite=TRUE))
})
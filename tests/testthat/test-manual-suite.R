manual_test_files <- list.files(
  testthat::test_path("manual"),
  pattern = "^test.*\\.[rR]$",
  full.names = TRUE
)

for (manual_test_file in manual_test_files) {
  sys.source(manual_test_file, envir = environment())
}

rm(manual_test_file, manual_test_files)

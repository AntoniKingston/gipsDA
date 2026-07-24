unit_test_files <- list.files(
  testthat::test_path("unit"),
  pattern = "^test.*\\.[rR]$",
  full.names = TRUE
)

for (unit_test_file in unit_test_files) {
  sys.source(unit_test_file, envir = environment())
}

rm(unit_test_file, unit_test_files)

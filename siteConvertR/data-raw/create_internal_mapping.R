internal_mapping <- read.table(
  "data-raw/mapping.tsv",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_cols <- c("old_ID", "old_chr", "new_chr", "new_ID")
missing_cols <- setdiff(required_cols, colnames(internal_mapping))

if (length(missing_cols) > 0) {
  stop(
    "mapping.tsv 缺少这些列: ",
    paste(missing_cols, collapse = ", ")
  )
}

usethis::use_data(
  internal_mapping,
  internal = TRUE,
  overwrite = TRUE
)

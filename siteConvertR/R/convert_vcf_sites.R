#' Convert VCF sites between old and new genome coordinates
#'
#' This function converts VCF marker IDs, chromosome names, positions,
#' and optionally VCF contig header lines between old and new genome
#' coordinate systems.
#'
#' The mapping table must contain four columns:
#' `old_ID`, `old_chr`, `new_chr`, and `new_ID`.
#'
#' @param vcf_file Path to input VCF file. It can be `.vcf` or `.vcf.gz`.
#' @param output_file Path to output VCF file. `vcfR::write.vcf()` writes gzipped VCF.
#' @param direction Conversion direction. Either `"old_to_new"` or `"new_to_old"`.
#' @param mapping_file Optional path to marker correspondence table.
#' If `NULL`, the package internal mapping table will be used.
#' @param keep_unmapped Logical. Whether to keep markers not found in the mapping table.
#' If `FALSE`, unmapped markers are removed.
#' @param update_contig_header Logical. Whether to update VCF `##contig` header lines.
#' @param sort_by_position Logical. Whether to sort variants by chromosome and position.
#'
#' @return Invisibly returns the converted `vcfR` object.
#'
#' @examples
#' \dontrun{
#' convert_vcf_sites(
#'   vcf_file = "old.vcf",
#'   output_file = "new.vcf.gz",
#'   direction = "old_to_new",
#'   mapping_file = "mapping.tsv"
#' )
#'
#' convert_vcf_sites(
#'   vcf_file = "new.vcf.gz",
#'   output_file = "old.vcf.gz",
#'   direction = "new_to_old",
#'   mapping_file = "mapping.tsv"
#' )
#' }
#'
#' @export
convert_vcf_sites <- function(
    vcf_file,
    output_file = "converted.vcf.gz",
    direction = c("old_to_new", "new_to_old"),
    mapping_file = NULL,
    keep_unmapped = FALSE,
    update_contig_header = TRUE,
    sort_by_position = TRUE
) {
  direction <- match.arg(direction)

  if (!requireNamespace("vcfR", quietly = TRUE)) {
    stop("Package 'vcfR' is required. Please install it first.", call. = FALSE)
  }

  message("Reading VCF: ", vcf_file)
  vcf <- vcfR::read.vcfR(vcf_file)

  # 读取 mapping
  if (is.null(mapping_file)) {
    if (!exists("internal_mapping", envir = asNamespace("siteConvertR"))) {
      stop(
        "No internal mapping found in this package. ",
        "Please provide mapping_file.",
        call. = FALSE
      )
    }

    mapping <- get("internal_mapping", envir = asNamespace("siteConvertR"))
    message("Using internal mapping table.")
  } else {
    mapping <- read.table(
      mapping_file,
      header = TRUE,
      sep = "\t",
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    message("Using user-provided mapping table: ", mapping_file)
  }

  # 检查 mapping 列名
  required_cols <- c("old_ID", "old_chr", "new_chr", "new_ID")
  missing_cols <- setdiff(required_cols, colnames(mapping))

  if (length(missing_cols) > 0) {
    stop(
      "Mapping table is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  # 根据方向设置 from/to
  if (direction == "old_to_new") {
    if (anyDuplicated(mapping$old_ID)) {
      dup <- unique(mapping$old_ID[duplicated(mapping$old_ID)])
      stop(
        "Mapping contains duplicated old_ID, for example: ",
        paste(head(dup, 10), collapse = ", "),
        call. = FALSE
      )
    }

    from_id <- mapping$old_ID
    to_id   <- mapping$new_ID
    to_chr  <- mapping$new_chr
    chr_prefix <- "chr"

    message("Conversion direction: old genome -> new genome")
  }

  if (direction == "new_to_old") {
    if (anyDuplicated(mapping$new_ID)) {
      dup <- unique(mapping$new_ID[duplicated(mapping$new_ID)])
      stop(
        "Mapping contains duplicated new_ID, for example: ",
        paste(head(dup, 10), collapse = ", "),
        call. = FALSE
      )
    }

    from_id <- mapping$new_ID
    to_id   <- mapping$old_ID
    to_chr  <- mapping$old_chr
    chr_prefix <- "LG"

    message("Conversion direction: new genome -> old genome")
  }

  # 匹配 VCF 中的 ID
  idx <- match(vcf@fix[, "ID"], from_id)
  matched <- !is.na(idx)

  message("VCF total variants: ", nrow(vcf@fix))
  message("Matched variants: ", sum(matched))
  message("Unmatched variants: ", sum(!matched))

  # 替换 CHROM / ID / POS
  vcf@fix[matched, "CHROM"] <- to_chr[idx[matched]]
  vcf@fix[matched, "ID"]    <- to_id[idx[matched]]
  vcf@fix[matched, "POS"]   <- sub(".*_", "", to_id[idx[matched]])

  # 是否保留无法转换的位点
  if (!keep_unmapped) {
    vcf <- vcf[matched, ]
    message("Unmatched variants removed.")
  } else {
    message("Unmatched variants kept.")
  }

  # 保证 fix 是字符矩阵
  vcf@fix <- as.matrix(vcf@fix)
  storage.mode(vcf@fix) <- "character"

  # 排序，减少 PLINK split chromosome 问题
  if (sort_by_position) {
    vcf <- sort_vcf_by_chr_pos(vcf, chr_prefix = chr_prefix)
    message("Variants sorted by chromosome and position.")
  }

  # 更新 ##contig header
  if (update_contig_header) {
    vcf <- update_vcf_contig_header(vcf, contig_prefix = chr_prefix)
    message("Updated ##contig header lines.")
  }

  # write.vcf 默认写 gzipped VCF
  vcfR::write.vcf(vcf, file = output_file)

  message("Conversion finished: ", output_file)
  message("Output is gzipped VCF. Use zcat or zless to view it.")

  invisible(vcf)
}

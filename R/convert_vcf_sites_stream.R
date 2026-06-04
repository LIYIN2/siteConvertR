#' Convert large VCF files by streaming
#'
#' This function converts VCF CHROM, POS, and ID columns between old and new
#' genome coordinate systems using line-by-line streaming. It is recommended
#' for large VCF files.
#'
#' @param vcf_file Input VCF file. Can be `.vcf` or `.vcf.gz`.
#' @param output_file Output VCF file. If it ends with `.gz`, gzipped output is written.
#' @param direction Either `"old_to_new"` or `"new_to_old"`.
#' @param mapping_file Optional mapping table. If NULL, internal mapping is used.
#' @param keep_unmapped Whether to keep variants not found in mapping.
#' @param update_contig_header Whether to replace `##contig` lines.
#' @param chunk_size Number of VCF lines processed per chunk.
#'
#' @return Invisibly returns output file path.
#' @export
convert_vcf_sites_stream <- function(
    vcf_file,
    output_file = "converted.vcf.gz",
    direction = c("old_to_new", "new_to_old"),
    mapping_file = NULL,
    keep_unmapped = FALSE,
    update_contig_header = TRUE,
    chunk_size = 100000
) {
  direction <- match.arg(direction)

  # 读取 mapping
  if (is.null(mapping_file)) {
    if (!exists("internal_mapping", envir = asNamespace("siteConvertR"))) {
      stop(
        "No internal mapping found in this package. Please provide mapping_file.",
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

  required_cols <- c("old_ID", "old_chr", "new_chr", "new_ID")
  missing_cols <- setdiff(required_cols, colnames(mapping))

  if (length(missing_cols) > 0) {
    stop(
      "Mapping table is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (direction == "old_to_new") {
    if (anyDuplicated(mapping$old_ID)) {
      stop("Mapping contains duplicated old_ID.", call. = FALSE)
    }

    from_id <- mapping$old_ID
    to_id   <- mapping$new_ID
    to_chr  <- mapping$new_chr
    contig_prefix <- "chr"

    message("Conversion direction: old genome -> new genome")
  } else {
    if (anyDuplicated(mapping$new_ID)) {
      stop("Mapping contains duplicated new_ID.", call. = FALSE)
    }

    from_id <- mapping$new_ID
    to_id   <- mapping$old_ID
    to_chr  <- mapping$old_chr
    contig_prefix <- "LG"

    message("Conversion direction: new genome -> old genome")
  }

  # named vector 查询，比反复 match 快
  id_map  <- stats::setNames(to_id, from_id)
  chr_map <- stats::setNames(to_chr, from_id)
  pos_map <- stats::setNames(sub(".*_", "", to_id), from_id)

  # 生成新的 contig header
  contigs <- unique(to_chr)
  contigs <- natural_order_chr(contigs, prefix = contig_prefix)

  contig_lengths <- sapply(contigs, function(chr) {
    ids <- to_id[to_chr == chr]
    pos <- suppressWarnings(as.integer(sub(".*_", "", ids)))
    max(pos, na.rm = TRUE)
  })

  new_contig_meta <- paste0(
    "##contig=<ID=",
    contigs,
    ",length=",
    contig_lengths,
    ">"
  )

  # 输入连接
  if (grepl("\\.gz$", vcf_file)) {
    con_in <- gzfile(vcf_file, open = "rt")
  } else {
    con_in <- file(vcf_file, open = "rt")
  }

  # 输出连接
  if (grepl("\\.gz$", output_file)) {
    con_out <- gzfile(output_file, open = "wt")
  } else {
    con_out <- file(output_file, open = "wt")
  }

  on.exit({
    close(con_in)
    close(con_out)
  }, add = TRUE)

  total_variants <- 0L
  matched_variants <- 0L
  unmatched_variants <- 0L
  header_written <- FALSE
  header_buffer <- character()

  repeat {
    lines <- readLines(con_in, n = chunk_size, warn = FALSE)

    if (length(lines) == 0) {
      break
    }

    # 处理 header
    is_header <- grepl("^#", lines)

    if (any(is_header)) {
      h <- lines[is_header]

      if (update_contig_header) {
        # 去掉旧 contig
        h <- h[!grepl("^##contig=", h)]

        # 在 #CHROM 前插入新 contig
        chrom_line_idx <- grep("^#CHROM", h)

        if (length(chrom_line_idx) == 1) {
          h <- c(
            h[seq_len(chrom_line_idx - 1)],
            new_contig_meta,
            h[chrom_line_idx:length(h)]
          )
        }
      }

      writeLines(h, con_out)
    }

    # 处理变异正文
    body <- lines[!is_header]

    if (length(body) > 0) {
      fields <- strsplit(body, "\t", fixed = TRUE)

      ids <- vapply(fields, function(x) x[3], character(1))
      new_ids <- id_map[ids]

      matched <- !is.na(new_ids)

      total_variants <- total_variants + length(ids)
      matched_variants <- matched_variants + sum(matched)
      unmatched_variants <- unmatched_variants + sum(!matched)

      if (!keep_unmapped) {
        fields <- fields[matched]
        ids <- ids[matched]
        new_ids <- new_ids[matched]
        matched <- matched[matched]
      }

      if (length(fields) > 0) {
        ids2 <- vapply(fields, function(x) x[3], character(1))

        replace_idx <- !is.na(id_map[ids2])

        for (i in seq_along(fields)) {
          old_id <- fields[[i]][3]

          if (!is.na(id_map[old_id])) {
            fields[[i]][1] <- chr_map[old_id]
            fields[[i]][2] <- pos_map[old_id]
            fields[[i]][3] <- id_map[old_id]
          }
        }

        out_lines <- vapply(
          fields,
          function(x) paste(x, collapse = "\t"),
          character(1)
        )

        writeLines(out_lines, con_out)
      }
    }

    message(
      "Processed variants: ", total_variants,
      " | matched: ", matched_variants,
      " | unmatched: ", unmatched_variants
    )
  }

  message("Finished: ", output_file)
  message("Total variants: ", total_variants)
  message("Matched variants: ", matched_variants)
  message("Unmatched variants: ", unmatched_variants)

  invisible(output_file)
}

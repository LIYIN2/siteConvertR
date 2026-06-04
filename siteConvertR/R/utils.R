natural_order_chr <- function(chr_names, prefix) {
  chr_num <- suppressWarnings(
    as.integer(gsub(paste0("^", prefix), "", chr_names))
  )

  if (all(!is.na(chr_num))) {
    chr_names[order(chr_num)]
  } else {
    sort(chr_names)
  }
}


sort_vcf_by_chr_pos <- function(vcf, chr_prefix = "chr") {
  chr <- vcf@fix[, "CHROM"]
  pos <- suppressWarnings(as.integer(vcf@fix[, "POS"]))

  chr_num <- suppressWarnings(
    as.integer(gsub(paste0("^", chr_prefix), "", chr))
  )

  if (all(!is.na(chr_num))) {
    ord <- order(chr_num, pos)
  } else {
    ord <- order(chr, pos)
  }

  vcf[ord, ]
}


update_vcf_contig_header <- function(vcf, contig_prefix = "chr") {
  old_meta <- vcf@meta

  # 删除旧的 ##contig 行
  non_contig_meta <- old_meta[!grepl("^##contig=", old_meta)]

  # 从正文 CHROM 列获取 contig
  contigs <- unique(vcf@fix[, "CHROM"])
  contigs <- natural_order_chr(contigs, prefix = contig_prefix)

  # 用当前 VCF 每条染色体最大 POS 作为 length
  # 注意：这不是参考基因组真实长度，只是当前 VCF 中最大坐标
  contig_lengths <- sapply(contigs, function(chr) {
    pos <- suppressWarnings(
      as.integer(vcf@fix[vcf@fix[, "CHROM"] == chr, "POS"])
    )

    max(pos, na.rm = TRUE)
  })

  new_contig_meta <- paste0(
    "##contig=<ID=",
    contigs,
    ",length=",
    contig_lengths,
    ">"
  )

  # 插回 header，放在 ##fileformat 后面
  fileformat_line <- grep("^##fileformat=", non_contig_meta)

  if (length(fileformat_line) == 1) {
    vcf@meta <- c(
      non_contig_meta[1:fileformat_line],
      new_contig_meta,
      non_contig_meta[-(1:fileformat_line)]
    )
  } else {
    vcf@meta <- c(new_contig_meta, non_contig_meta)
  }

  vcf
}

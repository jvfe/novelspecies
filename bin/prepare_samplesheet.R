#!/usr/bin/env Rscript
#
# Merge a minimal samplesheet (sample, fasta) with GTDB-Tk taxonomy and write a
# novelspecies-ready samplesheet (sample, fasta, genus).
#
# Rows are dropped when:
#   - sample is missing from either input file
#   - taxonomy is unclassified or has no parseable GTDB genus (g__...)
#   - genus is empty / Unknown
#   - FASTA path is missing or has an unsupported extension
#   - (optional) FASTA file does not exist on disk
#
# Usage:
#   Rscript bin/prepare_samplesheet.R \
#     --gtdbtk 00_qc/gtdbtk_taxonomy.tsv \
#     --samplesheet samplesheet_minimal.csv \
#     --output samplesheet_full.csv \
#     --dropped dropped_samples.tsv
#

parse_args <- function(args) {
  if (length(args) == 0 || args[1] %in% c("-h", "--help")) {
    cat(
      "Usage: Rscript prepare_samplesheet.R --gtdbtk FILE --samplesheet FILE --output FILE",
      "       [--dropped FILE] [--require-fasta-exists]",
      sep = "\n"
    )
    quit(status = 0)
  }

  out <- list(
    gtdbtk = NA_character_,
    samplesheet = NA_character_,
    output = NA_character_,
    dropped = NA_character_,
    require_fasta_exists = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--require-fasta-exists") {
      out$require_fasta_exists <- TRUE
      i <- i + 1L
      next
    }
    if (i == length(args)) {
      stop("Missing value for argument: ", key)
    }
    val <- args[[i + 1L]]
    switch(
      key,
      "--gtdbtk" = out$gtdbtk <- val,
      "--samplesheet" = out$samplesheet <- val,
      "--output" = out$output <- val,
      "--dropped" = out$dropped <- val,
      stop("Unknown argument: ", key)
    )
    i <- i + 2L
  }

  if (is.na(out$gtdbtk) || is.na(out$samplesheet) || is.na(out$output)) {
    stop("Required arguments: --gtdbtk, --samplesheet, --output")
  }

  out
}


read_table_auto <- function(path) {
  first_line <- readLines(path, n = 1L, warn = FALSE)
  sep <- if (grepl("\t", first_line, fixed = TRUE)) "\t" else ","
  df <- read.table(
    path,
    header = TRUE,
    sep = sep,
    quote = "\"",
    comment.char = "",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  names(df) <- trimws(names(df))
  df
}


pick_column <- function(df, candidates, label) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0L) {
    stop(label, " must contain one of: ", paste(candidates, collapse = ", "))
  }
  hit[[1L]]
}


parse_gtdb_genus <- function(taxonomy) {
  taxonomy <- trimws(taxonomy)
  if (!nzchar(taxonomy)) {
    return(NA_character_)
  }
  if (grepl("^unclassified", taxonomy, ignore.case = TRUE)) {
    return(NA_character_)
  }

  match <- regexpr("g__[^;]+", taxonomy, perl = TRUE)
  if (match[[1L]] < 0L) {
    return(NA_character_)
  }

  genus <- sub("^g__", "", regmatches(taxonomy, match)[[1L]])
  genus <- trimws(genus)
  if (!nzchar(genus) || tolower(genus) %in% c("unknown", "na")) {
    return(NA_character_)
  }

  genus
}


valid_fasta <- function(path) {
  grepl("\\.(fasta|fa|fna|fas|seq)(\\.gz)?$", path, ignore.case = TRUE, perl = TRUE)
}


valid_genus <- function(genus) {
  grepl("^[A-Za-z][A-Za-z0-9_-]*$", genus, perl = TRUE)
}


main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))

  gtdb <- read_table_auto(args$gtdbtk)
  sheet <- read_table_auto(args$samplesheet)

  id_col <- pick_column(gtdb, c("genome_id", "user_genome", "genome"), "GTDB taxonomy file")
  tax_col <- pick_column(
    gtdb,
    c("classification", "gtdb_taxonomy", "taxonomy"),
    "GTDB taxonomy file"
  )

  if (!all(c("sample", "fasta") %in% names(sheet))) {
    stop("Samplesheet must contain columns: sample, fasta")
  }

  gtdb$sample <- trimws(gtdb[[id_col]])
  gtdb$gtdb_taxonomy_raw <- trimws(gtdb[[tax_col]])
  gtdb$genus <- vapply(gtdb$gtdb_taxonomy_raw, parse_gtdb_genus, character(1L))

  sheet$sample <- trimws(sheet$sample)
  sheet$fasta <- trimws(sheet$fasta)

  merged <- merge(
    sheet[, c("sample", "fasta"), drop = FALSE],
    gtdb[, c("sample", "gtdb_taxonomy_raw", "genus"), drop = FALSE],
    by = "sample",
    all = FALSE,
    sort = FALSE
  )

  audit <- data.frame(
    sample = unique(c(sheet$sample, gtdb$sample)),
    stringsAsFactors = FALSE
  )
  audit <- merge(audit, sheet[, c("sample", "fasta"), drop = FALSE], by = "sample", all.x = TRUE)
  audit <- merge(
    audit,
    gtdb[, c("sample", "gtdb_taxonomy_raw", "genus"), drop = FALSE],
    by = "sample",
    all.x = TRUE
  )

  audit$drop_reason <- NA_character_

  missing_sheet <- is.na(audit$fasta) | !nzchar(audit$fasta)
  audit$drop_reason[missing_sheet] <- "missing_in_samplesheet"

  missing_gtdb <- is.na(audit$gtdb_taxonomy_raw) | !nzchar(audit$gtdb_taxonomy_raw)
  audit$drop_reason[missing_gtdb & is.na(audit$drop_reason)] <- "missing_in_gtdbtk"

  no_genus <- is.na(audit$genus) | !nzchar(audit$genus)
  audit$drop_reason[no_genus & is.na(audit$drop_reason)] <- "unclassified_or_no_genus"

  bad_ext <- !is.na(audit$fasta) & nzchar(audit$fasta) & !valid_fasta(audit$fasta)
  audit$drop_reason[bad_ext & is.na(audit$drop_reason)] <- "unsupported_fasta_extension"

  bad_genus <- !is.na(audit$genus) & nzchar(audit$genus) & !valid_genus(audit$genus)
  audit$drop_reason[bad_genus & is.na(audit$drop_reason)] <- "invalid_genus_for_pipeline"

  if (args$require_fasta_exists) {
    missing_file <- !is.na(audit$fasta) & nzchar(audit$fasta) & !file.exists(audit$fasta)
    audit$drop_reason[missing_file & is.na(audit$drop_reason)] <- "fasta_not_found"
  }

  kept <- merged
  kept <- kept[!is.na(kept$genus) & nzchar(kept$genus), , drop = FALSE]
  kept <- kept[valid_fasta(kept$fasta), , drop = FALSE]
  kept <- kept[valid_genus(kept$genus), , drop = FALSE]
  if (args$require_fasta_exists) {
    kept <- kept[file.exists(kept$fasta), , drop = FALSE]
  }

  kept <- kept[, c("sample", "fasta", "genus"), drop = FALSE]
  kept <- kept[order(kept$sample), , drop = FALSE]

  dropped <- audit[!audit$sample %in% kept$sample, , drop = FALSE]
  dropped <- dropped[order(dropped$sample), , drop = FALSE]

  write.table(
    kept,
    file = args$output,
    sep = ",",
    quote = TRUE,
    row.names = FALSE,
    col.names = TRUE
  )

  if (!is.na(args$dropped)) {
    write.table(
      dropped[, c("sample", "fasta", "gtdb_taxonomy_raw", "genus", "drop_reason"), drop = FALSE],
      file = args$dropped,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE,
      col.names = TRUE
    )
  }

  cat("Wrote ", nrow(kept), " samples to ", args$output, "\n", sep = "")
  if (!is.na(args$dropped)) {
    cat("Wrote ", nrow(dropped), " dropped rows to ", args$dropped, "\n", sep = "")
  }
  if (nrow(dropped) > 0L) {
    cat("\nDropped samples:\n")
    print(table(dropped$drop_reason, useNA = "ifany"))
  }
}

main()

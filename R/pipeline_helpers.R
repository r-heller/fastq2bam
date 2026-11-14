# Checked external commands shared by the tutorial and its failure-path tests.
.run_tool <- function(command, args = character(), stdout = "", stderr = "") {
  tool <- Sys.which(command)
  if (!nzchar(tool)) stop("Systemtool fehlt: ", command, call. = FALSE)
  status <- system2(tool, args = shQuote(as.character(args)),
                    stdout = stdout, stderr = stderr)
  if (!identical(as.integer(status), 0L)) {
    stop(command, " fehlgeschlagen (Exit-Code ", status,
         "). Eingaben und Zwischenstände bleiben erhalten.", call. = FALSE)
  }
  invisible(status)
}

# ShortRead's heatmap uses the range of bin counts as colour breaks. Equal
# counts produce duplicate breaks and a rendered error instead of a plot.
# Show the same complete count table as points when that range is constant.
.constant_count_quality_plot <- function(qa_result) {
  counts <- qa_result[["perCycle"]][["quality"]]
  if (!nrow(counts) || length(unique(counts$Count)) != 1L) return(NULL)
  ggplot2::ggplot(counts, ggplot2::aes(x = Cycle, y = Score, size = Count)) +
    ggplot2::geom_point(colour = "#4393C3") +
    ggplot2::facet_wrap(~lane) +
    ggplot2::labs(x = "Cycle", y = "Phred quality", size = "Base count",
                  title = "Base quality by cycle",
                  caption = "Equal bin counts: points show the complete quality-count table.") +
    ggplot2::theme_minimal(base_size = 11)
}

.update_constant_quality_caption <- function(dest) {
  path <- file.path(dest, "index.html")
  html <- paste(readLines(path, warn = FALSE), collapse = "\n")
  pattern <- paste0("(?s)<p>\\s*The reddish lines are quartiles.*?</p>\\s*",
                    "<pre>\\s*perCycle.*?ShortRead:::\\.plotCycleQuality",
                    "\\(perCycle\\$quality\\)\\s*</pre>")
  matches <- gregexpr(pattern, html, perl = TRUE)[[1]]
  if (length(matches) != 1L || matches[1] < 0L) {
    stop("ShortRead-Berichtsformat geändert: Beschreibung der Qualitätsgrafik prüfen.")
  }
  replacement <- paste0(
    "<p>Blue points show every cycle and Phred quality bin. Point size gives ",
    "the number of bases in that bin. All bin counts are equal; the complete ",
    "original count table is shown.</p>\n",
    "<pre>.constant_count_quality_plot(qa_result)</pre>")
  writeLines(sub(pattern, replacement, html, perl = TRUE), path)
  invisible(TRUE)
}

.require_output <- function(path) {
  if (!file.exists(path) || dir.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0) {
    stop("Erwartete Ausgabedatei fehlt oder ist leer: ", path, call. = FALSE)
  }
  invisible(path)
}

.hisat_index_complete <- function(prefix) {
  any(vapply(c("ht2", "ht2l"), function(ext) {
    paths <- paste0(prefix, ".", 1:8, ".", ext)
    all(file.exists(paths)) && all(file.info(paths)$size > 0)
  }, logical(1)))
}

# HISAT2's Perl wrapper constructs a second shell command from its arguments.
# Give that shell only fixed relative paths, including for quoted run/input names.
# Links avoid duplicating large references or reads; originals remain untouched.
.align_hisat2 <- function(index_prefix, fastq, alignment_dir, threads, type) {
  alignment_dir <- normalizePath(alignment_dir, mustWork = TRUE)
  fastq <- normalizePath(fastq, mustWork = TRUE)
  index_prefix <- file.path(normalizePath(dirname(index_prefix), mustWork = TRUE),
                            basename(index_prefix))
  if (!.hisat_index_complete(index_prefix)) stop("HISAT2-Index ist unvollständig.")
  small <- paste0(index_prefix, ".", 1:8, ".ht2")
  ext <- if (all(file.exists(small)) && all(file.info(small)$size > 0)) "ht2" else "ht2l"
  index_files <- paste0(index_prefix, ".", 1:8, ".", ext)
  if (!all(file.exists(index_files))) stop("HISAT2-Index ist unvollständig.")
  scratch <- tempfile("hisat-paths-", alignment_dir)
  if (!dir.create(scratch)) stop("HISAT2-Arbeitsordner konnte nicht erstellt werden.")
  old_wd <- getwd()
  on.exit({ setwd(old_wd); unlink(scratch, recursive = TRUE) }, add = TRUE)
  read_alias <- if (grepl("\\.gz$", fastq, ignore.case = TRUE)) "reads.fastq.gz" else "reads.fastq"
  targets <- c(index_files, fastq, alignment_dir)
  aliases <- file.path(scratch, c(paste0("index.", 1:8, ".", ext), read_alias, "output"))
  if (!all(file.symlink(targets, aliases))) {
    stop("HISAT2-Pfadverknüpfungen fehlgeschlagen; Eingaben bleiben erhalten.")
  }
  setwd(scratch)
  args <- c("-x", "index", "-U", read_alias, "-S", "output/aligned.sam",
            "--threads", threads, "--summary-file", "output/hisat2_summary.txt")
  if (type == "dna") args <- c(args, "--no-spliced-alignment")
  .run_tool("hisat2", args)
}

# A small matched dataset for verifying file processing without organism claims.
.create_demo <- function(fastq_dir, reference_dir) {
  set.seed(42)
  reference <- paste(sample(c("A", "C", "G", "T"), 20000, replace = TRUE), collapse = "")
  starts <- seq(101L, 18020L, length.out = 100L)
  starts <- as.integer(starts)
  ids <- sprintf("synthetic_%03d", seq_along(starts))
  fasta <- file.path(reference_dir, "synthetic.fasta")
  fastq <- file.path(fastq_dir, "synthetic.fastq.gz")
  if (any(file.exists(c(fasta, fastq)))) stop("Synthetic input destination already exists.")
  offsets <- seq.int(1L, nchar(reference), by = 80L)
  writeLines(c(">synthetic_reference", substring(reference, offsets, pmin(offsets + 79L, nchar(reference)))), fasta)
  con <- gzfile(fastq, "wt")
  on.exit(close(con))
  for (i in seq_along(starts)) {
    writeLines(c(paste0("@", ids[i]), substr(reference, starts[i], starts[i] + 74L),
                 "+", paste(rep("?", 75L), collapse = "")), con)
  }
  truth <- data.frame(read_id = ids, reference = "synthetic_reference", position = starts,
                      width = 75L, stringsAsFactors = FALSE)
  write.table(truth, file.path(reference_dir, "expected_alignments.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  list(fastq = fastq, reference = fasta, truth = truth)
}


.require_fresh_directory <- function(path) {
  link <- Sys.readlink(path)
  if (file.exists(path) || (!is.na(link) && nzchar(link))) {
    stop("Ausgabeordner existiert bereits; neuen FASTQ2BAM_RUN_DIR wählen.", call. = FALSE)
  }
  invisible(path)
}


.rsubread_index_complete <- function(prefix) {
  required <- paste0(prefix, c(".files", ".reads", ".00.b.array", ".00.b.tab"))
  all(file.exists(required)) && all(file.info(required)$size > 0)
}

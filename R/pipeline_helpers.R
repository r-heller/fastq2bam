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

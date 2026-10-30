source("R/pipeline_helpers.R")
root <- tempfile("matched-demo-")
dir.create(root)
fastq_dir <- file.path(root, "fastq"); dir.create(fastq_dir)
reference_dir <- file.path(root, "reference"); dir.create(reference_dir)
demo <- .create_demo(fastq_dir, reference_dir)
reference_lines <- readLines(demo$reference)
stopifnot(all(nchar(reference_lines[-1]) <= 80L))
reference <- paste(reference_lines[-1], collapse = "")
reads <- readLines(gzfile(demo$fastq))
stopifnot(length(reads) == 400L, nrow(demo$truth) == 100L,
          length(unique(demo$truth$read_id)) == 100L)
for (i in seq_len(nrow(demo$truth))) {
  stopifnot(reads[(i-1)*4+1] == paste0("@", demo$truth$read_id[i]),
            reads[(i-1)*4+2] == substr(reference, demo$truth$position[i], demo$truth$position[i]+74L),
            nchar(reads[(i-1)*4+4]) == 75L)
}
stopifnot(inherits(try(.create_demo(fastq_dir, reference_dir), silent=TRUE), "try-error"))
unlink(root, recursive=TRUE)
cat("Synthetic read identities, positions, qualities and no-overwrite checks passed.\n")

source("R/pipeline_helpers.R")
scratch <- tempfile("command tests "); dir.create(scratch)
script <- file.path(scratch, "tool.sh")
output <- file.path(scratch, "literal argument.txt")
sentinel <- file.path(scratch, "must-not-exist")
writeLines(c("#!/bin/sh", 'printf "%s" "$2" > "$1"'), script)
literal <- paste0("spaces ; $(touch ", shQuote(sentinel), ") 'quotes'")
.run_tool("sh", c(script, output, literal))
stopifnot(identical(readChar(output, file.info(output)$size), literal), !file.exists(sentinel))
writeLines(c("#!/bin/sh", "exit 23"), script)
error <- tryCatch({.run_tool("sh", script); NULL}, error = conditionMessage)
stopifnot(grepl("Exit-Code 23", error), file.exists(output))
stopifnot(inherits(try(.require_output(file.path(scratch, "missing")), silent=TRUE), "try-error"))
prefix <- file.path(scratch, "index")
for (i in 1:7) writeLines("part", paste0(prefix, ".", i, ".ht2"))
stopifnot(!.hisat_index_complete(prefix))
writeLines("part", paste0(prefix, ".8.ht2"))
stopifnot(.hisat_index_complete(prefix))
file.create(paste0(prefix, ".8.ht2"))
stopifnot(!.hisat_index_complete(prefix))
unlink(scratch, recursive=TRUE)
cat("Command quoting, failure retention and complete-index checks passed.\n")

# Exercise link cleanup and working-directory restoration without requiring HISAT2.
# The integration workflow separately checks the actual wrapper and aligner.
local({
  root <- tempfile("alignment's input "); dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  alignment <- file.path(root, "output's directory"); dir.create(alignment)
  prefix <- file.path(root, "reference's index")
  for (i in 1:8) writeLines("index fixture", paste0(prefix, ".", i, ".ht2l"))
  # A partial small index must not override a complete large index.
  file.create(paste0(prefix, ".", 1:8, ".ht2"))
  reads <- file.path(root, "sample's reads.fastq.gz")
  writeLines("read fixture", reads)
  before <- tools::md5sum(c(reads, paste0(prefix, ".", 1:8, ".ht2l")))
  original_wd <- getwd()
  original_tool <- .run_tool
  on.exit(assign(".run_tool", original_tool, envir = .GlobalEnv), add = TRUE)
  for (fail in c(FALSE, TRUE)) {
    assign(".run_tool", function(command, args, ...) {
      stopifnot(command == "hisat2", !any(grepl("[' ]", args)),
                identical(readLines("reads.fastq.gz"), "read fixture"),
                .hisat_index_complete("index"),
                "--no-spliced-alignment" %in% args)
      writeLines("retained SAM fixture", "output/aligned.sam")
      if (fail) stop("simulated aligner failure")
      invisible(0L)
    }, envir = .GlobalEnv)
    result <- try(.align_hisat2(prefix, reads, alignment, 1L, "dna"), silent = TRUE)
    stopifnot(inherits(result, "try-error") == fail,
              identical(getwd(), original_wd),
              identical(tools::md5sum(names(before)), before),
              identical(readLines(file.path(alignment, "aligned.sam")), "retained SAM fixture"),
              !length(list.files(alignment, pattern = "^hisat-paths-")))
  }
})
cat("HISAT2 relative paths, large index selection and failure cleanup passed.\n")

new <- tempfile("new-output-")
stopifnot(identical(.require_fresh_directory(new), new))
dir.create(new)
stopifnot(inherits(try(.require_fresh_directory(new), silent = TRUE), "try-error"))
unlink(new, recursive = TRUE)
link <- tempfile("broken-output-")
if (file.symlink(paste0(new, "-missing"), link)) {
  stopifnot(inherits(try(.require_fresh_directory(link), silent = TRUE), "try-error"))
  unlink(link)
}

index_root <- tempfile("rsubread-index-"); dir.create(index_root)
prefix <- file.path(index_root, "reference")
writeLines("failed build", paste0(prefix, ".log"))
stopifnot(!.rsubread_index_complete(prefix))
for (suffix in c(".files", ".reads", ".00.b.array", ".00.b.tab")) writeLines("fixture", paste0(prefix, suffix))
stopifnot(.rsubread_index_complete(prefix))
writeLines(character(), paste0(prefix, ".00.b.tab"))
stopifnot(!.rsubread_index_complete(prefix))
unlink(index_root, recursive = TRUE)

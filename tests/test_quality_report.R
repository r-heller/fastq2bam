source("R/pipeline_helpers.R")

# Equal histogram counts can occur with several quality scores, not just Q30.
counts <- expand.grid(Cycle = 1:3, Score = c(20, 30), lane = c("A", "B"))
counts$Count <- 5L
qa <- list(perCycle = list(quality = counts))
plot <- .constant_count_quality_plot(qa)
stopifnot(inherits(plot, "ggplot"), identical(plot$data, counts))
points <- ggplot2::ggplot_build(plot)$data[[1]]
stopifnot(nrow(points) == nrow(counts),
          identical(points$x, as.numeric(counts$Cycle)),
          identical(points$y, as.numeric(counts$Score)))

qa$perCycle$quality$Count[1] <- 6L
stopifnot(is.null(.constant_count_quality_plot(qa)))
qa$perCycle$quality <- counts[FALSE, ]
stopifnot(is.null(.constant_count_quality_plot(qa)))
cat("Constant-count quality panels preserve every bin; ordinary panels are unchanged.\n")

local({
  root <- tempfile("quality-caption-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  path <- file.path(root, "index.html")
  writeLines(c("<h2>Keep other results</h2>", "<p>",
               "The reddish lines are quartiles (solid: median, dotted: 25, 75), the green line",
               "is the mean. Shading is proporitional to number of reads.",
               "</p>", "<pre>", '  perCycle &lt;- qa[["perCycle"]]',
               "  ShortRead:::.plotCycleQuality(perCycle$quality)", "</pre>",
               '<a href="./image/perCycleQuality.pdf">Quality plot</a>'), path)
  .update_constant_quality_caption(root)
  html <- paste(readLines(path), collapse = "\n")
  stopifnot(grepl("Blue points", html, fixed = TRUE),
            !grepl("quartiles|Shading|ShortRead:::", html),
            grepl("<h2>Keep other results</h2>", html, fixed = TRUE),
            grepl('./image/perCycleQuality.pdf', html, fixed = TRUE))
  before <- tools::md5sum(path)
  stopifnot(inherits(try(.update_constant_quality_caption(root), silent = TRUE), "try-error"),
            identical(tools::md5sum(path), before))
})
cat("Fallback HTML explanation matches the plot; unknown templates stop without rewriting.\n")

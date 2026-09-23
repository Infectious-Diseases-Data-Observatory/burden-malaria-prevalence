#!/usr/bin/env Rscript
# Run first: catalogues every MICS .sav file (names, labels) into data/derived_mics/ (git-ignored).
# Variable-level metadata for every MICS .sav file (names, labels, value-label counts).
suppressPackageStartupMessages(library(haven))
root <- "data/MICS_extracted"
rows <- list(); k <- 0
for (s in sort(list.dirs(root, recursive = FALSE, full.names = FALSE))) {
  for (f in sort(list.files(file.path(root, s), pattern = "\\.sav$", full.names = TRUE))) {
    d <- tryCatch(read_sav(f, n_max = 0), error = function(e) NULL)
    if (is.null(d)) { k <- k + 1; rows[[k]] <- data.frame(survey = s, file = basename(f), variable = NA,
      label = "READ_ERROR", n_value_labels = NA); next }
    lab <- vapply(d, function(x) { l <- attr(x, "label"); if (is.null(l)) "" else as.character(l)[1] }, "")
    nv  <- vapply(d, function(x) length(attr(x, "labels")), 0L)
    k <- k + 1
    rows[[k]] <- data.frame(survey = s, file = basename(f), variable = names(d), label = unname(lab),
                            n_value_labels = unname(nv), stringsAsFactors = FALSE)
  }
}
m <- do.call(rbind, rows)
write.csv(m, "data/derived_mics/variable_metadata.csv", row.names = FALSE)
cat("variables catalogued:", nrow(m), "| files:", nrow(unique(m[, c("survey","file")])),
    "| read errors:", sum(m$label == "READ_ERROR", na.rm = TRUE), "\n")

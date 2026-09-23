#!/usr/bin/env Rscript
# Step 1: extract the SPSS files of each MICS download, data/MICS_Datasets/<survey>/**/*.zip,
# into data/MICS_extracted/<survey>/. Only .sav members are written, flat (junkpaths) and with
# lower-case names, the layout 00-11 read (bh.sav, wm.sav, hh.sav, ...). A zip that wraps another
# zip (Botswana 2000) is opened one level down. A survey whose output folder already holds .sav
# files is skipped, so an existing extraction is never overwritten. No network access.
#   Rscript R_mics/00a_extract.R          extract
#   Rscript R_mics/00a_extract.R --list   print what would be extracted; writes nothing
args <- commandArgs(trailingOnly = TRUE)
if (!all(args %in% "--list")) stop("Usage: Rscript R_mics/00a_extract.R [--list]")
dry <- "--list" %in% args
src <- "data/MICS_Datasets"; dst <- "data/MICS_extracted"
if (!dir.exists(src)) stop(src, " not found; run from the project root with the MICS downloads in place.")
# Member names in older zips are CP437 bytes, not UTF-8 ("C\x93te d'Ivoire"): match on bytes.
is_sav <- function(x) grepl("[.]sav$", x, ignore.case = TRUE, useBytes = TRUE)
is_zip <- function(x) grepl("[.]zip$", x, ignore.case = TRUE, useBytes = TRUE)
leaf <- function(x) { x <- sub("^.*/", "", x, useBytes = TRUE); bad <- !validUTF8(x)
  x[bad] <- iconv(x[bad], "CP437", "UTF-8", sub = "?"); x }
members <- function(zip) { l <- utils::unzip(zip, list = TRUE)$Name
  list(sav = l[is_sav(l) & !grepl("(^|/)(__MACOSX/|[.]_)", l, useBytes = TRUE)], inner = l[is_zip(l)]) }
# Extract the .sav members of `zip` (and of any zip inside it) into fresh folders under tmp.
unpack <- function(zip, tmp, depth = 0L) {
  m <- members(zip); got <- character()
  if (length(m$sav)) {
    got <- utils::unzip(zip, files = m$sav, junkpaths = TRUE, exdir = tempfile("sav", tmp))
    if (length(got) != length(m$sav)) stop(zip, ": extracted ", length(got), " of ", length(m$sav), " .sav members")
  }
  for (z in m$inner) {
    if (depth > 0L) stop(zip, ": archives nested more than one level deep are not handled")
    inner <- utils::unzip(zip, files = z, junkpaths = TRUE, exdir = tempfile("zip", tmp))
    if (length(inner) != 1L) stop(zip, ": could not extract the nested archive ", leaf(z))
    got <- c(got, unpack(inner, tmp, depth + 1L))
  }
  got
}
surveys <- sort(list.dirs(src, recursive = FALSE, full.names = FALSE))
n_skip <- n_todo <- n_nozip <- 0L
for (s in surveys) {
  out <- file.path(dst, s)
  have <- list.files(out, "[.]sav$", ignore.case = TRUE)
  if (length(have)) { n_skip <- n_skip + 1L; cat(sprintf("skip  %s: %d .sav files already in %s\n", s, length(have), out)); next }
  zips <- list.files(file.path(src, s), "[.]zip$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (!length(zips)) { n_nozip <- n_nozip + 1L; cat(sprintf("none  %s: no zip found\n", s)); next }
  n_todo <- n_todo + 1L
  if (dry) {
    for (z in zips) {
      m <- members(z)
      cat(sprintf("list  %s <- %s\n", s, z))
      if (length(m$sav)) cat(sprintf("        %s -> %s\n", leaf(m$sav), file.path(out, tolower(leaf(m$sav)))), sep = "")
      if (length(m$inner)) cat(sprintf("        nested archive %s: its .sav members are extracted into %s\n", leaf(m$inner), out), sep = "")
      if (!length(m$sav) && !length(m$inner)) cat("        no .sav member\n")
    }
    next
  }
  tmp <- tempfile("mics_extract_"); dir.create(tmp)
  files <- unlist(lapply(zips, unpack, tmp = tmp))
  if (!length(files)) { unlink(tmp, recursive = TRUE); stop(s, ": no .sav member in ", paste(zips, collapse = ", ")) }
  target <- tolower(leaf(files))
  if (anyDuplicated(target)) stop(s, ": two .sav members share the name ", target[duplicated(target)][1])
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(files, file.path(out, target), overwrite = FALSE)
  unlink(tmp, recursive = TRUE)
  if (!all(ok)) stop(s, ": could not write ", paste(target[!ok], collapse = ", "), " into ", out)
  cat(sprintf("done  %s: %d .sav files -> %s\n", s, length(target), out))
}
cat(sprintf("surveys: %d | already extracted (skipped): %d | %s: %d | no zip: %d\n", length(surveys), n_skip,
  if (dry) "would extract" else "extracted", n_todo, n_nozip))

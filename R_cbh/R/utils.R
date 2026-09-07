cbh_num <- function(x) {
  if (is.factor(x) || is.character(x)) return(suppressWarnings(as.numeric(as.character(x))))
  suppressWarnings(as.numeric(x))
}

cbh_label <- function(x) {
  lab <- attr(x, "labels")
  if (is.numeric(x) && length(lab)) {
    out <- as.character(cbh_num(x))
    hit <- match(cbh_num(x), unname(lab))
    out[!is.na(hit)] <- names(lab)[hit[!is.na(hit)]]
    return(out)
  }
  as.character(x)
}

cbh_column <- function(d, name, numeric = FALSE) {
  if (!name %in% names(d)) return(if (numeric) rep(NA_real_, nrow(d)) else rep(NA_character_, nrow(d)))
  if (numeric) cbh_num(d[[name]]) else cbh_label(d[[name]])
}

cbh_require <- function(d, columns, name) {
  missing <- setdiff(columns, names(d))
  if (length(missing)) stop(name, " is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
}

cbh_key <- function(d, columns) {
  do.call(paste, c(lapply(d[columns], as.character), sep = "\034"))
}

cbh_unique <- function(d, columns, name) {
  cbh_require(d, columns, name)
  bad <- vapply(d[columns], function(x) anyNA(x) || any(!nzchar(trimws(as.character(x)))), logical(1))
  if (any(bad) || anyDuplicated(cbh_key(d, columns))) {
    stop(name, " has missing or duplicate keys (", paste(columns, collapse = ", "), ").", call. = FALSE)
  }
  invisible(d)
}

cbh_read_csv <- function(path, required = TRUE) {
  if (!file.exists(path)) {
    if (required) stop("Missing input: ", path, call. = FALSE)
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE, na.strings = c("", "NA"))
}

cbh_cmc_year <- function(x) 1900L + (as.integer(x) - 1L) %/% 12L
cbh_cmc_time <- function(x) 1900 + (x - 1) / 12

cbh_hash <- function(x) {
  if (!requireNamespace("digest", quietly = TRUE)) stop("The installed digest package is required; no packages are installed automatically.")
  digest::digest(x, algo = "sha256")
}

cbh_file_hash <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  unname(tools::md5sum(path))
}

cbh_atomic_rds <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp <- tempfile(".cbh-", tmpdir = dirname(path))
  on.exit(unlink(temp), add = TRUE)
  saveRDS(object, temp, compress = "gzip")
  if (!file.rename(temp, path)) stop("Could not finalize generated RDS.")
}

cbh_atomic_csv <- function(object, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  temp <- tempfile(".cbh-", tmpdir = dirname(path))
  on.exit(unlink(temp), add = TRUE)
  write.csv(object, temp, row.names = FALSE, na = "")
  if (!file.rename(temp, path)) stop("Could not finalize generated CSV.")
}

cbh_bind <- function(rows) {
  rows <- Filter(function(x) !is.null(x), rows)
  if (!length(rows)) return(data.frame())
  template <- rows[[1]][0, , drop = FALSE]
  rows <- Filter(function(x) nrow(x) > 0, rows)
  if (!length(rows)) return(template)
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Read local recodes only. Select a whitelist immediately and never print rows.
cbh_read_recode <- function(path, extra = character()) {
  if (grepl("[.]rds$", path, ignore.case = TRUE)) {
    d <- readRDS(path)
  } else if (grepl("[.]dta$", path, ignore.case = TRUE)) {
    if (!requireNamespace("haven", quietly = TRUE)) stop("Reading .dta requires the installed haven package.")
    d <- as.data.frame(haven::read_dta(path))
  } else stop("Supported local recode formats are .rds and .dta.")
  if (!is.data.frame(d)) stop("The recode is not a data frame.")
  names(d) <- tolower(names(d))
  if (anyDuplicated(names(d))) stop("Duplicate recode column names.")
  keep <- c("caseid", "v001", "v002", "v003", "v005", "v008", "v011", "v021", "v022",
            "v023", "v024", "v025", "v133", "v190", "bidx", "bord", "b0", "b3", "b4",
            "b5", "b6", "b7", "b10", "b11", "b13", extra)
  d[intersect(unique(keep), names(d))]
}

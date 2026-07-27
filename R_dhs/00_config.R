# =============================================================================
# 00_config.R
# Shared configuration and functions for the rebuilt DHS/MIS-only analysis.
#
# Run all scripts from the repository root. Raw/access-controlled DHS records
# remain under data/ and are never written to results/. Only aggregate
# survey-region outputs and non-disclosive summaries are logged.
# =============================================================================

options(stringsAsFactors = FALSE)

REPO_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!file.exists(file.path(REPO_ROOT, "R_dhs", "00_config.R"))) {
  stop("Run the rebuilt pipeline from the repository root.")
}

DATA_DIR       <- file.path(REPO_ROOT, "data")
DERIVED_DIR    <- file.path(DATA_DIR, "derived_dhs")
RESULTS_DIR    <- file.path(REPO_ROOT, "results", "dhs_rebuild")
DHS_LOCAL_DIR  <- file.path(DATA_DIR, "dhs")
DHS_CACHE_DIR  <- path.expand("~/.rdhs_cache/datasets_reformatted")
BOUNDARY_DIR   <- file.path(DATA_DIR, "dhs_boundaries")
MAP_RASTER_DIR <- file.path(DATA_DIR, "map_annual")
MAP_REGION_DIR <- file.path(DATA_DIR, "map_region_cache")
UNICEF_GLOBAL_CSV <- file.path(
  DATA_DIR, "fusion_GLOBAL_DATAFLOW_UNICEF_1.0_all.csv"
)
UNICEF_IMMUNISATION_CSV <- file.path(
  DERIVED_DIR, "unicef_immunisation_country_year.csv"
)
UNICEF_IMMUNISATION_SUMMARY_CSV <- file.path(
  RESULTS_DIR, "unicef_immunisation_extraction_summary.csv"
)

for (d in c(DERIVED_DIR, RESULTS_DIR, BOUNDARY_DIR, MAP_RASTER_DIR, MAP_REGION_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

DHS_START_YEAR <- 2000L
DHS_END_YEAR   <- 2025L
PFPR_FLOOR     <- 1
MAX_MISSING    <- 0.05
AF_REFERENCE   <- 1

MAP_DATASET_ID <- "Malaria__202508_Global_Pf_Parasite_Rate"
MAP_AFRICA_EXTENT <- matrix(c(-18, -35, 52, 38), nrow = 2)
GPW_TIF <- file.path(
  DATA_DIR, "pop",
  "gpw_v4_population_density_rev11_2020_2.5m.tif"
)

SURVEY_REGISTRY_CSV <- file.path(DERIVED_DIR, "survey_registry.csv")
MAP_REGION_CSV      <- file.path(DERIVED_DIR, "map_pfpr_by_survey_region.csv")
MAP_STATUS_CSV      <- file.path(RESULTS_DIR, "map_extraction_status.csv")
ANALYSIS_CSV        <- file.path(DERIVED_DIR, "dhs_analysis_dataset.csv")
ANALYSIS_RDS        <- file.path(DERIVED_DIR, "dhs_analysis_dataset.rds")
COVARIATE_CSV       <- file.path(RESULTS_DIR, "covariate_missingness.csv")
MODEL_BUNDLE_RDS    <- file.path(DERIVED_DIR, "main_model_bundle.rds")

required_packages <- function(packages) {
  missing <- packages[
    !vapply(packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing)) {
    stop(
      "Missing required packages: ", paste(missing, collapse = ", "),
      ". Install only after confirming they are on the approved package list."
    )
  }
  invisible(TRUE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

rkey <- function(x) {
  gsub(
    "[^a-z0-9]", "",
    tolower(iconv(as.character(x), "", "ASCII//TRANSLIT"))
  )
}

survey_key <- function(filename) {
  x <- toupper(sub("\\..*$", "", basename(filename)))
  paste0(substr(x, 1, 2), substr(x, 5, 8))
}

best_region_var <- function(br, target_keys, prefer = "v024") {
  target_keys <- unique(target_keys[nzchar(target_keys)])
  coverage <- function(v) {
    if (!v %in% names(br)) return(0)
    mean(target_keys %in% rkey(unique(as.character(br[[v]]))))
  }
  if (prefer %in% names(br) && coverage(prefer) >= 0.5) return(prefer)

  candidates <- names(br)[vapply(
    br,
    function(x) is.character(x) || is.factor(x),
    logical(1)
  )]
  if (!length(candidates)) return(prefer)
  scores <- vapply(
    candidates,
    function(v) sum(target_keys %in% rkey(unique(as.character(br[[v]])))),
    integer(1)
  )
  best <- candidates[which.max(scores)]
  if (max(scores) > 0 && coverage(best) >= 0.6) best else prefer
}

weighted_mean_by_region <- function(value, weight, region, eligible) {
  keep <- eligible & is.finite(value) & is.finite(weight) &
    !is.na(region) & nzchar(region)
  if (!any(keep)) return(numeric(0))
  numerator <- tapply(weight[keep] * value[keep], region[keep], sum)
  denominator <- tapply(weight[keep], region[keep], sum)
  numerator / denominator
}

pick_named <- function(x, keys) {
  if (!length(x)) return(rep(NA_real_, length(keys)))
  unname(x[keys])
}

mortality_by_region <- function(br, region_var) {
  required_packages("DHS.rates")
  if (!region_var %in% names(br)) return(NULL)
  rates <- tryCatch(
    suppressMessages(DHS.rates::chmort(br, Class = region_var)),
    error = function(e) NULL
  )
  if (is.null(rates)) return(NULL)

  u5 <- rates[
    grepl("^U5MR", rownames(rates)),
    c("Class", "R", "WN"),
    drop = FALSE
  ]
  nn <- rates[
    grepl("^NNMR", rownames(rates)),
    c("Class", "R"),
    drop = FALSE
  ]
  if (!nrow(u5) || !nrow(nn)) return(NULL)
  names(u5)[2:3] <- c("u5mr", "exposure")
  names(nn)[2] <- "nnmr"
  out <- merge(u5, nn, by = "Class")
  data.frame(
    regkey = rkey(out$Class),
    u5mr = as.numeric(out$u5mr),
    nnmr = as.numeric(out$nnmr),
    postneonatal_mortality = as.numeric(out$u5mr - out$nnmr),
    exposure = as.numeric(out$exposure)
  )
}

local_recode_path <- function(no_extension) {
  candidates <- c(
    file.path(DHS_LOCAL_DIR, paste0(no_extension, ".rds")),
    file.path(DHS_CACHE_DIR, paste0(no_extension, ".rds"))
  )
  hit <- candidates[file.exists(candidates)]
  if (length(hit)) normalizePath(hit[1], winslash = "/") else NA_character_
}

wb_fetch <- function(indicator, date = "2000:2025") {
  required_packages(c("httr", "jsonlite"))
  url <- sprintf(
    paste0(
      "https://api.worldbank.org/v2/country/all/indicator/%s",
      "?date=%s&format=json&per_page=20000"
    ),
    indicator, date
  )
  for (attempt in seq_len(4)) {
    response <- tryCatch(
      httr::GET(url, httr::timeout(180)),
      error = function(e) NULL
    )
    if (!is.null(response) && httr::status_code(response) == 200) {
      payload <- jsonlite::fromJSON(
        httr::content(response, "text", encoding = "UTF-8"),
        simplifyDataFrame = TRUE
      )[[2]]
      out <- data.frame(
        iso3 = payload$countryiso3code,
        year = as.integer(payload$date),
        value = as.numeric(payload$value)
      )
      return(out[is.finite(out$value) & nchar(out$iso3) == 3, ])
    }
    message("World Bank retry ", attempt, " for ", indicator)
  }
  stop("World Bank fetch failed: ", indicator)
}

read_or_fetch_wb <- function(filename, indicator, value_name) {
  path <- file.path(DATA_DIR, filename)
  if (file.exists(path)) {
    out <- read.csv(path, stringsAsFactors = FALSE)
  } else {
    message("Fetching World Bank indicator ", indicator, " -> ", filename)
    out <- wb_fetch(indicator)
    names(out)[names(out) == "value"] <- value_name
    write.csv(out, path, row.names = FALSE)
  }
  if (!value_name %in% names(out) && "value" %in% names(out)) {
    names(out)[names(out) == "value"] <- value_name
  }
  out
}

nearest_panel_value <- function(panel, iso3, year, value_name) {
  rows <- panel[
    panel$iso3 == iso3 & is.finite(panel[[value_name]]),
    ,
    drop = FALSE
  ]
  if (!nrow(rows)) return(NA_real_)
  rows[[value_name]][which.min(abs(rows$year - year))]
}

single_impute <- function(x, country) {
  x[!is.finite(x)] <- NA_real_
  country_median <- ave(
    x, country,
    FUN = function(z) {
      value <- suppressWarnings(median(z, na.rm = TRUE))
      if (is.finite(value)) value else NA_real_
    }
  )
  overall <- suppressWarnings(median(x, na.rm = TRUE))
  out <- x
  missing <- !is.finite(out)
  out[missing] <- country_median[missing]
  out[!is.finite(out)] <- overall
  out
}

apply_missingness_rule <- function(data, catalog, threshold = MAX_MISSING) {
  stopifnot(all(c("variable", "level", "type") %in% names(catalog)))
  rows <- vector("list", nrow(catalog))

  for (i in seq_len(nrow(catalog))) {
    variable <- catalog$variable[i]
    if (!variable %in% names(data)) data[[variable]] <- NA_real_
    raw <- suppressWarnings(as.numeric(data[[variable]]))
    raw[!is.finite(raw)] <- NA_real_
    missing <- !is.finite(raw)
    missing_prop <- mean(missing)
    include <- is.finite(missing_prop) && missing_prop <= threshold &&
      any(is.finite(raw))

    analysis_name <- paste0(variable, "_analysis")
    imputed_name <- paste0(variable, "_imputed")
    data[[imputed_name]] <- missing
    data[[analysis_name]] <- if (include) {
      single_impute(raw, data$iso3)
    } else {
      rep(NA_real_, nrow(data))
    }

    rows[[i]] <- data.frame(
      variable = variable,
      level = catalog$level[i],
      type = catalog$type[i],
      missing_n = sum(missing),
      total_n = length(raw),
      missing_prop = missing_prop,
      included_in_main = include,
      imputation = if (include && any(missing)) {
        "country median; overall median fallback"
      } else if (include) {
        "none required"
      } else {
        "not imputed; excluded"
      }
    )
  }
  list(data = data, catalog = do.call(rbind, rows))
}

logit_percent <- function(x) {
  p <- pmin(pmax(x / 100, 0.005), 0.995)
  log(p / (1 - p))
}

make_ridge_matrix <- function(data, catalog, preprocessing = NULL) {
  included <- catalog$variable[as.logical(catalog$included_in_main)]
  if (!length(included)) stop("No covariates passed the missingness rule.")

  transformed <- lapply(included, function(variable) {
    x <- data[[paste0(variable, "_analysis")]]
    type <- catalog$type[match(variable, catalog$variable)]
    if (identical(type, "proportion")) logit_percent(x) else as.numeric(x)
  })
  names(transformed) <- included

  if (is.null(preprocessing)) {
    means <- vapply(transformed, mean, numeric(1), na.rm = TRUE)
    sds <- vapply(transformed, stats::sd, numeric(1), na.rm = TRUE)
    keep <- is.finite(means) & is.finite(sds) & sds > 0
    means <- means[keep]
    sds <- sds[keep]
    included <- names(means)
    preprocessing <- list(
      variables = included,
      means = means,
      sds = sds,
      types = setNames(
        catalog$type[match(included, catalog$variable)],
        included
      )
    )
  } else {
    included <- preprocessing$variables
  }

  matrix_out <- vapply(included, function(variable) {
    x <- data[[paste0(variable, "_analysis")]]
    if (identical(preprocessing$types[[variable]], "proportion")) {
      x <- logit_percent(x)
    }
    (x - preprocessing$means[[variable]]) / preprocessing$sds[[variable]]
  }, numeric(nrow(data)))
  matrix_out <- as.matrix(matrix_out)
  colnames(matrix_out) <- included

  list(matrix = matrix_out, preprocessing = preprocessing)
}

MODEL_SPECS <- data.frame(
  specification = c(
    "linear_no_interaction",
    "spline_no_interaction",
    "linear_time_interaction",
    "spline_time_interaction",
    "full_te_surface"
  ),
  prevalence_form = c("linear", "spline", "linear", "spline", "tensor"),
  time_interaction = c(FALSE, FALSE, TRUE, TRUE, TRUE)
)

model_formula <- function(
    specification,
    include_country_slope = TRUE,
    spline_k = 6,
    year_k = 8) {
  prevalence <- switch(
    specification,
    linear_no_interaction = "pfpr10",
    spline_no_interaction = sprintf("s(pfpr10, k=%d)", spline_k),
    linear_time_interaction = "pfpr10 + pfpr10:year_c",
    spline_time_interaction = sprintf(
      "s(pfpr10, k=%d) + ti(pfpr10, year_c, k=c(%d,%d))",
      spline_k, spline_k, min(year_k, 6)
    ),
    full_te_surface = sprintf(
      "te(pfpr10, year_c, k=c(%d,%d))",
      spline_k, min(year_k, 6)
    ),
    stop("Unknown model specification: ", specification)
  )
  random_terms <- c(
    "s(country, bs='re')",
    if (include_country_slope) "s(country, pfpr10, bs='re')"
  )
  terms <- c(
    prevalence,
    "G",
    if (!identical(specification, "full_te_surface")) {
      sprintf("s(year_c, k=%d)", year_k)
    },
    random_terms,
    "offset(log(exposure))"
  )
  as.formula(paste("deaths ~", paste(terms, collapse = " + ")))
}

fit_ridge_gam <- function(
    data,
    outcome,
    catalog,
    specification,
    method = "ML",
    preprocessing = NULL,
    include_country_slope = TRUE,
    spline_k = 6,
    year_k = 8) {
  required_packages("mgcv")
  keep <- is.finite(data[[outcome]]) & data[[outcome]] > 0 &
    is.finite(data$exposure) & data$exposure > 0 &
    is.finite(data$pfpr10) & is.finite(data$year_c) &
    !is.na(data$iso3)
  dd <- data[keep, , drop = FALSE]
  dd$country <- factor(dd$iso3)
  dd$deaths <- round(dd[[outcome]] / 1000 * dd$exposure)

  ridge <- make_ridge_matrix(dd, catalog, preprocessing)
  dd$G <- ridge$matrix
  penalty <- list(G = list(diag(ncol(ridge$matrix))))

  model <- mgcv::gam(
    model_formula(
      specification,
      include_country_slope = include_country_slope,
      spline_k = spline_k,
      year_k = year_k
    ),
    family = mgcv::nb(),
    method = method,
    paraPen = penalty,
    data = dd
  )
  list(
    model = model,
    data = dd,
    preprocessing = ridge$preprocessing,
    specification = specification,
    outcome = outcome,
    include_country_slope = include_country_slope,
    spline_k = spline_k,
    year_k = year_k
  )
}

model_summary_row <- function(fit) {
  model <- fit$model
  summary_model <- summary(model)
  linear <- "pfpr10" %in% rownames(summary_model$p.table)
  spline_row <- grep(
    "^s\\(pfpr10\\)",
    rownames(summary_model$s.table),
    value = TRUE
  )
  tensor_row <- grep(
    "^ti\\(pfpr10,year_c\\)",
    rownames(summary_model$s.table),
    value = TRUE
  )
  full_tensor_row <- grep(
    "^te\\(pfpr10,year_c\\)",
    rownames(summary_model$s.table),
    value = TRUE
  )
  interaction_linear <- "pfpr10:year_c" %in% rownames(summary_model$p.table)
  beta <- se <- p <- edf <- interaction_beta <- interaction_se <-
    interaction_edf <- interaction_p <- full_surface_edf <-
    full_surface_p <- NA_real_
  if (linear) {
    beta <- summary_model$p.table["pfpr10", "Estimate"]
    se <- summary_model$p.table["pfpr10", "Std. Error"]
    p <- summary_model$p.table[
      "pfpr10",
      grep("^Pr\\(", colnames(summary_model$p.table), value = TRUE)[1]
    ]
  } else if (length(spline_row)) {
    edf <- summary_model$s.table[spline_row[1], "edf"]
    p <- summary_model$s.table[spline_row[1], "p-value"]
  }
  if (interaction_linear) {
    interaction_beta <- summary_model$p.table["pfpr10:year_c", "Estimate"]
    interaction_se <- summary_model$p.table["pfpr10:year_c", "Std. Error"]
    interaction_p <- summary_model$p.table[
      "pfpr10:year_c",
      grep("^Pr\\(", colnames(summary_model$p.table), value = TRUE)[1]
    ]
  } else if (length(tensor_row)) {
    interaction_edf <- summary_model$s.table[tensor_row[1], "edf"]
    interaction_p <- summary_model$s.table[tensor_row[1], "p-value"]
  }
  if (length(full_tensor_row)) {
    full_surface_edf <- summary_model$s.table[full_tensor_row[1], "edf"]
    full_surface_p <- summary_model$s.table[full_tensor_row[1], "p-value"]
  }
  data.frame(
    specification = fit$specification,
    outcome = fit$outcome,
    n = nrow(model$model),
    countries = nlevels(model$model$country),
    edf_total = sum(model$edf),
    AIC = AIC(model),
    pfpr_beta = beta,
    pfpr_se = se,
    pct_change_per_10 = if (is.finite(beta)) 100 * (exp(beta) - 1) else NA,
    pct_change_lo = if (is.finite(beta)) 100 * (exp(beta - 1.96 * se) - 1) else NA,
    pct_change_hi = if (is.finite(beta)) 100 * (exp(beta + 1.96 * se) - 1) else NA,
    pfpr_smooth_edf = edf,
    pfpr_p = p,
    time_interaction_beta = interaction_beta,
    time_interaction_se = interaction_se,
    time_interaction_edf = interaction_edf,
    time_interaction_p = interaction_p,
    full_surface_edf = full_surface_edf,
    full_surface_p = full_surface_p
  )
}

newdata_at_mean <- function(model, pfpr10, year_c = 0) {
  n <- length(pfpr10)
  country_level <- levels(model$model$country)[1]
  g_names <- colnames(model$model$G)
  out <- data.frame(
    pfpr10 = pfpr10,
    year_c = rep(year_c, length.out = n),
    exposure = 1,
    country = factor(rep(country_level, n), levels = levels(model$model$country))
  )
  out$G <- matrix(0, nrow = n, ncol = length(g_names))
  colnames(out$G) <- g_names
  out
}

population_lpmatrix <- function(model, newdata) {
  X <- predict(model, newdata, type = "lpmatrix")
  random_columns <- grep(
    "^s\\(country\\)|^s\\(country,pfpr10\\)",
    colnames(X)
  )
  if (length(random_columns)) X[, random_columns] <- 0
  X
}

link_prediction <- function(model, newdata) {
  X <- population_lpmatrix(model, newdata)
  beta <- coef(model)
  fit <- as.numeric(X %*% beta)
  se <- sqrt(rowSums((X %*% vcov(model)) * X))
  data.frame(fit = fit, se = se)
}

link_difference <- function(model, high, reference) {
  Xh <- population_lpmatrix(model, high)
  Xr <- population_lpmatrix(model, reference)
  dX <- Xh - Xr
  beta <- coef(model)
  fit <- as.numeric(dX %*% beta)
  se <- sqrt(rowSums((dX %*% vcov(model)) * dX))
  data.frame(fit = fit, se = se)
}

af_from_model <- function(model, prevalence, year_c = 0, reference = AF_REFERENCE) {
  high <- newdata_at_mean(model, prevalence / 10, year_c)
  low <- newdata_at_mean(model, rep(reference / 10, length(prevalence)), year_c)
  delta <- link_difference(model, high, low)
  data.frame(
    prevalence = prevalence,
    af = pmax(1 - exp(-delta$fit), 0),
    lo = pmax(1 - exp(-(delta$fit - 1.96 * delta$se)), 0),
    hi = pmax(1 - exp(-(delta$fit + 1.96 * delta$se)), 0)
  )
}

read_analysis_data <- function() {
  if (file.exists(ANALYSIS_RDS)) {
    readRDS(ANALYSIS_RDS)
  } else if (file.exists(ANALYSIS_CSV)) {
    out <- read.csv(ANALYSIS_CSV, stringsAsFactors = FALSE)
    for (v in intersect(
      c("main_sample", "pfpr_region_ge_1", "country_mean_pfpr_gt_1",
        "complete_case_eligible"),
      names(out)
    )) out[[v]] <- as.logical(out[[v]])
    out
  } else {
    stop("Analysis dataset not found. Run R_dhs/03_build_analysis_dataset.R.")
  }
}

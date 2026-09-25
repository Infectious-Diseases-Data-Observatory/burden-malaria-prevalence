# Shared helpers for the tensor-product sensitivity: saved PfPR/time components and the log hazard
# ratio for PfPR `from` -> `to` at calendar time `year` (vectors recycled), with its SE. Covariates,
# random effects and (in the separate and ti models) s(calendar_year) cancel in the contrast.
cbh_tensor_components <- function(private = "data/derived_cbh/models/tensor_pfpr_year_dhsmics_map_gamma2_v1")
  readRDS(file.path(private, "pfpr_year_components.rds"))
cbh_tensor_last_entry <- 2023 + 11 / 12   # last observed band entry and last calendar-year knot
cbh_tensor_eval_year <- function(year) pmin(year + .5, cbh_tensor_last_entry)
cbh_tensor_lhr <- function(z, from, to, year) {
  n <- max(length(from), length(to), length(year)); from <- rep_len(from, n); to <- rep_len(to, n); year <- rep_len(year, n)
  L <- matrix(0, n, length(z$ix))
  for (k in seq_along(z$smooths)) {
    s <- z$smooths[[k]]; cols <- match(z$index[[k]], z$ix)
    L[, cols] <- L[, cols] + mgcv::PredictMat(s, data.frame(pfpr_pct = to, calendar_year = year)) - mgcv::PredictMat(s, data.frame(pfpr_pct = from, calendar_year = year))
  }
  list(est = drop(L %*% z$coef), se = sqrt(pmax(rowSums((L %*% z$Vp) * L), 0)))
}

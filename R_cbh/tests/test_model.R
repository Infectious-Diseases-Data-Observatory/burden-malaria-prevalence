source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
set.seed(46)
n <- 7000L
spec <- cbh_trial_spec()
age_levels <- cbh_config()$age_bands$age_band
s <- sample(1:12, n, TRUE)
a <- sample(1:7, n, TRUE)
d <- data.frame(age_band = age_levels[a], survey = paste0("s", s),
  country = paste0("c", (s - 1) %/% 4), region = paste0("s", s, "r", sample(1:3, n, TRUE)),
  pfpr_pct = runif(n, 0, 80), calendar_year = runif(n, 2000, 2024),
  band_years = c(1, 5, 6, 12, 12, 12, 12)[a] / 12,
  sex = sample(c("female", "male"), n, TRUE), multiple_birth = rbinom(n, 1, .08),
  birth_order = sample(1:8, n, TRUE), maternal_age_birth = runif(n, 15, 45),
  maternal_education_years = sample(0:16, n, TRUE), wealth_quintile = sample(1:5, n, TRUE),
  urban = rbinom(n, 1, .4), log_hiv_prev = rnorm(n), log_gdp_pc = rnorm(n, 7),
  log_health_expenditure_pc = rnorm(n, 4), political_stability = rnorm(n))
eta <- -1.7 - .18 * a + .012 * d$pfpr_pct + .2 * d$multiple_birth
d$death <- rbinom(n, 1, -expm1(-d$band_years * exp(eta)))
p <- cbh_trial_prepare(d, age_levels, spec)
stopifnot(nrow(p$data) == n, !is.ordered(p$data$age_band), nlevels(p$data$wealth_quintile) == 5,
          all(abs(vapply(p$data[paste0("z_", spec$scaled)], mean, numeric(1))) < 1e-10))
bad <- d; bad$log_hiv_prev[1] <- NA
stopifnot(inherits(try(cbh_trial_prepare(bad, age_levels, spec), silent = TRUE), "try-error"))
res <- cbh_trial_fit(p$data, spec, trace = FALSE)
fit <- res$fit
stopifnot(!identical(fit$converged, FALSE), all(is.finite(coef(fit))), length(fit$fitted.values) == n,
          sum(vapply(fit$smooth, function(s) identical(s$term, "pfpr_pct"), logical(1))) == 7,
          sum(vapply(fit$smooth, function(s) identical(s$term, "calendar_year"), logical(1))) == 7)
z <- cbh_trial_contrasts(fit, p$data)
identity_contrast <- cbh_trial_contrasts(fit, p$data, from = 20, to = 20)
stopifnot(nrow(z) == 7, all(z$lower_95 <= z$hazard_ratio & z$upper_95 >= z$hazard_ratio),
          all(identity_contrast$hazard_ratio == 1), all(identity_contrast$lower_95 == 1),
          all(identity_contrast$upper_95 == 1))
cat("Synthetic model checks passed: seven age-specific PfPR/time curves, complete input, and prediction contrasts.\n")

# Adebayo & Fahrmeir (2005), Statistics in Medicine 24:709-728.
# DOI: 10.1002/sim.1842. See docs/ADEBAYO_2005_MGCV.md.
#
# Model specification only: sourcing this file reads no data and fits no models.
# Input: one row per child-month at risk, with event = 0/1 and month = 1, 2, ...
# (month 1 denotes [0, 1) months of age). Stop each history at death/censoring.
# nb: named, symmetric district adjacency list, as required by mgcv's mrf basis.
# The template assumes a connected graph; review islands/disconnected components
# before extending it to multiple countries.
#
# *_ec variables below are numeric effect codes (+1 / -1), not factors.
# bft is numeric 0/1 breastfeeding status IN THAT MONTH.
# maternal_age_birth is mother's age when the child was born, in years.
# See the companion note for the covariate dictionary and censoring convention.

adebayo_2005_formula <- function(nb, model = c("M7", "M6"),
                                 include_socioeconomic = FALSE,
                                 k_age = 10L, k_maternal = 10L) {
  model <- match.arg(model)
  # Keep the formula self-contained without attaching mgcv to the search path.
  s <- mgcv::s

  if (model == "M7") {
    form <- event ~
      s(month, bs = "ps", k = k_age, m = c(2, 2)) +
      s(month, by = bft, bs = "ps", k = k_age, m = c(2, 2)) +
      s(month, by = ma1, bs = "ps", k = k_age, m = c(2, 2)) +
      s(month, by = ma2, bs = "ps", k = k_age, m = c(2, 2)) +
      s(district, bs = "mrf", xt = list(nb = nb)) +
      s(district_iid, bs = "re") +
      urban_ec + male_ec + primary_or_less_ec + assisted_ec +
      hospital_ec + long_interval_ec + assisted_ec:antenatal_ec
  } else {
    form <- event ~
      s(month, bs = "ps", k = k_age, m = c(2, 2)) +
      s(month, by = bft, bs = "ps", k = k_age, m = c(2, 2)) +
      s(maternal_age_birth, bs = "ps", k = k_maternal, m = c(2, 2)) +
      s(district, bs = "mrf", xt = list(nb = nb)) +
      s(district_iid, bs = "re") +
      urban_ec + male_ec + primary_or_less_ec + assisted_ec +
      hospital_ec + long_interval_ec + assisted_ec:antenatal_ec
  }

  if (include_socioeconomic) {
    form <- stats::update.formula(form, . ~ . + working_ec + household_le5_ec +
                                   toilet_ec + electricity_ec +
                                   quality_house_ec + water_residence_ec)
  }
  # The formula retains nb and the basis dimensions in its environment.
  form
}

fit_adebayo_2005 <- function(person_months, nb, model = c("M7", "M6"),
                            include_socioeconomic = FALSE,
                            k_age = 10L, k_maternal = 10L) {
  model <- match.arg(model)
  d <- as.data.frame(person_months)
  required <- c("event", "month", "bft", "maternal_age_birth", "district",
                "urban_ec", "male_ec", "primary_or_less_ec", "assisted_ec",
                "hospital_ec", "long_interval_ec", "antenatal_ec")
  if (include_socioeconomic) {
    required <- c(required, "working_ec", "household_le5_ec", "toilet_ec",
                  "electricity_ec", "quality_house_ec", "water_residence_ec")
  }
  missing <- setdiff(required, names(d))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
  if (anyNA(d[required])) stop("Resolve missing model inputs before fitting.")
  if (!is.numeric(d$event) || !all(d$event %in% c(0, 1))) {
    stop("event must be numeric 0/1 for each child-month.")
  }
  if (!is.numeric(d$bft) || !all(d$bft %in% c(0, 1))) {
    stop("bft must be numeric 0/1, not a factor.")
  }
  if (!is.numeric(d$month) || any(!is.finite(d$month)) ||
      any(d$month < 1 | d$month != floor(d$month))) {
    stop("month must be an integer-valued index starting at 1.")
  }
  if (!is.numeric(d$maternal_age_birth) ||
      any(!is.finite(d$maternal_age_birth))) {
    stop("maternal_age_birth must be finite and numeric.")
  }
  ec <- required[endsWith(required, "_ec")]
  if (!all(vapply(d[ec], function(x) is.numeric(x) && all(x %in% c(-1, 1)),
                  logical(1)))) {
    stop("All *_ec columns must contain numeric -1/+1 effect codes.")
  }
  if (is.null(names(nb)) || anyDuplicated(names(nb)) ||
      !all(as.character(d$district) %in% names(nb))) {
    stop("nb must have unique district names covering all observed districts.")
  }

  # The exact maternal-age contrasts described for M7 on page 718.
  d$ma1 <- as.numeric(d$maternal_age_birth < 22) -
    as.numeric(d$maternal_age_birth > 35)
  d$ma2 <- as.numeric(d$maternal_age_birth >= 22 & d$maternal_age_birth <= 35) -
    as.numeric(d$maternal_age_birth > 35)
  d$district <- factor(as.character(d$district), levels = names(nb))
  d$district_iid <- d$district

  form <- adebayo_2005_formula(nb, model, include_socioeconomic, k_age, k_maternal)
  mgcv::gam(form, data = d, family = stats::binomial(link = "logit"),
            method = "REML", na.action = stats::na.fail,
            drop.unused.levels = FALSE,
            knots = list(district = names(nb)))
}

# Example, after explicitly constructing the required data and adjacency list:
# source("R_dhs/specifications/adebayo_2005_mgcv.R")
# fit_m7 <- fit_adebayo_2005(person_months, nb)
# fit_m6 <- fit_adebayo_2005(person_months, nb, model = "M6")
# fit_m7_extended <- fit_adebayo_2005(person_months, nb,
#                                   include_socioeconomic = TRUE)

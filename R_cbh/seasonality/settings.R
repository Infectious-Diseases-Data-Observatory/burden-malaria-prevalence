# Settings for the Sahel calendar-month seasonality analysis. No reads or fits here.
cbh_seasonality_settings <- function() {
  list(
    # Sahelian survey regions: boundary centroid at or above 11 degrees north, west of
    # the Horn of Africa. 11 rather than 12 degrees so that Nigeria's North West and
    # North East zones, whose centroids fall just below 12 degrees, are included.
    latitude_min = 11,
    longitude_max = 36,
    exclude_countries = c("ETH", "ERI", "SOM", "DJI"),
    centroids = "data/derived_dhs/dhs_region_centroids.csv",
    # months before the interview month that contribute exposure (the interview month
    # itself is partial and excluded), and the oldest completed month of age analysed
    lookback_months = 60L,
    max_age_months = 23L,
    bands = data.frame(band = c("<1 month", "1-5 months", "6-11 months", "12-23 months"),
                       age_lo = c(0L, 1L, 6L, 12L), age_hi = c(1L, 6L, 12L, 24L),
                       stringsAsFactors = FALSE),
    # surveys need at least this many deaths in a band to be drawn as their own line
    min_deaths_for_survey_line = 50,
    out = "results/cbh/seasonality_sahel_v1"
  )
}

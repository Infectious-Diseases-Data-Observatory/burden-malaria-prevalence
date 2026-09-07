# Complete-birth-history dataset pipeline

Stage 1 of the new pipeline for the [seven-band PfPR model](../docs/PFPR_AGE_BAND_MODEL.md). It reads local Births Recodes and external-data snapshots, creates child–age-band records, and checks eligibility and joins. It does not fit a model. No legacy fitting script or fitted object is required.

## Run

From the project root, using R and the already installed `digest` package (`haven` is needed only for `.dta` inputs):

```sh
Rscript R_cbh/tests/test_dataset.R
Rscript R_cbh/01_make_analysis_data.R --check-inputs
Rscript R_cbh/01_make_analysis_data.R
```

For a subset, use a separate output directory so its manifest does not replace the full build's manifest:

```sh
Rscript R_cbh/01_make_analysis_data.R --surveys=NG7BFL,ET81FL,ET8AFL --output=data/derived_cbh_validation
```

Output defaults to **`data/derived_cbh/`**. `child_bands/<svkey>.rds` contains each survey's records and audits. These are local microdata, retained under the project's ignored `data/` directory. Original mother/household/cluster identifiers are replaced by survey-specific pseudonymous identifiers. These remain research microdata, not public release files.

`manifest.rds` is the authoritative list of shards for the last completed build. Readers refuse interrupted builds and inconsistent shard signatures. Old shards not named in the manifest are ignored. A subset run represents only that subset; it does not add surveys to an earlier manifest. Do not run two builders against the same output directory concurrently.

`--force` rebuilds selected surveys. Otherwise, matching cached shards are reused. Signatures cover the raw recode (and geographic donor, if used), settings, code, registry, rules, overrides and external snapshots. Each manifest also records file hashes, build times and the R session. Missing optional panels remain missing; invalid/duplicate join keys stop the build before processing records. A failed survey is recorded and makes the command exit unsuccessfully; the reader requires explicit permission to use such a partial build.

## Eligibility and outcome

Seven completed-month bands: **<1, 1–5, 6–11, 12–23, 24–35, 36–47, 48–59**. For each child and band:

1. The child must have reached the band alive.
2. Band entry must be in the 60 months before interview: `interview_cmc - 60 <= band_entry_cmc < interview_cmc`. Birth itself need not fall within that window.
3. The full potential band must end by interview, **for deaths and survivors alike**. A death in an incomplete potential band is not included, and no band after death is included.
4. The entry calendar year must be within the configured MAP period, currently 2000–2024. Out-of-range entries are counted separately.

The binary `death` outcome is 1 for a death in that band and 0 for survival of the band. `band_years` is its predetermined width: 1/12, 5/12, 6/12 or 1 year. The cloglog offset is `log_band_years`. **This is not observed person-time until death**, and these records are not the former aggregate Poisson/negative-binomial dataset. Calendar-month precision cannot resolve exact interview/birthday ordering within a month.

Reported death-age units in B6 are used first. Month reports identify a completed-month interval; year reports identify a whole-year interval. B7 imputed months provide a flagged fallback when B6 is unusable. Reports spanning a chosen band boundary are excluded and counted. B6/B7 band disagreements are retained as flags, with B6 taking precedence; review them before fitting. Unknown outcomes, invalid dates/weights, impossible deaths and duplicate child identifiers are counted and excluded.

The neonatal definition follows the DHS month-based convention: reports of **0–30 days** enter the first band, and 31–96 days the 1–5-month band. Reported days are retained. This differs from the exact 28-day boundary in Burstein's text. See the [DHS mortality guide](https://www.dhsprogram.com/Data/Guide-to-DHS-Statistics/Early_Childhood_Mortality.htm) and the [DHS statistics manual](https://www.dhsprogram.com/pubs/pdf/DHSG1/Guide_to_DHS_Statistics_DHS-7_v2.pdf).

## Exposure, geography and confounders

**Use values at band entry**, as requested: MAP and all annual national covariates are matched to the calendar year containing entry. Values are fixed across that band and do not depend on the child's death time. This uses annual estimates, not a measurement made on the entry date. There is no within-band averaging, nearest-year borrowing or interpolation. `calendar_year` is fractional calendar year at entry; `entry_year`/`exposure_year` are the integer source year.

PfPR is population-weighted regional MAP PfPR₂–₁₀, in **percentage points** (0–100), from the existing annual extraction snapshot. It is not weighted by this child's person-time. The child's region is the mother's region at interview; historical residence/migration is not reconstructed. Missing regional exposure is retained with a status flag. `model_ready` means core eligibility plus available region and PfPR; it does **not** mean the adjustment set has been selected or all confounders are observed.

Geographic matching uses exact names, canonical language normalization, curated synonyms and explicit rules/overrides. Fuzzy matching and matching by elimination are not accepted. Donor-based coarse groupings are declared in `config/survey_rules.csv`. Review `region_crosswalk.csv`; unmatched records remain outside the default model input. By default `region` contains country, survey and boundary key: it distinguishes boundary versions. Stable region identifiers across surveys require reviewed `region_id` overrides. A shared longitudinal regional effect must not be assumed from similar names alone.

Older Ethiopian recodes ET41FL–ET81FL use the explicit **+92 CMC** approximation, applied to birth and interview dates. ET8AFL (2024) already uses Gregorian CMCs and is not shifted. Calendar provenance is retained, and disagreement between adjusted interview dates and registry year stops that survey. The approximation does not resolve all Ethiopian/Gregorian month boundaries.

Candidate X variables include sex, multiple birth, birth order, maternal age at birth, preceding birth interval, maternal education, household wealth and urban residence; national HIV prevalence, GDP per capita, health expenditure per capita, political stability and vaccine coverage; and survey-region/national WASH and wasting measures. Household conditions and maternal education were measured at interview, not retrospectively at band entry. Wasting and some other variables may be downstream of malaria. Their inclusion as causal confounders needs a separate decision.

There is **no missing-value imputation or automatic adjustment-set selection**. First births have a structural missing preceding-birth interval; do not require that column without a deliberate first-birth parameterization. Survey covariates are clearly named `_survey_region_pct` and `_survey_national_pct`; national values never silently replace missing regional ones. Immunisation values are retained only when the snapshot identifies a WUENIC estimate; previously assumed pre-series/no-series zeros are treated as missing. The existing HIV snapshot retains its own source assumptions (including midpoint conversion of censored UNAIDS counts); it measures ages 0–14 and does not supply Nigeria or Comoros. Log HIV is missing at zero, with the raw zero retained.

## Inputs and current scope

Paths are declared in `00_config.R`; no network access or package installation occurs.

| Input | Purpose |
|---|---|
| `data/derived_dhs/survey_registry.csv` | Survey inventory and paths to local Births Recodes (`.rds` or `.dta`) |
| `data/derived_dhs/map_pfpr_by_survey_region.csv` | Survey-specific boundary names/keys |
| `data/derived_dhs/map_pfpr_window_years.csv` | Annual regional MAP exposure |
| `data/derived_dhs/hiv_prevalence_country_year.csv` | Annual national HIV prevalence (%) |
| `data/derived_dhs/unicef_immunisation_country_year.csv` | Annual WUENIC coverage and source-status flags |
| `data/wb_gdp_pc.csv`, `data/wb_hexp_pc.csv`, `data/wb_polstab.csv` | Annual World Bank GDP and health expenditure per capita in current US dollars; political stability estimate on the WGI scale |
| `data/derived_dhs/statcompiler_covariates.csv` | Optional survey-level WASH and wasting candidates |
| `R_cbh/config/*.csv` | Explicit survey calendar/geography choices and overrides |

This first implementation reuses **external-data extraction snapshots** already present locally; it rebuilds child histories from the local recodes. Migration of ZIP ingestion, registry discovery, raster extraction and external-panel downloading into this new folder is subsequent work. No old R code is sourced at runtime. Four Lesotho surveys in the registry currently lack MAP geography and are reported as unavailable, not assigned zero exposure.

## Load for a model

```r
source("R_cbh/load_pipeline.R")

# Examples of covariates, not a final causal adjustment set.
d <- cbh_load_analysis(
  "data/derived_cbh",
  required_covariates = c("sex", "maternal_age_birth", "hiv_prev_pct")
)
attr(d, "selection_report")  # aggregate exclusions and selected deaths by survey
attr(d, "skipped_surveys")
```

The reader uses available PfPR/region by default, filters only the explicitly requested covariates, and supplies unordered `age_band`, `survey`, `country`, `region` and `country_age` factors for `mgcv`. Set `surveys` to reduce memory use. `model_ready_only = FALSE` returns eligible audit rows including missing exposure; such rows are not directly fit-ready. Loading all surveys combines several million records and requires substantially more memory than the per-survey builder.

`survey_weight = V005 / 10^6` is retained. `analysis_weight` is normalized to mean 1 within survey **after the reader's row selection**. If you later filter rows, normalize again. These are candidate likelihood weights, not a survey-design variance estimator; mother, PSU and stratum IDs are retained for a future clustering/bootstrap strategy. The stratum variable is recorded, but its survey-specific design interpretation still requires review. Use `na.action = na.fail` when fitting to prevent further silent exclusions.

## Audit files and variable guide

All generated files are under the same ignored output directory:

| File | Content |
|---|---|
| `survey_manifest.csv` | Build/cache/failure status, rows, deaths, exposure availability and signatures |
| `eligibility_flow.csv` | Counts by survey/band: reach alive, five-year window, incomplete bands, year coverage, eligible deaths |
| `child_checks.csv` | Mutually exclusive validity/exclusion reasons before band expansion |
| `region_crosswalk.csv` | Geographic labels, mappings, matching methods and versioned region identifiers |
| `variable_missingness.csv` | Missing counts by survey and output column |
| `input_checks.csv`, `source_files.csv` | Availability and provenance |

| Variables | Meaning |
|---|---|
| `child_id`, `mother_id`, `psu`, `stratum` | Pseudonymous within-snapshot identifiers; survey-qualified |
| `country`, `survey`, `survey_year`, `region`, `regkey` | Country ISO3, survey key and geographic keys |
| `age_band`, `age_band_index`, `age_lo`, `age_hi` | Unordered model group; numeric order; half-open completed-month bounds |
| `death`, `band_years`, `log_band_years` | Bernoulli outcome and predetermined interval offset |
| `birth_cmc`, `interview_cmc`, `band_entry_cmc`, `band_end_cmc` | Gregorian/adjusted century-month codes |
| `entry_year`, `exposure_year`, `calendar_year` | Integer annual join year and fractional entry time |
| `pfpr_pct`, `pfpr_status`, `model_ready` | Exposure in percentage points, join status and core inclusion flag |
| `sex`, `urban`, `wealth_quintile`, `multiple_birth`, `first_birth` | Male/female; urban 1/rural 0; wealth 1–5; binary birth indicators |
| `maternal_age_birth`, `maternal_education_years`, `birth_order`, `preceding_birth_interval_months` | Candidate individual covariates in named units |
| `hiv_prev_pct`, `hib3_pct`, `pcv3_pct`, `rotavirus_pct` | National entry-year percentages, each with a status column |
| `gdp_pc`, `health_expenditure_pc`, `political_stability` | National entry-year current US dollars per capita for the first two; WGI estimate scale for political stability; each with status |
| `log_hiv_prev`, `log_gdp_pc`, `log_health_expenditure_pc` | Natural logarithms for positive observations only |
| `imp_water_*_pct`, `imp_sanit_*_pct`, `wasting_*_pct` | Survey-time candidates, separated into region/national measurements with status |
| `birth_date_flag`, `death_age_flag`, `death_age_source`, `death_age_unit`, `neonatal_death_days`, `death_band_b6_b7_disagree` | DHS date/age precision and reconciliation flags; reported day values are retained for all day-reported deaths |
| `cmc_offset_months`, `calendar_conversion`, `stratum_variable`, `months_before_interview_at_entry`, `survey_covariate_year` | Calendar/design/measurement provenance |

The synthetic tests cover death-age boundaries, complete-band censoring, exact five-year eligibility, year alignment, calendar conversion, invalid histories, deterministic geography, missing/duplicate joins, cache invalidation and post-selection weights. Aggregate verification against local surveys is documented in [BUILD_VALIDATION.md](BUILD_VALIDATION.md).

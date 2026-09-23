# Imputed-covariate sensitivity on the DHS and MICS sample (v6)

This folder refits the 17-covariate primary specification after imputing every remaining covariate gap, so that all
MAP-eligible DHS and MICS records are retained. The version is `primary_map_regional17_dhsmics_imputed_gamma2_v6`.
It is compared with the complete-case DHS+MICS primary, `primary_map_regional17_dhsmics_gamma2_v5`.

The methods are those of the DHS-only version, `primary_map_regional17_imputed_gamma2_v4`, which runs through
`R_cbh/primary/run_regional.R --imputed`. This version lives in its own folder so that v5's code provenance is not
touched: that provenance covers `R_cbh/primary/settings.R` and `00_prepare_regional.R`.

Run from the project root, in order. The fits and the multiple-imputation refits are long, so run them detached.

1. `01_impute_national.R` fills the national series on the 40-country panel:
   - the 2001 WGI round, by interpolation;
   - South Sudan's political stability before independence and its GDP outside 2008–2015, from GAMs;
   - health expenditure for Zimbabwe 2000–2009, Somalia 2000–2012, South Sudan outside 2017–2023 and every country in
     2024, from the v4 GAM.

   Draws are saved for multiple-imputation propagation. Outputs go to
   `data/derived_cbh/regional_adjustment/imputed_dhsmics_v6/` and `results/cbh/covariate_imputation_dhsmics_v6/`.
2. `02_impute_regional.R` imputes the combined DHS and MICS survey-region overlay jointly with `mice`: predictive
   mean matching, m = 10, with the v4 specification.
3. `03_prepare.R` builds the analysis dataset from both survey manifests. It uses:
   - the imputed overlay and national panel, joined by country and band entry year;
   - the extended HIV panel.

   It stops unless every MAP-eligible record is retained: 8,797,963 records, 123,419 deaths, 1,457 survey-regions,
   166 surveys and 40 countries. The DHS part must equal v4.
4. `04_fit.R` fits the seven age-band models with the primary formula, the reference knots and gamma = 2.
5. `05_propagate.R` runs the multiple-imputation check. It refits in each of the 10 imputed datasets, with the
   smoothing parameters fixed at the point fit, and pools with Rubin's rules.
6. `06_report.R` writes the comparison with v5, the supplementary manuscript figure and caption, and `REPORT.md`.

The South Sudan pre-independence fills are a modelling choice. The alternatives are whole-Sudan values, carrying the
2011 values back, or omitting MC_SSD2010, whose 57,509 records are 0.65% of the sample. The imputed values are listed
in `results/cbh/covariate_imputation_dhsmics_v6/somalia_south_sudan_zimbabwe_series.csv`.

# bam (fREML, discretised) versus gam (REML)

Added 24 September 2026. The exposure and every covariate are constant within a survey-region × entry-month cell and band width is constant within a band, so the child-level Bernoulli likelihood equals a binomial likelihood on cells (60,000–74,000 cells per band). `gam(method = "REML")` and `bam(method = "fREML", discrete = FALSE)` were fitted to the cells with the same formula, reference knots and gamma = 2, and compared with the saved child-level `bam(method = "fREML", discrete = TRUE)` fits.

**Model compared:** the primary without the survey-region random intercept (`no_region_re_dhsmics_map_gamma2_v1`, about 230 coefficients per band), which gives the same hazard ratios as v7 to within 0.023. The full v7 model (about 1,450 coefficients, 1,226 of them region levels) was started with `gam` on the 48–59 month band but stopped after 65 minutes without finishing (single-threaded reference BLAS); no result was written.

**Result:** across the seven bands, `gam` and `bam` differ by at most 0.0004 in the 40%→20% hazard ratio and 0.0005 in the 20%→0% hazard ratio; covariate log hazard ratios by at most 0.0009; PfPR effective degrees of freedom by at most 0.01. Turning discretisation off changes the hazard ratios by at most 0.0003. `gam`'s REML reaches a slightly higher log-likelihood (0.3–1.0 units on 24,600–172,600) and takes 185–458 seconds per band on cells, against 14–17 seconds for `bam` on cells. For the 1–5 month band `gam` reached directly the curved PfPR solution that `bam` reached only after its strict restart.

Per-band values: [comparison_no_region.csv](comparison_no_region.csv). Reproduce: `Rscript R_cbh/sensitivity/bam_vs_gam/01_compare.R --no-region`.

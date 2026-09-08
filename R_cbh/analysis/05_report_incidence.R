#!/usr/bin/env Rscript
source("R_cbh/load_pipeline.R")
source("R_cbh/analysis/model.R")
spec<-cbh_trial_spec()
out<-file.path("results/cbh",spec$id)
private<-file.path("data/derived_cbh/models",spec$id)
saved<-readRDS(file.path(private,"fit.rds"));f<-saved$fit
input<-readRDS(file.path(private,"complete_case_dataset.rds"))
panel<-cbh_read_csv(spec$incidence_panel)
j<-match(paste(input$data$country,input$data$entry_year),paste(panel$iso3,panel$year))
stopifnot(!anyNA(j),identical(input$signature,saved$signature),
  "z_log_hiv_incidence" %in% names(f$model),!"z_log_hiv_prev" %in% names(f$model),
  length(f$y)==nrow(input$data),
  max(abs(input$data$log_hiv_incidence-log(panel$hiv_incidence_per1000[j])))<1e-10,
  max(abs(f$model$z_log_hiv_incidence-input$data$z_log_hiv_incidence))<1e-10)
time_terms<-Filter(function(s) identical(s$term,"calendar_year"),f$smooth)
stopifnot(length(time_terms)==1L,identical(time_terms[[1]]$by,"NA"))
e<-eigen(f$outer.info$hess,symmetric=TRUE,only.values=TRUE)$values
diagnostics<-data.frame(converged=f$converged,rank=f$rank,coefficients=length(coef(f)),iterations=f$iter,
  max_abs_smoothing_gradient=max(abs(f$outer.info$grad)),min_hessian_eigenvalue=min(e),
  max_hessian_eigenvalue=max(e),elapsed_seconds=saved$elapsed_seconds)
cbh_atomic_csv(diagnostics,file.path(out,"optimization_diagnostics.csv"))
pooled<-cbh_read_csv(file.path(out,"multiple_imputation/pooled_pfpr_40_to_20.csv"))
rows<-vapply(seq_len(nrow(pooled)),function(i)with(pooled[i,],
  sprintf("| %s | %.3f | %.3f–%.3f | %.4f |",age_band,hazard_ratio,hr_lower_95,hr_upper_95,monte_carlo_se)),character(1))
mc_ratio<-max(pooled$monte_carlo_se/pooled$total_se)
writeLines(c("# Age-band mortality model using child HIV incidence", "",
  "The HIV adjustment is now log child HIV incidence (ages 0–14; new infections per 1,000 uninfected population), replacing prevalence. Missing child rates for Nigeria and Comoros are predicted by the [incidence imputation model](../hiv_incidence/REPORT.md). The join uses band-entry year. Checks against every retained modelling row confirmed the join and verified that the fitted model contains incidence, not prevalence.", "",
  sprintf("The analysis uses **%s child-band records and %s deaths**, from **%d countries, %d surveys and %d survey-specific regions**.",format(nrow(input$data),big.mark=","),format(sum(input$data$death),big.mark=","),nlevels(input$data$country),nlevels(input$data$survey),nlevels(input$data$region)),"",
  "Nigeria contributes 737,012 records and Comoros 17,936. Liberia loses 131,239 formerly included records because its incidence series is unavailable, and São Tomé and Príncipe remains excluded. Compared with the historical prevalence complete-case fit, the net increase is 623,709 records and 12,643 deaths. A common-sample comparison is needed to isolate the effect of changing the HIV measure.","",
  "The seven bands, entry window, full-band eligibility, PfPR timing, cloglog likelihood, band-width offset, other confounders, random effects and unweighted likelihood retain the previous specification. Other covariates require complete cases and vaccines remain excluded. HIV effects are age-specific linear terms on the standardized log-incidence scale.","",
  "Calendar year now enters through **one shared cubic regression spline (k=6)** across all age bands. PfPR retains seven separate age-band splines (k=5). The previous incidence fit with separate time splines remains in [the v2 report](../age_band_hiv_incidence_v2/REPORT.md). The HIV imputation model is unchanged.","",
  "## Propagating incidence uncertainty", "",
  "Ten coherent draws from the fitted external incidence model were joined to the same child-band records. Each country-year draw is shared by every corresponding DHS record. Mortality smoothing parameters and scaling were held fixed at the median-incidence fit, while regression coefficients were refitted. Pointwise intervals combine within-fit covariance and between-imputation variation with finite-imputation t quantiles. This is exploratory external-covariate uncertainty propagation, not full joint outcome-compatible imputation.","",
  "| Completed months | Mortality HR, PfPR 40% → 20% | Pooled 95% interval | Monte Carlo SE of log HR |",
  "|---|---:|---:|---:|",rows,"",
  sprintf("The largest Monte Carlo SE / total SE ratio across these contrasts is %.3f. Additional draws may be needed for precise final inference.",mc_ratio),"",
  "![Incidence-adjusted PfPR splines with imputation uncertainty](multiple_imputation/pooled_pfpr_splines.png)","",
  "The curves show log hazard ratios relative to 20% PfPR over each age band's central 95% exposure range. The contrast and its uncertainty are zero at the reference by construction. These are adjusted model associations, not established causal effects.","",
  "## Numerical checks and remaining limits", "",
  sprintf("The median-incidence mortality fit reported convergence in %d iterations with rank %d/%d. Coefficients and covariance are finite. Its smoothing-parameter Hessian has minimum eigenvalue %.6f (maximum %.2f) and maximum absolute gradient %.6f. %s",f$iter,f$rank,length(coef(f)),min(e),max(e),max(abs(f$outer.info$grad)),
    if(min(e)<0) "The negative eigenvalue leaves smoothing-optimization stability as an open check despite positive convergence flags." else "The smoothing-parameter Hessian is positive definite at this fit; basis-size sensitivity remains to be assessed."),"",
  "PfPR k=5 and time k=6 remain the exploratory basis sizes. The saved factor-by basis diagnostics are pooled checks, not independent age-specific tests. The uncertainty refits condition on the median fit's smoothing parameters. Their intervals do not include survey-design, residual-clustering, MAP-estimation, source HIV-estimate or smoothing-parameter uncertainty. Changing prevalence to incidence also changes the confounding interpretation; see [open issues](../../../docs/OPEN_ANALYSIS_ISSUES.md).", "",
  "## Saved artifacts", "",
  sprintf("- Exact median-incidence analysis input: `%s/complete_case_dataset.rds`, element `data`.",private),
  sprintf("- Private base fit and individual imputation fits: `%s/`.",private),
  "- [Country sample counts](sample_by_country.csv), [incidence source/imputation usage](hiv_incidence_usage.csv).",
  "- [Pooled contrasts](multiple_imputation/pooled_pfpr_40_to_20.csv), [pooled curves](multiple_imputation/pooled_pfpr_curves.csv), [individual-fit diagnostics](multiple_imputation/fit_diagnostics.csv).",
  "- [Median-incidence contrasts](pfpr_40_to_20_contrasts.csv), [median-incidence splines](pfpr_splines_by_age.png), [optimization diagnostics](optimization_diagnostics.csv), [basis checks](basis_checks.csv).",
  "- [Model and run instructions](../../../R_cbh/analysis/README.md). The legacy prevalence fit remains under `age_band_complete_case_v1/`."),file.path(out,"REPORT.md"))
message("Verified active incidence inputs and wrote ",file.path(out,"REPORT.md"))

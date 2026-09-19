#!/usr/bin/env Rscript
# Impute the national annual covariate gaps for the imputed-covariate sensitivity
# version (plan section 2.4). Two structural gaps remain after the base build:
#   political stability : no WGI round in 2001 -> linear interpolation 2000/2002
#   health expenditure  : Zimbabwe 2000-2009 and every country's 2024 -> GAM on the
#                         observed log panel with year smooth, log GDP and country effects
# GDP per capita has no gaps and is copied unchanged. Aggregate outputs only.
source("R_cbh/load_pipeline.R")
library(mgcv)
out_private <- "data/derived_cbh/regional_adjustment/imputed_v4"
out <- "results/cbh/covariate_imputation_v4"
for(p in c(out_private,out)) dir.create(p,recursive=TRUE,showWarnings=FALSE)
cfg <- cbh_config()
years <- cfg$first_entry_year:cfg$last_entry_year
manifest <- cbh_read_csv("data/derived_cbh/survey_manifest.csv")
countries <- sort(unique(manifest$country[manifest$status %in% c("built","cached")]))
polstab <- cbh_read_csv(cfg$annual_panels$governance$path)
hexp <- cbh_read_csv(cfg$annual_panels$health_spending$path)
gdp <- cbh_read_csv(cfg$annual_panels$gdp$path)
grid <- expand.grid(iso3=countries,year=years,stringsAsFactors=FALSE)
pick <- function(panel,col) panel[[col]][match(paste(grid$iso3,grid$year),paste(panel$iso3,panel$year))]
grid$political_stability <- pick(polstab,"polstab")
grid$gdp_pc <- pick(gdp,"gdp_pc"); grid$health_expenditure_pc <- pick(hexp,"hexp_pc")
stopifnot(all(is.finite(grid$gdp_pc) & grid$gdp_pc>0))
grid$log_gdp_pc <- log(grid$gdp_pc)
grid$log_health_expenditure_pc <- ifelse(is.finite(grid$health_expenditure_pc) & grid$health_expenditure_pc>0,
  log(grid$health_expenditure_pc),NA_real_)
grid$political_stability_source <- ifelse(is.finite(grid$political_stability),"wgi_observed",NA)
grid$log_health_expenditure_pc_source <- ifelse(is.finite(grid$log_health_expenditure_pc),"ghed_observed",NA)

## ---- political stability: interpolate the missing 2001 round -----------------------
miss_ps <- which(!is.finite(grid$political_stability))
stopifnot(all(grid$year[miss_ps]==2001L))
for(i in miss_ps) {
  a <- grid$political_stability[grid$iso3==grid$iso3[i] & grid$year==2000L]
  b <- grid$political_stability[grid$iso3==grid$iso3[i] & grid$year==2002L]
  stopifnot(is.finite(a),is.finite(b))
  grid$political_stability[i] <- (a+b)/2
  grid$political_stability_source[i] <- "wgi_2001_interpolated"
}

## ---- health expenditure: GAM on the observed log panel -------------------------------
obs <- grid[is.finite(grid$log_health_expenditure_pc),]
obs$iso3 <- factor(obs$iso3,levels=countries); grid$iso3_f <- factor(grid$iso3,levels=countries)
set.seed(20260918)
fit <- gam(log_health_expenditure_pc ~ s(year,k=6)+log_gdp_pc+s(iso3,bs="re")+s(iso3,year,bs="re"),
  data=obs,method="REML")
miss_he <- which(!is.finite(grid$log_health_expenditure_pc))
nd <- data.frame(year=grid$year[miss_he],log_gdp_pc=grid$log_gdp_pc[miss_he],iso3=grid$iso3_f[miss_he])
Xp <- predict(fit,nd,type="lpmatrix")
point <- drop(Xp %*% coef(fit))
n_draws <- 10L
beta_draws <- MASS::mvrnorm(n_draws,coef(fit),fit$Vp)
residual_sd <- sqrt(fit$sig2)
draws <- Xp %*% t(beta_draws) + matrix(rnorm(length(miss_he)*n_draws,0,residual_sd),ncol=n_draws)
grid$log_health_expenditure_pc[miss_he] <- point
last_observed_year <- max(obs$year)
grid$log_health_expenditure_pc_source[miss_he] <- ifelse(grid$year[miss_he]>last_observed_year,
  "gam_one_year_extrapolation","gam_within_series_model")
stopifnot(all(grid$year[miss_he]<=last_observed_year+1L))
grid$health_expenditure_pc <- exp(grid$log_health_expenditure_pc)
grid$iso3_f <- NULL
stopifnot(all(is.finite(grid$political_stability)),all(is.finite(grid$log_health_expenditure_pc)))

## ---- outputs ----------------------------------------------------------------------------
cbh_atomic_csv(grid,file.path(out_private,"national_covariates_imputed.csv"))
cbh_atomic_csv(grid,file.path(out,"national_covariates_imputed.csv"))
cbh_atomic_rds(list(keys=grid[miss_he,c("iso3","year")],log_health_expenditure_pc=draws,
  political_stability_note="deterministic interpolation; no draws",n_draws=n_draws,
  seed=20260918L,model=deparse(formula(fit))),file.path(out_private,"national_imputations.rds"))
# In-sample fit check for the health-expenditure model: leave-one-year-out is not
# needed for a descriptive fill; report residual SD and R^2 on the observed panel.
check <- data.frame(model="log_health_expenditure_pc",n_observed=nrow(obs),n_imputed=length(miss_he),
  n_imputed_2024=sum(grid$year[miss_he]==2024L),n_imputed_zwe_2000_2009=sum(grid$iso3[miss_he]=="ZWE"),
  residual_sd=residual_sd,deviance_explained=summary(fit)$dev.expl,
  political_stability_interpolated=length(miss_ps))
cbh_atomic_csv(check,file.path(out,"national_imputation_summary.csv"))
zwe <- grid[grid$iso3=="ZWE",c("year","health_expenditure_pc","log_health_expenditure_pc_source")]
cbh_atomic_csv(zwe,file.path(out,"zimbabwe_health_expenditure.csv"))
writeLines(c("# National covariate imputation (imputed_v4)","",
  sprintf("Countries: %d built survey countries; years %d-%d.",length(countries),min(years),max(years)),"",
  sprintf("**Political stability (WGI):** %d country-years missing, all in 2001 (no WGI round). Filled by the mean of each country's 2000 and 2002 estimates; deterministic, flagged `wgi_2001_interpolated`.",length(miss_ps)),"",
  sprintf("**Health expenditure per capita (log, current US$):** %d country-years missing: Zimbabwe 2000-2009 (%d) and 2024 for all %d countries. Filled from `%s` fitted by REML on %d observed country-years (residual SD %.3f on the log scale, deviance explained %.1f%%). Point values are fitted means; %d draws from the coefficient posterior plus residual noise are saved for multiple-imputation propagation. Flags: `gam_within_series_model` (Zimbabwe) and `gam_one_year_extrapolation` (2024).",
    length(miss_he),sum(grid$iso3[miss_he]=="ZWE"),length(countries),paste(deparse(formula(fit)),collapse=""),nrow(obs),residual_sd,100*summary(fit)$dev.expl,n_draws),"",
  "GDP per capita has no gaps and is unchanged. Child HIV incidence for Liberia and Sao Tome and Principe is imputed separately by the extended incidence model (`R_cbh/hiv/01_fit_incidence.R --extended`).","",
  "[Imputed panel](national_covariates_imputed.csv) · [Summary](national_imputation_summary.csv) · [Zimbabwe series](zimbabwe_health_expenditure.csv)"),
  file.path(out,"NATIONAL_IMPUTATION.md"))
paths <- c(cfg$annual_panels$governance$path,cfg$annual_panels$health_spending$path,cfg$annual_panels$gdp$path,
  "data/derived_cbh/survey_manifest.csv","R_cbh/covariates/15_impute_national.R")
cbh_atomic_csv(data.frame(file=paths,md5=vapply(paths,cbh_file_hash,"")),file.path(out,"national_imputation_provenance.csv"))
print(check); message("National covariate imputation complete")

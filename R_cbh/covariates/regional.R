# Regional DHS summaries. Source after R_cbh/load_pipeline.R.
# No mortality outcome is used to choose the population for these means.
cbh_regional_binary <- function(x) cbh_decode_category(x,c("no","yes"),c(0,1))

cbh_regional_recode <- function(br, survey, rule, geo) {
  n <- nrow(br)
  num <- function(v) cbh_column(br,v,TRUE)
  column <- function(v) if(v %in% names(br)) br[[v]] else rep(NA_real_,n)
  age <- num("v008")-num("b3")
  recent <- is.finite(age) & age>=0 & age<60
  w <- num("v005")/1e6
  mother <- trimws(cbh_column(br,"caseid"))
  if(all(is.na(mother))) mother <- cbh_key(br,c("v001","v002","v003"))
  mother_ok <- !is.na(mother) & nzchar(mother)
  first_mother <- recent & mother_ok & !duplicated(ifelse(recent & mother_ok,mother,NA_character_))
  sex <- cbh_decode_category(column("b4"),c("male","female"),c(1,2))
  urban <- cbh_decode_category(column("v025"),c("urban","rural"),c(1,2))
  wealth <- cbh_decode_category(column("v190"),c("poorest","poorer","middle","richer","richest"),1:5)
  order <- num("bord"); order[!is.finite(order) | order<1 | order>40] <- NA_real_
  mage <- (num("b3")-num("v011"))/12; mage[!is.finite(mage)|mage<10|mage>55] <- NA_real_
  educ <- num("v133"); educ[!is.finite(educ)|educ<0|educ>30] <- NA_real_
  mult <- num("b0"); mult[!is.finite(mult)|mult<0|mult>8] <- NA_real_
  mult <- as.numeric(mult>0)
  lab <- tolower(cbh_column(br,"b0"))
  mult[lab %in% "single birth"] <- 0
  mult[!is.na(lab) & grepl("multiple|twin|triplet",lab)] <- 1
  values <- list(male_pct=100*(sex==1),multiple_birth_pct=100*mult,
    mean_birth_order=order,mean_maternal_age_birth=mage,
    mean_maternal_education_years=educ,mean_wealth_quintile=wealth,urban_pct=100*(urban==1))
  populations <- c(rep("live_births_last_60_months",4),rep("mothers_with_birth_last_60_months",3))
  names(populations) <- names(values)
  masks <- lapply(names(values),function(v) if(v %in% names(values)[1:4]) recent else first_mother)
  names(masks) <- names(values)

  # DHS IYCF denominator: youngest living child under 24 months residing with
  # their mother, then restrict to age 0-5 months. Never repeat maternal feeding
  # answers across siblings. Use m4 (child breastfeeding status), not v404.
  use_b19 <- "b19" %in% names(br) && any(is.finite(num("b19")))
  infant_age <- if(use_b19) num("b19") else age
  alive <- cbh_decode_category(column("b5"),c("no","yes"),c(0,1))
  home <- cbh_decode_category(column("b9"),c("respondent","lives elsewhere"),c(0,1))
  idx <- num("bidx")
  eligible <- which(mother_ok & alive==1 & home==0 & is.finite(infant_age) &
    infant_age>=0 & infant_age<24 & is.finite(idx))
  eligible <- eligible[order(mother[eligible],idx[eligible])]
  eligible <- eligible[!duplicated(mother[eligible])]
  infant <- rep(FALSE,n); infant[eligible[infant_age[eligible]<6]] <- TRUE
  m4 <- num("m4"); m4lab <- tolower(cbh_column(br,"m4"))
  bf <- rep(NA_real_,n)
  bf[is.finite(m4) & ((m4>=0 & m4<=93)|m4==94)] <- 0
  bf[m4lab %in% c("never breastfed","ever breastfed, not currently breastfeeding","not breastfeeding")] <- 0
  bf[which(m4==95 | m4lab %in% c("still breastfeeding","currently breastfeeding"))] <- 1
  # Recodes contain all-NA placeholders for questionnaire items not administered.
  # Only active standard items are used. Require all four food/liquid domains;
  # missing/don't-know answers in an active item cannot establish exclusivity.
  domains <- list(water="v409",liquids=c("v409a","v410","v410a","v412c",paste0("v413",c("","a","b","c","d"))),
    milk=c("v411","v411a","v412"),solids=c("v412a","v412b",paste0("v414",letters[1:23]),"m39a"))
  decoded <- lapply(unique(unlist(domains)),function(v) cbh_regional_binary(column(v)))
  names(decoded) <- unique(unlist(domains))
  active <- names(decoded)[vapply(decoded,function(x) any(is.finite(x[eligible])),logical(1))]
  module_ok <- all(vapply(domains,function(v) any(v %in% active),logical(1)))
  ebf <- rep(NA_real_,n)
  if(module_ok) {
    a <- do.call(cbind,decoded[active])
    other_yes <- rowSums(a==1,na.rm=TRUE)>0
    food_complete <- rowSums(is.na(a))==0
    ebf[which(bf==0 | other_yes)] <- 0
    ebf[which(bf==1 & !other_yes & food_complete)] <- 1
  }
  values$exclusive_breastfeeding_pct <- 100*ebf
  masks$exclusive_breastfeeding_pct <- infant
  populations <- c(populations,exclusive_breastfeeding_pct="youngest_resident_child_age_0_5_months")
  # For coverage, any nonempty valid denominator is estimable; flag <25 / <50
  # separately instead of introducing an undeclared sample-size exclusion.
  regions <- unique(na.omit(geo$regkey))
  rows <- list()
  for(v in names(values)) for(r in regions) {
    target <- which(geo$regkey==r & masks[[v]] & is.finite(w) & w>0)
    observed <- target[is.finite(values[[v]][target])]
    val <- if(length(observed)) weighted.mean(values[[v]][observed],w[observed]) else NA_real_
    rows[[length(rows)+1L]] <- data.frame(survey=survey$svkey,country=survey$iso3,regkey=r,
      variable=v,value=val,eligible_n=length(target),observed_n=length(observed),
      missing_n=length(target)-length(observed),weighted_n=sum(w[observed]),
      population=populations[[v]],source="local_BR_weighted_summary",
      small_denominator=length(observed)<25,low_precision=length(observed)<50)
  }
  list(values=cbh_bind(rows),feeding=data.frame(survey=survey$svkey,
    age_source=if(use_b19) "b19" else "v008_minus_b3",module_available=module_ok,
    active_items=paste(active,collapse=";"),eligible_infants=sum(infant),
    classified_infants=sum(infant & is.finite(ebf)),
    national_ebf_pct=if(any(infant & is.finite(ebf) & is.finite(w) & w>0))
      weighted.mean(100*ebf[infant & is.finite(ebf) & is.finite(w) & w>0],w[infant & is.finite(ebf) & is.finite(w) & w>0]) else NA_real_))
}

cbh_regional_spec <- function() {
  regional <- c("male_pct","multiple_birth_pct","mean_birth_order","mean_maternal_age_birth",
    "mean_maternal_education_years","mean_wealth_quintile","urban_pct","dtp3_pct","measles_pct",
    "facility_delivery_pct","exclusive_breastfeeding_pct","short_birth_interval_pct",
    "improved_water_pct","improved_sanitation_pct","electricity_pct")
  annual <- c("log_hiv_incidence","log_gdp_pc","log_health_expenditure_pc","political_stability",
    "hib3_pct","pcv3_pct","rotavirus_pct")
  list(id="regional_adjustment_v1",regional=regional,annual=annual,covariates=c(regional,annual))
}

cbh_regional_formula <- function() {
  as.formula(paste("death ~ s(pfpr_pct,bs='cr',k=5) + s(calendar_year,bs='cr',k=6) +",
    paste(paste0("z_",cbh_regional_spec()$covariates),collapse=" + "),
    "+ s(survey,bs='re') + s(country,bs='re') + s(region,bs='re') + offset(log(band_years))"))
}

# One row per survey-region, before repeating covariates over child-band rows.
# Each variable uses its own observed donor regions from the same survey.
cbh_regional_mean_fill <- function(wide,variables=cbh_regional_spec()$regional) {
  wide <- as.data.frame(wide)
  cbh_require(wide,c("survey","regkey",variables),"Regional mean substitution")
  cbh_unique(wide,c("survey","regkey"),"Regional mean substitution")
  for(v in variables) {
    original <- wide[[v]]
    missing <- !is.finite(original)
    donor_n <- integer(nrow(wide)); donor_mean <- rep(NA_real_,nrow(wide))
    for(s in unique(wide$survey)) {
      ix <- which(wide$survey==s)
      donors <- ix[is.finite(original[ix])]
      donor_n[ix] <- length(donors)
      if(length(donors))donor_mean[ix] <- mean(original[donors])
    }
    fill <- missing & is.finite(donor_mean)
    wide[[v]][fill] <- donor_mean[fill]
    wide[[paste0(v,"_before_regional_mean")]] <- original
    wide[[paste0(v,"_regional_mean_imputed")]] <- fill
    wide[[paste0(v,"_donor_regions")]] <- donor_n
    wide[[paste0(v,"_regional_mean_source")]] <- ifelse(fill,"mean_available_regions_same_survey",
      ifelse(missing,"no_available_region_same_survey","retained_input_estimate"))
    stopifnot(identical(wide[[v]][!missing],original[!missing]),all(!fill | donor_n>0L))
  }
  wide
}

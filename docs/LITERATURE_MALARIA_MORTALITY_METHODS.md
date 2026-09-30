# How IHME (GBD) and UN IGME (CA CODE) estimate under-5 malaria deaths

Compiled 29 September 2026 from the primary sources below. Each summary was checked by a second reviewer who opened every citation. Statements marked *inferred* are our reading of the method, not text from the sources.

## IHME / Global Burden of Disease

**Main methods papers.**
- Gething PW, Casey DC, Weiss DJ, et al. Mapping *Plasmodium falciparum* mortality in Africa between 1990 and 2015. *N Engl J Med* 2016; 375: 2435–45. doi:10.1056/NEJMoa1606701. This paper joined the Malaria Atlas Project (MAP) to GBD 2015. It set out the incidence × case-fatality approach for Africa.
- Weiss DJ, Lucas TCD, Nguyen M, et al. Mapping the global prevalence, incidence, and mortality of *P falciparum*, 2000–17. *Lancet* 2019; 394: 322–31. doi:10.1016/S0140-6736(19)31097-9. This paper added a geospatial case-fatality model and extended the method worldwide. Its results "form the malaria estimates for the Global Burden of Disease 2017 study".

**Most recent.**
- GBD 2023 Causes of Death Collaborators. Global burden of 292 causes of death in 204 countries and territories and 660 subnational locations, 1990–2023. *Lancet* 2025; 406: 1811–72. doi:10.1016/S0140-6736(25)01917-8. The malaria methods are in appendix 1, "Input data and methodological summary for malaria", PDF pp. 574–79.
- Companion MAP input paper: Weiss DJ, Dzianach PA, Saddler A, et al. Mapping the global prevalence, incidence, and mortality of *P falciparum* and *P vivax* malaria, 2000–22. *Lancet* 2025; 405: 979–90. doi:10.1016/S0140-6736(25)00038-8.
- The GBD 2021 causes-of-death appendix (*Lancet* 2024; 403: 2100–32) lists a malaria write-up but omits it from the PDF.

**Model (P. falciparum in sub-Saharan Africa, GBD 2023).**
1. MAP maps PfPR2–10 geostatistically from household surveys (DHS, MIS and others) and published prevalence surveys. It converts prevalence to clinical incidence with an ensemble of three microsimulation models (EMOD, OpenMalaria, malariasimulation).
2. Effective-treatment surfaces combine treatment seeking, the mix of drug classes and drug efficacy. Incidence without effective treatment gives *untreated incidence*.
3. An untreated case-fatality rate (uCFR) is estimated from geo-referenced cause-of-death data, mostly verbal autopsy (VA). For each site-year, malaria deaths are the VA malaria cause fraction × national all-cause mortality. Dividing by MAP untreated incidence gives the site-year uCFR.
   - A geostatistical model is fitted to logit all-ages uCFR.
   - Its covariates are travel time to cities, the proportions of adults and of infants, log all-cause mortality and sickle-cell trait prevalence.
4. Pixel deaths = uCFR × untreated incidence × population. These are aggregated to countries. This is our paraphrase of appendix 1, PDF p. 576: "Pixel-year predictions of CFR were then multiplied by the untreated incidence rate rasters from the MAP cube to yield pixel-year mortality rate estimates, which were then multiplied by pixel-year population to derive pixel-year malaria death counts." Gething 2016 (Methods, "Malaria case fatality rate") gives the same step without population: "applying the estimated case fatality rate to the estimate of untreated incidence for each 5-km2 grid cell". Because the GBD 2023 uCFR is all-ages, this step gives all-age deaths; under-5 deaths come from step 5.
5. Deaths are split by age and sex using separate CODEm "fatal age pattern" models (covariates: incidence and effective treatment). *Inferred:* the under-5 share is therefore set by VA/VR age patterns, not by the case-fatality model.
6. Cause fractions above an MR-BRT ceiling fitted to the maximum observed African cause fractions are capped. This step is new in GBD 2023.
7. CoDCorrect rescales all causes to the GBD all-cause envelope (+5.4% for global malaria in 2023).

GBD 2023 estimates 670,000 malaria deaths in 2023 (95% UI 261,000–1,258,000).

**Primary inputs.**
- The GBD all-cause envelope, from VR, sample VR, birth histories and censuses. It is used three times: to convert cause fractions to rates, as a uCFR covariate, and in CoDCorrect.
- The GBD cause-of-death database. For malaria: 20,168 VR, 953 sample VR and 1,800 VA site-years. Only 280 of MAP's 4,750 location-years are African, mostly collected before 2010; the most recent is from 2019.
- MAP PfPR, incidence and effective-treatment surfaces, including ITN, IRS and environmental covariates.
- Population rasters.
- WHO PULSE surveys, used for a COVID-19 treatment adjustment in 26 African countries for 2020–22.

GBD 2023 notes that "systematic literature reviews for malaria were not conducted".

## UN IGME / CA CODE

UN IGME estimates all-cause mortality only. The under-5 malaria series on childmortality.org comes from the Child and Adolescent Causes of Death Estimation (CA CODE) project, formerly MCEE. It is applied to the UN IGME envelope.

**Main methods paper.**
- Perin J, Mulick A, Yeung D, et al. Global, regional, and national causes of under-5 mortality in 2000–19. *Lancet Child Adolesc Health* 2022; 6: 106–15. doi:10.1016/S2352-4642(21)00311-4.
- The statistical model is in Mulick AR, Oza S, Prieto-Merino D, et al. *J R Stat Soc A* 2022; 185: 2097–120. doi:10.1111/rssa.12853.
- PfPR replaced the earlier "malaria index" covariate in Liu L, Oza S, Hogan D, et al. *Lancet* 2016; 388: 3027–35. doi:10.1016/S0140-6736(16)31593-8.

**Most recent (the "CA-CODE 2026" series used here).**
- Perin J, Prieto-Merino D, Wahi A, et al. Systematic estimates of global causes of neonatal and under 5 mortality in 2000–24: secondary data analysis using bayesian multinomial logistic regression. *BMJ* 2026; 393: e088686. doi:10.1136/bmj-2025-088686. Released on childmortality.org on 18 March 2026. Code: github.com/JHU-CACODE/BMJ2026-MortUnder5.

**All-cause envelope.**
- Alkema L, New JR. Global estimation of child mortality using a Bayesian B-spline bias-reduction model. *Ann Appl Stat* 2014; 8: 2122–49. doi:10.1214/14-AOAS768.
- Most recent: Sharrow D, Hug L, Liu Y, Fell GW, You D. *BMJ* 2026; 393: e088684. doi:10.1136/bmj-2025-088684.
- UN IGME, *Levels & Trends in Child Mortality: Report 2025* (UNICEF, 2026).

**Model.**
1. **Envelope.** A B3 penalised spline is fitted through VR, censuses and household-survey birth histories (DHS, MICS, MIS), with source-specific bias terms and no covariates. HIV and crisis deaths are added back. Rates are extrapolated to 2024 for most countries; in the 2025 round, 162 of 200 countries had no 2024 data. Deaths come from WPP 2024 births.
2. **Cause fractions for ages 1–59 months.** Malaria is not a neonatal cause in CA CODE. High-mortality countries (U5MR above 35) use a Bayesian multinomial logistic regression with LASSO priors, fitted to VA studies.
   - Input: 411 data points, 241,048 deaths, 42 countries in total.
   - Since 2026, VA-reported causes are first calibrated for misclassification with CHAMPS MITS matrices.
   - Covariates: U5MR, **PfPR** (the only malaria-specific one), GNI, measles, Hib3, PCV and rotavirus coverage, sanitation, wasting and year.
   - Countries with U5MR between 25 and 35 average the VA and VR models; low-mortality countries use a VR model.
   - The response is the cause mix, not the death rate. Each covariate is a national country-year value, attached to a study at its mid-year and standardised with the study-data mean and SD (U5MR: mean 76.5, SD about 45 per 1,000; PfPR: mean 0.07, SD 0.14, on a 0–1 scale). Coefficients are log-odds of each cause relative to pneumonia, per SD.
   - Posterior means for malaria, from the repository's run log `Postneonates/code/Estimate_VA_stan.Rout`: U5MR +0.08, PfPR +0.52, GNI −0.48, year +0.24, sanitation −0.20, Hib3 +0.13. U5MR mainly moves other causes (injuries −0.22, congenital −0.27, diarrhoea +0.16).
   - For prediction, U5MR comes from `vavr_covariates_20250930.dta`, whose source is not documented (presumably UN IGME).
3. **Low-burden countries.** Malaria is replaced by WHO Global Malaria Programme estimates. This covers about 13 African countries, including Ethiopia, Madagascar, Zimbabwe, Eritrea and Somalia.
4. **Deaths.** The cause fractions are multiplied by the UN IGME 1–59-month envelope, after HIV (from UNAIDS), measles (from WHO) and crisis deaths are removed.

**Primary inputs.**
- The UN IGME envelope.
- VA studies (national surveys, HDSS and the literature) calibrated with CHAMPS.
- VR data (WHO Mortality Database).
- Country-year covariates, including national PfPR. The source of the PfPR covariate is not stated in the main text or code (probably MAP; unverified).
- WHO malaria estimates for low-burden countries.
- Separately estimated HIV, measles, TB and crisis deaths.

## What this means for comparisons with the PfPR-ACM model

- Both estimators rest on VA for the malaria share of deaths in high-burden Africa. IHME uses VA to calibrate a case-fatality rate applied to MAP incidence. CA CODE uses VA to fit cause fractions, with PfPR as one covariate.
- Neither estimator captures indirect malaria mortality.
- Their levels are set by different all-cause envelopes: GBD for IHME, UN IGME for CA CODE. This is the main reason for country gaps such as DRC, where UN IGME has 1.85 times IHME's all-cause under-5 deaths in 2024.
- WHO's World Malaria Report 2025 reuses the 2021 CA CODE fractions for 2022–24, a different vintage from the 2026 series.

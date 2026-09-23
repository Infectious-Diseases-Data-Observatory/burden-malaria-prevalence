#!/usr/bin/env python3
"""Write results/mics_inventory/ELIGIBILITY.md from eligibility.csv (aggregate only)."""
import csv
rows = list(csv.DictReader(open("results/mics_inventory/eligibility.csv")))
lab = lambda s: s.replace("_v01_M", "")
def section(v):
    rs = sorted((r for r in rows if r["verdict"] == v), key=lambda r: r["survey"])
    t = ["| Survey | Year | Regions | Under-5 deaths, last 5 y | Boundary | Note |", "|---|---:|---:|---:|---|---|"]
    t += [f"| {lab(r['survey'])} | {r['interview_year']} | {r['regions']} | {int(r['u5_deaths_last5y']):,} | "
          f"{r['boundary'][:95]} | {r['reason'][:170]} |" for r in rs]
    return t, len(rs), sum(int(r["u5_deaths_last5y"]) for r in rs), len({r["survey"][:3] for r in rs})
inc, ni, di, ci = section("include"); wk, nw, dw, cw = section("include_with_work"); ex, ne, de, ce = section("exclude")
md = ["# Can MICS surveys enter the primary analysis?", "",
 "Assessment of the 93 MICS surveys in `data/MICS_Datasets/` against the requirements of the primary child age-band analysis. "
 "Code: `R_mics/01_inventory.R` (inventory), `R_mics/02_inventory_tables.R` (tables), `R_mics/03_eligibility_report.py` (this page). "
 "Details: [inventory](INVENTORY.md), [eligibility table](eligibility.csv).", "",
 "Verdicts come from two independent routes: a scripted inventory with direct checks, and per-country audit agents each followed by "
 "an adversarial re-check. The two agreed on 38 of the 40 surveys both covered; the two corrections are applied below. "
 "Malawi, Mozambique and Zimbabwe were checked directly only.", "",
 "## Summary", "",
 "- **46 of 93 surveys have no complete birth history** (all 24 MICS2, 12 of 16 MICS3, 10 of 19 MICS4). They carry only summary "
 "questions (children ever born, children dead, date of last birth) and cannot support the age-band design.",
 "- **47 surveys in 28 countries have a complete birth history**, structurally equivalent to the DHS Births Recode. Seven use "
 "non-standard names (Portuguese in Mozambique 2008, lowercase in Zimbabwe 2009, month-year dates or unit-coded ages elsewhere) but are complete.",
 f"- **{ni} surveys in {ci} countries can enter with the existing method** ({di:,} under-five deaths among births in the 5 years before "
 f"interview), **{nw} in {cw} countries after extra work** ({dw:,} deaths), and **{ne} cannot** ({de:,} deaths). The eligible surveys would "
 "add roughly a third to the 75,726 deaths in the current primary sample, including two countries with no DHS survey "
 "(Central African Republic, Guinea-Bissau).", "",
 "## What each requirement needs for MICS", "",
 "- **Birth history.** Complete in all 47. Variable names differ across rounds and must be mapped (BH4C or CCDOB or month and year for "
 "birth date; BH9C, or units and number, for age at death; WDOI, WM6C or month and year for interview).",
 "- **Regional covariates.** For DHS, 4 of the 13 came from the recodes and 9 from published DHS survey-region indicator tables. MICS has "
 "no such tables, so all 13 must be computed from the microdata. Every source exists in the eligible surveys: water, toilet type and "
 "electricity in `hh.sav`; DTP3, measles and anthropometry in `ch.sav` or `who_z.sav`; place of delivery and education in `wm.sav`; age "
 "at first birth and birth intervals from `bh.sav`.",
 "- **Definition differences to resolve.** MICS records education as level and grade, not completed years, so a country-specific mapping "
 "is needed. It asks about place of delivery only for the last birth in the 2 years before interview, against the DHS 5-year window. "
 "MICS3 to MICS5 supply both old NCHS-referenced and WHO z-scores; only the WHO versions (`HAZ2`, `WHZ2`) match. MICS6 electricity has "
 "an off-grid category whose coding must match the DHS definition. Guinea 2016, Comoros 2022 and Chad 2019 lack usable recall doses, "
 "so DTP3 and measles are card-based there.",
 "- **Regions and MAP prevalence.** Most MICS regions match an existing DHS boundary file in the same country name for name, "
 "sometimes after merging MICS regions into a DHS grouping (Chad, Mauritania, Togo) or aggregating districts to regions (Malawi 2013). "
 "Guinea-Bissau and the Central African Republic have no DHS boundary, so new admin-1 polygons are needed. Nigeria's DHS surveys in the "
 "study use the six zones, which the MICS zone variable matches.",
 "- **National series.** Complete for every eligible country. They block Somalia (health expenditure absent before 2013), South Sudan "
 "(political stability absent before 2011, health expenditure before 2017) and Sao Tome (no child HIV incidence series), as they block "
 "the equivalent DHS records.", "",
 f"## Can enter with the existing method ({ni} surveys)", ""] + inc + ["",
 f"## Can enter after extra work ({nw} surveys)", ""] + wk + ["",
 f"## Cannot enter the complete-case primary ({ne} surveys)", ""] + ex + ["",
 "## Caveats", "",
 "- Adding MICS changes the primary sample definition from DHS and MIS to DHS, MIS and MICS, so it would be a new primary version rather "
 "than a sensitivity. The survey random effect absorbs level differences between programmes, but not differences in covariate definitions.",
 "- Sub-national MICS (five Kenyan county, province or settlement samples, Dakar city, southern Madagascar, two Somali zones) are excluded "
 "for consistency with the national-survey design. Their wealth quintiles are also ranked within the sub-national sample, not nationally, "
 "so that covariate would have to be rebuilt against a national reference before any could enter. Nyanza 2011 is the strongest candidate.",
 "- Several MICS surveys overlap DHS surveys in the same country and years (for example Gambia 2018 with DHS 2019, Sierra Leone 2017 with "
 "DHS 2019). They are independent samples, not duplicates, so both can enter.",
 "- Nigeria 2021 (no anthropometry) and Sao Tome 2014 and 2019 (no child HIV series) fit the methods the imputed-covariate sensitivity "
 "already uses, so they could enter that version. Somalia and South Sudan would need several years of national series extrapolated, "
 "which that version does not currently do.",
 "- The under-five death counts here are among births in the 5 years before interview. The primary analysis counts deaths in eligible "
 "complete age bands, which is a similar but not identical quantity.", ""]
open("results/mics_inventory/ELIGIBILITY.md", "w").write("\n".join(md) + "\n")
print(f"ELIGIBILITY.md: {ni} include, {nw} with work, {ne} exclude")

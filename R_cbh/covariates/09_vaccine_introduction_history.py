"""Summarise a cached public WHO introduction history; no imputation or fits.

Source: https://xmart-api-public.who.int/WIISE/MT_AD_INTRO_LONG
Query: (WHO_REGION eq 'AFRO' or ISO_3_CODE eq 'SOM' or ISO_3_CODE eq 'DJI'
or ISO_3_CODE eq 'SDN') and (ANTIGEN eq 'HIB' or ANTIGEN eq 'PNEUMO_CONJ'
or ANTIGEN eq 'ROTAVIRUS'). Select ISO_3_CODE,COUNTRYNAME,WHO_REGION,YEAR,
ANTIGEN,INTRO,NATIONWIDE,PARTIALLY,HIGHRISKGROUP; order by ISO_3_CODE,ANTIGEN,YEAR;
$count=true. Retrieved 2026-09-17. Public API identified from WHO portal script.
"""
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

src = Path("data/derived_cbh/regional_adjustment/who_introduction/africa_history_2026-09-17.json")
out = Path("results/cbh/regional_adjustment_unicef_v2")
obj = json.loads(src.read_text())
raw = obj["value"]
assert len(raw) == obj["@odata.count"] and "@odata.nextLink" not in obj
assert len({(r["ISO_3_CODE"], r["ANTIGEN"], r["YEAR"]) for r in raw}) == len(raw)
assert {r["INTRO"] for r in raw} <= {"Yes", "Yes (P)", "No"}
study = {r["country"] for r in csv.DictReader((out / "selection_by_country.csv").open())
         if r["baseline"] == "eligible_MAP"}
groups = defaultdict(list)
for r in raw:
    # 48 sub-Saharan African sovereign states; Sudan/Algeria retained in raw only.
    if r["ISO_3_CODE"] not in {"SDN", "DZA"}:
        groups[r["ISO_3_CODE"], r["ANTIGEN"]].append(r)
assert len(groups) == 48 * 3
assert study <= {c for c, a in groups}

def write_csv(name, rows):
    with (out / name).open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

rows = []
for (country, antigen), history in sorted(groups.items()):
    history.sort(key=lambda r: r["YEAR"])
    introduced = [r for r in history if r["INTRO"].startswith("Yes")]
    nationwide = [r for r in history if r["INTRO"] == "Yes"]
    discrepancy = country == "NGA" and antigen == "PNEUMO_CONJ"
    rows.append(dict(country=country, country_name=history[0]["COUNTRYNAME"],
        antigen=antigen, in_MAP_eligible_sample=country in study,
        first_reported_introduction_year=introduced[0]["YEAR"] if introduced else "",
        initial_status=introduced[0]["INTRO"] if introduced else "",
        first_reported_nationwide_year=nationwide[0]["YEAR"] if nationwide else "",
        last_reported_year=history[-1]["YEAR"], last_reported_status=history[-1]["INTRO"],
        source_url="https://immunizationdata.who.int/global/wiise-detail-page/vaccine-introduction-in-country_name?ISO_3_CODE="+country,
        documented_launch_date="2014-12-22" if discrepancy else "",
        launch_source="https://www.afro.who.int/countries/nigeria/news/nigeria-introduces-new-vaccine-pcv-10" if discrepancy else "",
        note="WHO launch report predates annual introduction status; do not classify 2014 as wholly pre-introduction." if discrepancy else ""))
write_csv("vaccine_introduction_history.csv", rows)
write_csv("vaccine_introduction_source_manifest.csv", [dict(source=str(src),
    sha256=hashlib.sha256(src.read_bytes()).hexdigest(), retrieved="2026-09-17",
    endpoint="https://xmart-api-public.who.int/WIISE/MT_AD_INTRO_LONG", source_rows=len(raw),
    summary_countries=48, study_countries=len(study))])

lookup = {(r["country"],r["antigen"]):r for r in rows}
lines = ["# Hib, PCV and rotavirus introduction histories", "",
    "WHO annual introduction-status records, retrieved 17 September 2026. "
    "This is a descriptive source audit; no missing coverage values or model inputs were changed.", "",
    "A single year means the first reported nationwide introduction. `P year → year` "
    "means first reported partial introduction followed by first reported nationwide introduction. "
    "These are annual reporting dates, not necessarily exact launch dates, and nationwide introduction "
    "does not mean 100% coverage. `No through 2025` means no introduction recorded in this snapshot.", "",
    "WHO cautions that introduction can precede the first year of consistent reporting. "
    "Nigeria PCV is a verified example: the [WHO launch announcement]"
    "(https://www.afro.who.int/countries/nigeria/news/nigeria-introduces-new-vaccine-pcv-10) "
    "dates launch to 22 December 2014, whereas the annual table first records partial introduction "
    "in 2015 and nationwide introduction in 2017. The raw annual values are preserved below. "
    "No other dates have been individually reconciled with launch announcements.", "",
    "Source: [WHO introduction portal]"
    "(https://immunizationdata.who.int/global/wiise-detail-page/vaccine-introduction-in-country_name). "
    "CSV contains each country's direct source link and separate partial/nationwide fields. "
    "The 48-country scope excludes Algeria and Sudan; it includes Djibouti and Somalia. "
    "The study column marks the 36 countries with MAP-eligible records before expanded covariate selection.", "",
    "| Country | Study | Hib | PCV | Rotavirus |", "|---|:---:|---|---|---|"]
for country in sorted({r["country"] for r in rows}, key=lambda c: lookup[c,"HIB"]["country_name"]):
    cells = []
    for antigen in ("HIB","PNEUMO_CONJ","ROTAVIRUS"):
        r = lookup[country,antigen]
        first, full = r["first_reported_introduction_year"], r["first_reported_nationwide_year"]
        cells.append(f"P {first} → {full}" if first and first != full else str(first) if first else "No through 2025")
    lines.append("| " + " | ".join([lookup[country,"HIB"]["country_name"],"Yes" if country in study else "—"]+cells) + " |")
lines += ["", "For analysis, documented pre-introduction years are candidates for structural zero "
    "routine-programme coverage. Partial-rollout and launch years must not automatically receive zero; "
    "missing values during or after introduction remain coverage-data gaps. "
    "This audit does not measure private vaccination or trial participation.", "",
    "Reproduce: `python3 R_cbh/covariates/09_vaccine_introduction_history.py`."]
(out / "VACCINE_INTRODUCTION_HISTORY.md").write_text("\n".join(lines)+"\n")
print(f"Saved introduction histories for {len(groups)//3} countries, including {len(study)} MAP-eligible countries.")

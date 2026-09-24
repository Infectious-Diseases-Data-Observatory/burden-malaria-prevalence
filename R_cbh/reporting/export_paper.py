#!/usr/bin/env python3
"""Export verified, plan-defined primary assets; never write a TeX file."""
import argparse
import csv
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument("--destination", type=Path, required=True)
parser.add_argument("--copy", action="store_true", help="Copy after validating; default is a local export manifest only")
parser.add_argument("--root", type=Path, default=Path("results/cbh/primary_map_regional17_dhsmics_gamma2_v5"))
parser.add_argument("--reference", type=Path, default=Path("results/cbh/primary_map_regional17_gamma2_v3"),
                    help="DHS-only primary that the DHS part of a DHS+MICS root must reproduce")
args = parser.parse_args()
root = args.root
out = root / "paper_refresh"
out.mkdir(parents=True, exist_ok=True)
destination = args.destination.expanduser().resolve()
assert destination.is_dir()
# Sensitivities refitted on the DHS and MICS sample; the DHS-only versions are kept as history.
subgroups = Path("results/cbh/subgroups_dhsmics_map_gamma2_v2")
nutrition = Path("results/cbh/nutrition_adjustment_dhsmics_map_gamma2_v2")
imputed = Path("results/cbh/primary_map_regional17_dhsmics_imputed_gamma2_v6")
manuscript_updates = Path("docs/MANUSCRIPT_UPDATE_DHS_MICS.md")

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def read_csv(path):
    with path.open() as f:
        return list(csv.DictReader(f))

figures = read_csv(root / "paper_figures/manifest.csv")
assert [int(x["figure"]) for x in figures] == list(range(1, 4))
assets = []
for row in figures:
    source = Path(row["source"])
    assert source.is_relative_to(root)
    assert hashlib.md5(source.read_bytes()).hexdigest() == row["md5"]
    assets.append((source, Path("figures") / row["filename"]))
assets.extend([
    (root / "paper_figures/CAPTIONS.md", Path("figures/primary_figure_captions.md")),
    (root / "burden_comparison/CAPTION.md", Path("figures/fig3_burden_comparison_caption.md")),
    (root / "burden_comparison/sfig_burden_comparison_counts.png", Path("Supplementary Figures/sfig_burden_comparison_counts.png")),
    (root / "burden_comparison/CAPTION_counts.md", Path("Supplementary Figures/sfig_burden_comparison_counts_caption.md")),
    (subgroups / "sfig_pfpr_splines_by_subgroup.png", Path("Supplementary Figures/sfig_pfpr_splines_by_subgroup.png")),
    (subgroups / "CAPTION.md", Path("Supplementary Figures/sfig_pfpr_splines_by_subgroup_caption.md")),
    (root / "study_flow/study_flow_diagram.png", Path("Supplementary Figures/sfig_study_flow.png")),
    (root / "pfpr_splines.png", Path("Supplementary Figures/sfig_pfpr_mortality_by_age.png")),
    (nutrition / "sfig_pfpr_splines_without_nutrition.png", Path("Supplementary Figures/sfig_pfpr_splines_without_nutrition.png")),
    (nutrition / "CAPTION.md", Path("Supplementary Figures/sfig_pfpr_splines_without_nutrition_caption.md")),
    (imputed / "sfig_pfpr_splines_imputed_covariates.png", Path("Supplementary Figures/sfig_pfpr_splines_imputed_covariates.png")),
    (imputed / "CAPTION_sfig_imputed_covariates.md", Path("Supplementary Figures/sfig_pfpr_splines_imputed_covariates_caption.md")),
    (root / "study_flow/CAPTION.md", Path("Supplementary Figures/sfig_study_flow_caption.md")),
])
for name in ("age_band_results.csv", "age_band_results.md", "age_band_results.latex.txt",
             "country_comparison_2024.csv", "country_comparison_2024_top10.md", "country_comparison_2024.latex.txt"):
    assets.append((root / "tables" / name, Path("tables") / name))
assets.extend([
    (root / "burden_comparison/source_comparison.latex.txt", Path("tables/source_comparison.latex.txt")),
    (root / "burden_comparison/source_comparison_2000_2024.csv", Path("tables/source_comparison_2000_2024.csv")),
    (root / "annual_comparison/annual_totals_2000_2024.csv", Path("tables/annual_totals_2000_2024.csv")),
    (root / "nigeria_states/state_totals_2024.csv", Path("tables/nigeria_state_totals_2024.csv")),
    (manuscript_updates, Path("MANUSCRIPT_TEXT_UPDATES.md")),
])
handoff = out / "OVERLEAF_UPDATE.md"
handoff_rel = Path("RESULTS_UPDATE.md")

def check(source, rel):
    assert source.is_file(), f"Missing export source: {source}"
    assert source.suffix.lower() != ".tex"
    assert not rel.is_absolute() and ".." not in rel.parts and rel.suffix.lower() != ".tex"
    assert (destination / rel).resolve().is_relative_to(destination)

# Validate every source before the handoff is written, so a failed export leaves no new handoff.
targets = [str(rel) for _, rel in assets] + [str(handoff_rel)]
assert len(set(targets)) == len(targets)
for source, rel in assets:
    check(source, rel)
# A sensitivity figure is current only if its report was rerun after the last change to its inputs or code
# (for example R_cbh/sensitivity/subgroups/02_report.R after an axis change).
for sensitivity in (subgroups, nutrition):
    stale = [row["file"] for row in read_csv(sensitivity / "report_provenance.csv")
             if not Path(row["file"]).is_file() or hashlib.md5(Path(row["file"]).read_bytes()).hexdigest() != row["md5"]]
    assert not stale, f"Stale report provenance in {sensitivity}: {', '.join(stale)}; rerun its report script"

sample =read_csv(root / "prepared_sample.csv")[0]
annual = next(row for row in read_csv(root / "annual_comparison/annual_totals_2000_2024.csv") if row["year"] == "2024")
n_covariates = len(read_csv(root / "covariate_scaling.csv"))
coverage = read_csv(root / "survey_map/survey_coverage.csv")
selection = {row["survey"]: row for row in read_csv(root / "prepared_selection_by_survey.csv")}
timeline = read_csv(root / "survey_map/survey_timeline_all.csv")
burden_countries = {row["iso3"] for row in read_csv(root / "annual_comparison/included_countries.csv")}
assert {row["type"] for row in coverage} <= {"DHS", "MICS"} and len(coverage) == int(sample["surveys"])

def part(kind):
    rows = [row for row in coverage if row["type"] == kind]
    return dict(surveys=len(rows), countries={row["country"] for row in rows},
                regions=sum(int(row["regions"]) for row in rows),
                records=sum(int(selection[row["survey"]]["records"]) for row in rows),
                deaths=sum(int(selection[row["survey"]]["deaths"]) for row in rows))

dhs, mics = part("DHS"), part("MICS")
assert dhs["records"] + mics["records"] == int(sample["records"])
assert dhs["deaths"] + mics["deaths"] == int(sample["deaths"])
assert dhs["regions"] + mics["regions"] == int(sample["regions"])
survey_countries = dhs["countries"] | mics["countries"]
assert len(survey_countries) == int(sample["countries"]) and int(annual["countries"]) == len(burden_countries)
fmt = lambda x: f"{float(x):,.0f}"
codes = lambda x: ", ".join(sorted(x))
mis_total = sum(row["type"] == "MIS" for row in timeline)
outside = survey_countries - burden_countries
if mics["surveys"]:
    reference = read_csv(args.reference / "prepared_sample.csv")[0]
    assert all(int(reference[k]) == dhs[k] for k in ("records", "deaths", "regions", "surveys")) and \
        int(reference["countries"]) == len(dhs["countries"]), "DHS part does not reproduce the DHS-only primary"
    added = mics["countries"] - dhs["countries"]
    composition = (
        f"It combines {dhs['surveys']} DHS and {mics['surveys']} UNICEF MICS surveys. "
        f"The DHS part reproduces the DHS-only primary `{args.reference.name}` exactly "
        f"({fmt(dhs['records'])} records, {fmt(dhs['deaths'])} deaths, {fmt(dhs['regions'])} survey-regions, "
        f"{dhs['surveys']} surveys, {len(dhs['countries'])} countries); the MICS surveys add "
        f"{fmt(mics['records'])} records, {fmt(mics['deaths'])} deaths, {fmt(mics['regions'])} survey-regions"
        + (f" and {len(added)} {'country' if len(added) == 1 else 'countries'} ({codes(added)})." if added else ".")
        + (f" None of the {mis_total} MIS surveys passes complete-case selection." if mis_total else ""))
    predictors = ("All 13 survey-level predictors are survey-region summaries: for DHS surveys, published "
                  "regional indicators or weighted recode summaries; for MICS surveys, values computed from "
                  "the microdata with the DHS definitions and the same WUENIC and within-survey fallbacks.")
else:
    composition = f"The analysis contains {dhs['surveys']} DHS surveys."
    predictors = "All 13 survey-level predictors are survey-region summaries."
overlap = (f"{len(survey_countries & burden_countries)} of the {len(survey_countries)} survey countries"
           + (f"; {codes(outside)} contributes survey data but is outside the burden set" if len(outside) == 1 else
              f"; {codes(outside)} contribute survey data but are outside the burden set" if outside else ""))
v6 = imputed / "prepared_sample.csv"
v6_text = ""
if v6.is_file():
    s6 = read_csv(v6)[0]
    flow = {row["item"]: int(row["count"]) for row in read_csv(root / "study_flow/flow_counts.csv")}
    retained = int(s6["records"]) == flow["eligible"] - flow["missing_map"]
    v6_text = (f"{', which retains every MAP-eligible record' if retained else ''}: {fmt(s6['records'])} records, {fmt(s6['deaths'])} deaths, "
               f"{s6['surveys']} surveys and {s6['countries']} countries")
title = "Updated primary results (DHS and MICS surveys)" if mics["surveys"] else "Updated primary results"
handoff.write_text(f"""# {title}

Figures 1–3 (Figure 3 combines the 2024 country comparison, the Nigerian state
comparison and the 2000–2024 annual trend as panels A–C, all as deaths per 1,000
under-five child-years; its death-count version is a supplementary figure) and the
supplementary inclusion flow use the {n_covariates}-variable
regional-adjustment MAP gamma=2 models in `{root.name}`. The analysis contains
{fmt(sample['children'])} children, {fmt(sample['records'])} child-band records and
{fmt(sample['deaths'])} deaths in {fmt(sample['regions'])} survey-regions,
{fmt(sample['surveys'])} surveys and {fmt(sample['countries'])} countries.

{composition}

The revised adjustment removes sex, multiple births and birth order; replaces
maternal age at each birth with regional mean age at first birth; and adds
regional wasting and stunting prevalence. {predictors} The four national annual
variables are unchanged.

No TeX file has been edited. Update manuscript text and captions manually using:

- [Manuscript text updates](MANUSCRIPT_TEXT_UPDATES.md) for the DHS and MICS
  sample (copied from `{manuscript_updates}` in the analysis project).
- [Age-band table source](tables/age_band_results.latex.txt), with death
  percentages summing to 100.0% and both prevalence contrasts.
- [Country table source](tables/country_comparison_2024.latex.txt).
- [Source comparison table source](tables/source_comparison.latex.txt) and its
  [values](tables/source_comparison_2000_2024.csv).
- [Main figure captions](figures/primary_figure_captions.md) and the
  [Figure 3 caption](figures/fig3_burden_comparison_caption.md); Figure 3 is
  `figures/fig3_burden_comparison.png`.
- [Inclusion-flow caption](<Supplementary Figures/sfig_study_flow_caption.md>).
- [Annual totals](tables/annual_totals_2000_2024.csv) and
  [Nigerian state totals](tables/nigeria_state_totals_2024.csv).

The 2024 totals are {fmt(annual['model_deaths'])} PfPR-ACM deaths,
{fmt(annual['ihme_malaria_deaths'])} IHME deaths and
{fmt(annual['who_cacode_deaths'])} UN IGME deaths, summed over the
{annual['countries']} countries in the national burden comparison ({overlap}).
These are conditional model-based estimates; source and imputation uncertainty
are not fully propagated.

The supplementary sensitivity figures are refitted on the DHS and MICS surveys:
subgroups (`{subgroups.name}`) and no nutrition covariates (`{nutrition.name}`)
on this sample, and imputed covariates on the larger MAP-eligible sample
(`{imputed.name}`{v6_text}). Their captions are exported to
`Supplementary Figures/`. The DHS-only sensitivity fits are kept in the analysis
project as history.

Previous figure versions are backed up in the analysis project. `.latex.txt`
files contain source for manual insertion; they do not change the compiled
manuscript automatically. Figure placement and TeX references remain under the
author's control.
""")
assets.append((handoff, handoff_rel))
check(handoff, handoff_rel)
rows = [dict(source=str(source.resolve()), destination=str(destination / rel),
             sha256=digest(source), previous_sha256=digest(destination / rel) if (destination / rel).is_file() else "")
        for source, rel in assets]
with (out / "export_manifest.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader(); writer.writerows(rows)
if not args.copy:
    print(f"Validated {len(assets)} assets for export; destination unchanged.")
    raise SystemExit(0)

tex_before = {str(p): digest(p) for p in destination.rglob("*.tex")}
backup = out / "previous_overleaf" / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
for source, rel in assets:
    target = destination / rel
    if target.is_file():
        saved = backup / rel
        saved.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(target, saved)
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, target)
    assert digest(target) == digest(source)
tex_after = {str(p): digest(p) for p in destination.rglob("*.tex")}
assert tex_before == tex_after, "TeX changed during export; inspect without overwriting the author's edits"
(out / "export_verification.json").write_text(json.dumps(dict(
    destination=str(destination), exported_assets=len(assets), all_copy_hashes_match=True,
    tex_unchanged_during_export=True, tex_hashes=tex_after, previous_asset_backup=str(backup),
    completed_at=datetime.now(timezone.utc).isoformat()), indent=2) + "\n")
print(f"Exported and verified {len(assets)} assets; {len(tex_after)} TeX files unchanged.")

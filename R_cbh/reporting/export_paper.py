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
args = parser.parse_args()
root = args.root
out = root / "paper_refresh"
out.mkdir(parents=True, exist_ok=True)
destination = args.destination.expanduser().resolve()
assert destination.is_dir()

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

with (root / "paper_figures/manifest.csv").open() as f:
    figures = list(csv.DictReader(f))
assert [int(x["figure"]) for x in figures] == list(range(1, 4))
assets = []
for row in figures:
    source = Path(row["source"])
    assert source.is_relative_to(root)
    assert hashlib.md5(source.read_bytes()).hexdigest() == row["md5"]
    assets.append((source, Path("figures") / row["filename"]))
assets.extend([
    (root / "paper_figures/CAPTIONS.md", Path("figures/primary_figure_captions.md")),
    (root / "burden_comparison/CAPTION.md", Path("figures/fig4_burden_comparison_caption.md")),
    (root / "burden_comparison/sfig_burden_comparison_counts.png", Path("Supplementary Figures/sfig_burden_comparison_counts.png")),
    (root / "burden_comparison/CAPTION_counts.md", Path("Supplementary Figures/sfig_burden_comparison_counts_caption.md")),
    (Path("results/cbh/subgroups_map_gamma2_v1/sfig_pfpr_splines_by_subgroup.png"), Path("Supplementary Figures/sfig_pfpr_splines_by_subgroup.png")),
    (Path("results/cbh/subgroups_map_gamma2_v1/CAPTION.md"), Path("Supplementary Figures/sfig_pfpr_splines_by_subgroup_caption.md")),
    (root / "study_flow/study_flow_diagram.png", Path("Supplementary Figures/sfig_study_flow.png")),
    (root / "pfpr_splines.png", Path("Supplementary Figures/sfig_pfpr_mortality_by_age.png")),
    (Path("results/cbh/nutrition_adjustment_map_gamma2_v1/sfig_pfpr_splines_without_nutrition.png"),
     Path("Supplementary Figures/sfig_pfpr_splines_without_nutrition.png")),
    (Path("results/cbh/nutrition_adjustment_map_gamma2_v1/CAPTION.md"),
     Path("Supplementary Figures/sfig_pfpr_splines_without_nutrition_caption.md")),
    (Path("results/cbh/primary_map_regional17_imputed_gamma2_v4/sfig_pfpr_splines_imputed_covariates.png"), Path("Supplementary Figures/sfig_pfpr_splines_imputed_covariates.png")),
    (Path("results/cbh/primary_map_regional17_imputed_gamma2_v4/CAPTION_sfig_imputed_covariates.md"), Path("Supplementary Figures/sfig_pfpr_splines_imputed_covariates_caption.md")),
    (root / "study_flow/CAPTION.md", Path("Supplementary Figures/sfig_study_flow_caption.md")),
])
for name in ("age_band_results.csv", "age_band_results.md", "age_band_results.latex.txt",
             "country_comparison_2024.csv", "country_comparison_2024_top10.md", "country_comparison_2024.latex.txt"):
    assets.append((root / "tables" / name, Path("tables") / name))
assets.extend([
    (root / "annual_comparison/annual_totals_2000_2024.csv", Path("tables/annual_totals_2000_2024.csv")),
    (root / "nigeria_states/state_totals_2024.csv", Path("tables/nigeria_state_totals_2024.csv")),
])
handoff = out / "OVERLEAF_UPDATE.md"
with (root / "prepared_sample.csv").open() as f:
    sample = next(csv.DictReader(f))
with (root / "annual_comparison/annual_totals_2000_2024.csv").open() as f:
    annual = next(row for row in csv.DictReader(f) if row["year"] == "2024")
with (root / "covariate_scaling.csv").open() as f:
    n_covariates = len(list(csv.DictReader(f)))
fmt = lambda x: f"{float(x):,.0f}"
handoff.write_text(f"""# Updated primary results

Figures 1–4 (Figure 4 combines the 2024 country comparison, the Nigerian state
comparison and the 2000–2024 annual trend as panels A–C, all as deaths per 1,000
under-five child-years; its death-count version is a supplementary figure) and the
supplementary inclusion flow use the {n_covariates}-variable
regional-adjustment MAP gamma=2 models in `{root.name}`. The analysis contains
{fmt(sample['children'])} children, {fmt(sample['records'])} child-band records and
{fmt(sample['deaths'])} deaths in {fmt(sample['regions'])} survey-regions,
{fmt(sample['surveys'])} surveys and {fmt(sample['countries'])} countries.

The revised adjustment removes sex, multiple births and birth order; replaces
maternal age at each birth with regional mean age at first birth; and adds
regional wasting and stunting prevalence. All retained DHS predictors are regional
summaries; the four existing national annual variables remain.

No TeX file has been edited. Update manuscript text and captions manually using:

- [Age-band table source](tables/age_band_results.latex.txt), with death
  percentages summing to 100.0% and both prevalence contrasts.
- [Country table source](tables/country_comparison_2024.latex.txt).
- [Main figure captions](figures/primary_figure_captions.md).
- [Inclusion-flow caption](<Supplementary Figures/sfig_study_flow_caption.md>).
- [Annual totals](tables/annual_totals_2000_2024.csv) and
  [Nigerian state totals](tables/nigeria_state_totals_2024.csv).

The 2024 totals are {fmt(annual['model_deaths'])} PfPR-ACM deaths,
{fmt(annual['ihme_malaria_deaths'])} IHME deaths and
{fmt(annual['who_cacode_deaths'])} UN IGME deaths across the same 42 countries.
These are conditional model-based estimates; source and imputation uncertainty
are not fully propagated. Sensitivity fits remain historical and are not
relabelled as this revised primary specification.

Previous figure versions are backed up in the analysis project. `.latex.txt`
files contain source for manual insertion; they do not change the compiled
manuscript automatically. Figure placement and TeX references remain under the
author's control.
""")
assets.append((handoff, Path("RESULTS_UPDATE.md")))
assert len({str(rel) for _, rel in assets}) == len(assets)
for source, rel in assets:
    assert source.is_file() and source.suffix.lower() != ".tex"
    assert not rel.is_absolute() and ".." not in rel.parts and rel.suffix.lower() != ".tex"
    assert (destination / rel).resolve().is_relative_to(destination)
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

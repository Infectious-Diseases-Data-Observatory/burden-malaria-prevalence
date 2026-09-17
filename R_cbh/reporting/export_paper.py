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
args = parser.parse_args()
root = Path("results/cbh/primary_map_regional18_gamma2_v2")
out = root / "paper_refresh"
out.mkdir(parents=True, exist_ok=True)
destination = args.destination.expanduser().resolve()
assert destination.is_dir()

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

with (root / "paper_figures/manifest.csv").open() as f:
    figures = list(csv.DictReader(f))
assert [int(x["figure"]) for x in figures] == list(range(1, 7))
assets = []
for row in figures:
    source = Path(row["source"])
    assert source.is_relative_to(root)
    assert hashlib.md5(source.read_bytes()).hexdigest() == row["md5"]
    assets.append((source, Path("figures") / row["filename"]))
assets.extend([
    (root / "paper_figures/CAPTIONS.md", Path("figures/primary_figure_captions.md")),
    (root / "annual_comparison/CAPTION.md", Path("figures/fig5_annual_malaria_mortality_caption.md")),
    (root / "nigeria_states/CAPTION.md", Path("figures/fig6_nigeria_states_vs_ihme_caption.md")),
    (root / "study_flow/study_flow_diagram.png", Path("Supplementary Figures/sfig_study_flow.png")),
    (root / "study_flow/CAPTION.md", Path("Supplementary Figures/sfig_study_flow_caption.md")),
])
for name in ("age_band_results.csv", "age_band_results.md", "age_band_results.latex.txt",
             "country_comparison_2024.csv", "country_comparison_2024_top10.md", "country_comparison_2024.latex.txt"):
    assets.append((root / "tables" / name, Path("tables") / name))
assets.extend([
    (root / "annual_comparison/annual_totals_2004_2024.csv", Path("tables/annual_totals_2004_2024.csv")),
    (root / "nigeria_states/state_totals_2024.csv", Path("tables/nigeria_state_totals_2024.csv")),
])
handoff = out / "OVERLEAF_UPDATE.md"
handoff.write_text("""# Updated primary results

Figures 1–6 and the supplementary inclusion flow use the revised 18-variable
regional-adjustment MAP gamma=2 models. The analysis contains 1,755,838 children,
5,680,117 child-band records and 78,634 deaths in 973 survey-regions, 100 surveys
and 34 countries. All 100 selected surveys are DHS surveys.

No TeX file has been edited. The tables embedded in `main.tex`, figure captions
and narrative numbers therefore still need the author's updates:

- [Age-band table source](tables/age_band_results.latex.txt), with death
  percentages summing to 100.0% and both prevalence contrasts.
- [Country table source](tables/country_comparison_2024.latex.txt), using UN IGME
  rather than the previous WHO proxy. The total covers all 42 countries, even
  though only the ten largest absolute model–IHME differences are shown.
- [Main figure captions](figures/primary_figure_captions.md).
- [Inclusion-flow caption](<Supplementary Figures/sfig_study_flow_caption.md>).
- [Annual totals](tables/annual_totals_2004_2024.csv) and
  [Nigerian state totals](tables/nigeria_state_totals_2024.csv).

The 2024 totals are 501,065 PfPR-ACM deaths, 428,147 IHME deaths and 440,123
UN IGME deaths across the same 42 countries. The corresponding rates are 282.4,
241.3 and 248.0 deaths per 100,000 under-five child-years.

Figure 6 is provided as specified in the plan; the author controls whether and
where it is included in TeX. Existing sensitivity and historical covariate-forest
figures have not been refitted or relabelled. They are outside this requested
primary reporting refresh and should not be described as results from the
revised primary adjustment. The conceptual causal-diagram assets are unchanged.

The source outputs, provenance and previous-version backups are retained in
the analysis project. `.latex.txt` files contain source for manual insertion;
they do not automatically change the compiled manuscript.
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

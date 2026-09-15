# Analysis-plan/code audit outputs

Read the [component audit](../../../docs/CODE_AUDIT.md) with the [analysis plan](../../../docs/ANALYSIS_PLAN.md).

- `code_inventory.csv`: per-file disposition, reason, original hash and archive/current location. For replaced entry points, the location records the preserved original; the replacement remains at the original path. The inventory includes the newly added audit verifier and results-index script.
- `archive_manifest.csv`: 68 retired scripts plus three preserved/replaced entry points; original hashes are unchanged.
- `primary_fit_checks.csv`, `primary_sample.csv`: independent exact primary model-frame checks and distinct-child counts.
- `burden_checks.csv`: saved national accounting and ages2–4 allocation checks.
- `verification_provenance.csv`, `session_info.txt`: model/data/source hashes and audit runtime.
- `parse_checks.csv`, `literal_source_checks.csv`, `CLEANUP_VALIDATION.md`: cleanup verification and its limits.

Repeat the saved-model checks with `Rscript run_all.R --audit`. This reads local prepared data and fitted objects, exports only aggregates, and does not fit models.

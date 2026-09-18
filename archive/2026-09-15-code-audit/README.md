# Historical code archived on 15 September 2026

These are the superseded country, survey-region and six-band person-time analyses, plus the previous run-all scripts and historical figure copier. They do not implement the seven-band MAP gamma=2 primary analysis in the current plan.

Files were moved without editing their contents. Original paths, reasons and SHA-256 hashes are in the [archive manifest](../../results/cbh/code_audit_2026_09_15/archive_manifest.csv). The [full inventory](../../results/cbh/code_audit_2026_09_15/code_inventory.csv) also explains what was retained and why. Existing historical results, fitted objects and raw data were not moved or deleted.

The archived code preserves its original project-root paths and is a reference snapshot, not a runnable pipeline from this directory. For historical reproduction, use an isolated checkout at commit `f44fa326b58735f6266315398be154233da09138` with the corresponding authorized input snapshots. Do not copy archived scripts back over the active workflow. Some historical scripts can fetch external data, fit models or copy figures outside this project; archiving does not endorse those actions.

The primary CBH builder, HIV model, seven MAP gamma=2 fits, Snow supplementary workflows (archived separately on 18 September 2026 in `../2026-09-18-snow-comparison/`), input bridges and their shared helpers remain active. Historical CBH fitting/burden code with still-needed data/model-frame dependencies was deliberately retained and is documented in [the audit](../../docs/CODE_AUDIT.md).

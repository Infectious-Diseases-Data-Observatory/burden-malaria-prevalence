# =============================================================================
# tests/test_rkey.R
#
# Guard the region-key normaliser in R_dhs/00_config.R.
#
# rkey() joins DHS Births Recode region labels (v024) to DHS boundary labels
# (DHSREGEN). Some recodes reach us as Latin-1 text decoded as Windows-1250, so
# accented letters arrive as the Central-European letter sharing the same byte
# and the join silently drops the region. rkey() repairs that. These checks
# pin the repair down and, just as importantly, pin down that it is a no-op for
# correctly encoded names.
#
#   Rscript tests/test_rkey.R
# =============================================================================

source("R_dhs/00_config.R")

failures <- 0L
check <- function(label, ok) {
  ok <- isTRUE(ok)
  message(sprintf("%-58s %s", label, if (ok) "ok" else "FAIL"))
  if (!ok) failures <<- failures + 1L
  invisible(ok)
}

# --- the corrupted names must reconcile with their boundary spellings --------
# DR Congo 2007 and Senegal 2005/2008/2010/2012/2014 are the affected surveys.
check("Kasai Occident: CP1250 d-caron reconciles with boundary",
      identical(rkey("kasaď occident"), rkey("Kasai Occident")))
check("Kasai Oriental: CP1250 d-caron reconciles with boundary",
      identical(rkey("kasaď oriental"), rkey("Kasai Oriental")))
check("Thies: CP1250 c-caron reconciles with boundary",
      identical(rkey("thičs"), rkey("Thies")))
# DHS spells Thies with e-grave in some rounds and e-acute in others; both must
# land on the same key as the unaccented boundary label.
check("Thies: e-grave, e-acute and plain ASCII all agree",
      length(unique(rkey(c("Thiès", "Thiés", "Thies")))) == 1L)

# --- the repair must not touch correctly encoded names -----------------------
# These are real region names from Angola, Gabon, Guinea, Cameroon, Sao Tome,
# Mozambique, Niger and Chad recodes, all of which decode correctly today.
correct <- c("Huíla", "Bié", "Ogooué-Lolo", "N'zérékoré",
             "Extrême-Nord", "São Tomé", "Côte d'Ivoire",
             "Zambézia", "Tillabéri", "Ouaddaï", "N'Djaména",
             "Mohéli", "Ségou", "Boké")
expected <- c("huila", "bie", "ogoouelolo", "nzerekore", "extremenord", "saotome",
              "cotedivoire", "zambezia", "tillaberi", "ouaddai", "ndjamena",
              "moheli", "segou", "boke")
check("correctly encoded names are unchanged by the repair",
      identical(rkey(correct), expected))

# --- structural properties ---------------------------------------------------
check("output is restricted to [a-z0-9]",
      all(grepl("^[a-z0-9]*$", rkey(c(correct, "Saint-Louis", "Ville d'Abidjan")))))
check("rkey is idempotent",
      identical(rkey(rkey(correct)), rkey(correct)))
check("empty and NA inputs do not error",
      length(rkey(c("", NA_character_))) == 2L)
# A character iconv cannot transliterate must degrade the key, not void it:
# returning NA here would drop the region silently, which is the very failure
# this normaliser exists to prevent.
check("an untransliterable character degrades rather than voiding the key",
      !is.na(rkey("Bo\u043ae")) && nzchar(rkey("Bo\u043ae")))

# --- the repair must not collapse two distinct regions into one key ----------
# Every region name in the legacy aggregate is already a clean key, so applying
# rkey() to it must be the identity: if this ever fails, the change would break
# the merge that builds the analysis panel.
aggregate_path <- file.path(REPO_ROOT, "results",
                            "component2_region_data_full.csv")
if (file.exists(aggregate_path)) {
  keys <- read.csv(aggregate_path, stringsAsFactors = FALSE)$regkey
  check("rkey() is the identity on every key in the legacy aggregate",
        identical(rkey(keys), keys))
} else {
  message("legacy aggregate absent; skipping the identity check")
}

if (failures > 0L) {
  stop(sprintf("%d rkey() check(s) failed", failures))
}
message("\nAll rkey() checks passed.")

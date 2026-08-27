# =============================================================================
# tests/test_region_match.R
#
# Guard match_region_keys() in R_dhs/00_config.R, which reconciles a DHS Births
# Recode's region labels with the DHS boundary's before the two are merged.
#
# The function exists because exact key matching silently dropped 199 region
# rows across the panel. It must recover the lexical mismatches WITHOUT ever
# inventing a pairing: an unmatched region is a visible gap, a wrong pairing
# attaches one region's mortality to another region's prevalence.
#
#   Rscript tests/test_region_match.R
# =============================================================================

source("R_dhs/00_config.R")

failures <- 0L
check <- function(label, ok) {
  ok <- isTRUE(ok)
  message(sprintf("%-62s %s", label, if (ok) "ok" else "FAIL"))
  if (!ok) failures <<- failures + 1L
}

# --- exact matching is preserved unchanged -----------------------------------
identity_map <- match_region_keys(c("north", "south"), c("south", "north"))
check("all-exact returns an identity crosswalk",
      nrow(identity_map) == 2L && all(identity_map$how == "exact") &&
        all(identity_map$from == identity_map$to))
check("empty and blank keys are ignored",
      nrow(match_region_keys(c("", NA_character_), c("north"))) == 0L)

# --- the real mismatches it must repair --------------------------------------
# Malawi 2000: recode says north/south, the boundary says Northern/Southern.
malawi <- match_region_keys(c("north", "central", "south"),
                            c("northern", "central", "southern"))
check("short form reconciles with long form (Malawi 2000)",
      nrow(malawi) == 3L && malawi$to[malawi$from == "north"] == "northern" &&
        malawi$to[malawi$from == "south"] == "southern")
# Cameroon 2004: older recodes truncate the label to twelve characters.
truncated <- match_region_keys(c("extremenor", "nord"),
                               c("extremenord", "nord"))
check("truncated recode label reconciles (Cameroon 2004)",
      truncated$to[truncated$from == "extremenor"] == "extremenord")
# Malawi 2015: the recode appends a generic noun the boundary omits.
noun <- match_region_keys(c("centralregion", "northernregion"),
                          c("central", "northern"))
check("appended generic noun reconciles (Malawi 2015)",
      nrow(noun) == 2L && noun$to[noun$from == "centralregion"] == "central")

# --- what it must REFUSE -----------------------------------------------------
# Where a boundary polygon covers several recode regions, a prefix pairing would
# attach one region's mortality to a polygon spanning several. Unequal unit
# counts are the signal, and they must veto prefix matching outright.
check("unequal counts refuse prefix pairing (Mali 2001)",
      nrow(match_region_keys(c("kidal", "gao", "timbouctou"),
                             c("kidalgaotimbouctou"))) == 0L)
check("unequal counts refuse prefix pairing (Tanzania 2004)",
      nrow(match_region_keys(c("zanzibarnorth", "zanzibarsouth"),
                             c("zanzibar"))) == 0L)
# One recode key with two plausible boundary keys is a guess, not a match.
ambiguous <- match_region_keys(c("north", "northx"),
                               c("northern", "northwest"))
check("a recode key with two candidates is refused",
      sum(ambiguous$how == "prefix") == 0L)
# And the same in the other direction: two recode keys competing for one
# boundary key must leave BOTH unmatched rather than awarding it to either.
contested <- match_region_keys(c("north", "northe"),
                               c("northern", "northeast"))
check("two recode keys contesting one boundary key are both refused",
      sum(contested$how == "prefix") == 0L)
check("a prefix shorter than the minimum is refused",
      nrow(match_region_keys("ab", "abcdef")) == 0L)

# --- elimination is allowed, guessing is not ---------------------------------
# With equal counts and every other region matched exactly, the last pair is
# determined rather than guessed, so it is accepted.
forced <- match_region_keys(c("kasai", "kasaioriental"),
                            c("kasaioccidental", "kasaioriental"))
check("the last remaining pair is settled by elimination",
      nrow(forced) == 2L && forced$to[forced$from == "kasai"] == "kasaioccidental")

# --- the crosswalk must never be many-to-one ---------------------------------
# Two recode regions mapping onto one boundary region would double-count that
# polygon's prevalence.
many <- match_region_keys(c("centre", "centresansouagadougou", "nord"),
                          c("centre", "ouagadougou", "nord"))
check("no two recode keys map onto the same boundary key",
      !any(duplicated(many$to)))

if (failures > 0L) stop(sprintf("%d match_region_keys() check(s) failed", failures))
message("\nAll match_region_keys() checks passed.")

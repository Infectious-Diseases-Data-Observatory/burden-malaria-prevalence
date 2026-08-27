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
# Two units a side, so the elimination stage cannot rescue a pairing the prefix
# rule declines: a fragment shorter than the minimum is not a prefix match.
check("a prefix shorter than the minimum is refused",
      nrow(match_region_keys(c("ab", "zzz"), c("abcdef", "yyy"))) == 0L)

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

# --- language, word order and qualifiers -------------------------------------
# These run on RAW labels, not collapsed keys: word order and language can only
# be normalised while the word boundaries survive.
sig_agree <- function(a, b) identical(region_signature(a), region_signature(b))

check("French direction words fold onto English (Rwanda 2008)",
      sig_agree("nord", "North") && sig_agree("ouest", "West") &&
        sig_agree("est", "East") && sig_agree("sud", "South"))
check("Portuguese direction words fold onto English (Sao Tome 2008)",
      sig_agree("Regi\u00e3o Norte", "North") &&
        sig_agree("Regi\u00e3o Sul", "South"))
check("a generic noun in Portuguese is dropped (Sao Tome 2008)",
      sig_agree("Regi\u00e3o do Pr\u00edncipe", "Principe"))
check("word order is ignored (Angola 2023)",
      sig_agree("Cuanza Norte", "North Cuanza") &&
        sig_agree("Lunda Sul", "South Lunda"))
check("written-together compass compounds split (Cote d'Ivoire 2012)",
      sig_agree("nord-est", "Northeast") && sig_agree("centre-nord", "North Central"))
check("an exclusion clause is dropped in either language (Togo 2013)",
      sig_agree("Maritime (sans agglom\u00e9ration de Lom\u00e9)",
                "Maritime excluding Greater Lom\u00e9 area"))
check("et and and are both dropped (Mauritania 2020)",
      sig_agree("Tiris Zemour et Inchiri", "Tiris Zemour and Inchiri"))
check("city survives while province does not (Mozambique 2003)",
      sig_agree("Maputo Cidade", "Maputo City") &&
        sig_agree("Maputo Prov\u00edncia", "Maputo") &&
        !sig_agree("Maputo Cidade", "Maputo"))

# An accented word must stay ONE token. macOS iconv TRANSLIT writes the
# diacritic as a separate ASCII character, which would split it in two.
check("an accented word is not split into two tokens",
      identical(region_tokens("Regi\u00e3o Centro"), "central"))

# --- the full crosswalk on real cases ----------------------------------------
ivory <- match_region_keys(
  c("centre", "centre-est", "centre-nord", "centre-ouest", "nord", "nord-est",
    "nord-ouest", "ouest", "sud-ouest", "sud sans abidjan", "ville d'abidjan"),
  c("Central", "East Central", "North Central", "West Central", "North",
    "Northeast", "Northwest", "West", "South", "Southwest", "Abidjan")
)
check("Cote d'Ivoire 2012 reconciles completely",
      nrow(ivory) == 11L && !any(duplicated(ivory$to)))
check("Cote d'Ivoire 2012 maps centre-nord to North Central, not Northeast",
      ivory$to[ivory$from == "centrenord"] == "northcentral")

# Spelling slips must be repaired, but only within a tight edit budget.
check("a one-letter spelling slip is repaired (Zimbabwe 2005)",
      match_region_keys(c("Matebeleland North", "Harare"),
                        c("Matabeleland North", "Harare"))$to[1] %in%
        c("matabelelandnorth", "harare"))
check("two unrelated names are not fuzzed together",
      nrow(match_region_keys(c("Kayes", "Mopti"), c("Kidal", "Segou"))) == 0L)

# A curated synonym is an explicit fact, so it applies even at unequal counts.
check("a curated rename applies (Namibia: Caprivi became Zambezi)",
      "zambezi" %in% match_region_keys(c("Caprivi", "Hardap", "Karas"),
                                       c("Zambezi", "Hardap"))$to)

# --- granularity must still be refused ---------------------------------------
# Tanzania 2004 reports 26 regions against an 8-zone boundary. The counts must
# differ for this to be the granularity case it really is.
check("many regions against few zones stay unmatched (Tanzania 2004)",
      nrow(match_region_keys(c("Arusha", "Dodoma", "Zanzibar North", "Mbeya"),
                             c("Northern", "Central", "Zanzibar"))) == 0L)
# And the Senegal continuous series: 14 regions against a 4-zone boundary.
check("14 regions against 4 zones stay unmatched (Senegal 2014)",
      nrow(match_region_keys(c("dakar", "thies", "kolda", "matam", "louga"),
                             c("North", "South", "West", "Center"))) == 0L)

if (failures > 0L) stop(sprintf("%d match_region_keys() check(s) failed", failures))
message("\nAll match_region_keys() checks passed.")

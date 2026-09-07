# Survey rules and geographic overrides

`survey_rules.csv` declares all 124 local surveys. Region-variable choices were migrated from the existing survey merge audit, then checked against the local recode schemas. Older Ethiopian surveys use +92 CMC; ET8AFL already uses Gregorian dates. Empty stratum choices mean prefer an observed V022, otherwise V023. This preserves design metadata, not a completed survey-design review.

SN6DFL/SN71FL use SN7AFL's V024-to-SZONE grouping; TZ4IFL uses TZ63FL's V024-to-V023 grouping. The builder verifies that each donor's fine area has only one coarse group. Donor files are included in the cache signature.

`region_overrides.csv` records survey-specific label equivalences reviewed against the local recode and cached MAP boundary labels. Most are spelling, language or abbreviated-name differences. GA41FL's regional labels spell out constituent provinces in parentheses. MD81FL explicitly separates Antananarivo from the other 22 regions; its Analamanga value maps to the boundary excluding the separately recorded capital. TZ4IFL's two spelling differences refer to the declared donor's Eastern and Zanzibar groups.

These are reproducible local crosswalk decisions, not evidence that regional boundaries are stable across surveys. All `region_id` fields currently remain empty. A nonempty value would explicitly assert that a region can share an identifier across surveys; review boundary equivalence before adding one.

Unverified historical/coarse equivalences are deliberately left unmatched, including the older Benin groupings, Burundi's Bujumbura group, Mali's northern group, Chad's Ennedi group, Uganda's South/Central groupings and Zimbabwe's Harare/Chitungwiza grouping. Their unmatched source labels appear in the generated crosswalk; candidate targets are in the configured boundary table. Some other regions have no corresponding entry in the cached MAP boundary/exposure table. Missing geography or exposure is never converted to zero.

Exclusion clauses also remain meaningful: a label such as “Centre excluding Yaoundé” is not automatically reduced to “Centre.” Where the boundary table uses a shorter name, establish that its polygon covers the same population before adding an explicit override. The output retains these rows for review and marks them outside the default model input.

Do not repair a missing region with a nearest name or by elimination. Check the intended survey geography and update this explicit table when the correspondence is established.

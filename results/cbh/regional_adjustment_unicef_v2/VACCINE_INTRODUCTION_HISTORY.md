# Hib, PCV and rotavirus introduction histories

WHO annual introduction-status records, retrieved 17 September 2026. This is a descriptive source audit; no missing coverage values or model inputs were changed.

A single year means the first reported nationwide introduction. `P year → year` means first reported partial introduction followed by first reported nationwide introduction. These are annual reporting dates, not necessarily exact launch dates, and nationwide introduction does not mean 100% coverage. `No through 2025` means no introduction recorded in this snapshot.

WHO cautions that introduction can precede the first year of consistent reporting. Nigeria PCV is a verified example: the [WHO launch announcement](https://www.afro.who.int/countries/nigeria/news/nigeria-introduces-new-vaccine-pcv-10) dates launch to 22 December 2014, whereas the annual table first records partial introduction in 2015 and nationwide introduction in 2017. The raw annual values are preserved below. No other dates have been individually reconciled with launch announcements.

Source: [WHO introduction portal](https://immunizationdata.who.int/global/wiise-detail-page/vaccine-introduction-in-country_name). CSV contains each country's direct source link and separate partial/nationwide fields. The 48-country scope excludes Algeria and Sudan; it includes Djibouti and Somalia. The study column marks the 36 countries with MAP-eligible records before expanded covariate selection.

| Country | Study | Hib | PCV | Rotavirus |
|---|:---:|---|---|---|
| Angola | Yes | 2006 | 2013 | 2014 |
| Benin | Yes | 2005 | 2011 | 2019 |
| Botswana | — | 2011 | 2012 | 2012 |
| Burkina Faso | Yes | 2006 | 2013 | 2013 |
| Burundi | Yes | 2004 | 2011 | 2013 |
| Cabo Verde | — | 2010 | No through 2025 | No through 2025 |
| Cameroon | Yes | 2009 | 2011 | 2014 |
| Central African Republic | — | 2008 | 2011 | No through 2025 |
| Chad | Yes | 2008 | 2024 | 2024 |
| Comoros | Yes | 2009 | No through 2025 | No through 2025 |
| Congo | Yes | 2009 | 2012 | 2014 |
| Côte d'Ivoire | Yes | 2009 | 2014 | 2017 |
| Democratic Republic of the Congo | Yes | 2009 | P 2011 → 2013 | 2019 |
| Djibouti | — | 2007 | 2012 | 2014 |
| Equatorial Guinea | — | 2013 | No through 2025 | No through 2025 |
| Eritrea | — | 2008 | 2015 | 2014 |
| Eswatini | Yes | 2009 | 2014 | 2015 |
| Ethiopia | Yes | 2007 | 2011 | P 2013 → 2014 |
| Gabon | Yes | 2010 | No through 2025 | No through 2025 |
| Gambia | Yes | P 1994 → 1997 | 2009 | 2013 |
| Ghana | Yes | P 2001 → 2002 | 2012 | 2012 |
| Guinea | Yes | 2008 | No through 2025 | No through 2025 |
| Guinea-Bissau | — | 2008 | 2015 | 2012 |
| Kenya | Yes | 2001 | 2011 | 2014 |
| Lesotho | — | 2008 | 2015 | 2017 |
| Liberia | Yes | 2008 | 2014 | 2016 |
| Madagascar | Yes | 2008 | 2012 | 2014 |
| Malawi | Yes | 2002 | 2011 | 2012 |
| Mali | Yes | P 2005 → 2007 | 2011 | P 2014 → 2015 |
| Mauritania | Yes | 2009 | 2013 | 2014 |
| Mauritius | — | 2006 | 2016 | 2015 |
| Mozambique | Yes | 2009 | 2013 | 2015 |
| Namibia | Yes | 2009 | 2014 | 2014 |
| Niger | Yes | 2008 | 2014 | 2014 |
| Nigeria | Yes | P 2012 → 2013 | P 2015 → 2017 | 2022 |
| Rwanda | Yes | 2002 | 2009 | 2012 |
| Sao Tome and Principe | Yes | 2009 | 2012 | 2016 |
| Senegal | Yes | 2005 | 2013 | 2014 |
| Seychelles | — | 2010 | 2018 | 2017 |
| Sierra Leone | Yes | 2007 | 2011 | 2014 |
| Somalia | — | 2013 | 2025 | 2025 |
| South Africa | Yes | P 1998 → 1999 | P 2008 → 2009 | P 2008 → 2009 |
| South Sudan | — | 2014 | 2025 | 2025 |
| Togo | Yes | 2008 | 2014 | 2014 |
| Uganda | Yes | 2002 | P 2013 → 2014 | 2018 |
| United Republic of Tanzania | Yes | 2009 | 2013 | 2013 |
| Zambia | Yes | 2004 | 2013 | P 2012 → 2013 |
| Zimbabwe | Yes | 2008 | 2012 | 2014 |

For analysis, documented pre-introduction years are candidates for structural zero routine-programme coverage. Partial-rollout and launch years must not automatically receive zero; missing values during or after introduction remain coverage-data gaps. This audit does not measure private vaccination or trial participation.

Reproduce: `python3 R_cbh/covariates/09_vaccine_introduction_history.py`.

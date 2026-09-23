# Inputs the MICS scripts read from outside this repository. Sourced by 04, 05 and 06.
# Africa admin-1 shapefile from the Snow prevalence project: the source of the Guinea-Bissau
# and Central African Republic analysis regions and of the SNOW: boundaries (Senegal, Somalia,
# South Sudan). Set MICS_ADMIN1_SHP to use a copy elsewhere.
mics_admin1_shp <- Sys.getenv("MICS_ADMIN1_SHP",
  "/Users/jameswatson/Documents/Claude Projects/Prevalence Model/Snow Prevalence Data/shape files/Africa_New_Admin.shp")

# Source from the project root, then use cbh_config(), cbh_build(), or cbh_load_analysis().
cbh_source <- function(root = getwd(), envir = parent.frame()) {
  paths <- c("00_config.R", file.path("R", c("utils.R", "geography.R", "child_bands.R",
                                            "external_data.R", "build.R", "load.R")))
  for (path in paths) sys.source(file.path(root, "R_cbh", path), envir = envir)
  invisible(NULL)
}
cbh_source()

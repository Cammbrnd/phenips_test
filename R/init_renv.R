# =============================================================================
# A lancer UNE FOIS, en local, depuis la racine du depot, pour generer renv.lock
# -----------------------------------------------------------------------------
#   source("R/init_renv.R")
#
# Ensuite : commiter renv.lock (et renv/activate.R, .Rprofile). Le workflow
# detecte renv.lock et installe EXACTEMENT ces versions (barrks, terra, ...)
# au lieu de prendre la derniere version CRAN a chaque run.
#
# Pour mettre a jour volontairement un package plus tard :
#   renv::update("barrks"); renv::snapshot()
# =============================================================================

if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv")

renv::init(bare = TRUE, restart = FALSE)

pkgs <- c("barrks", "terra", "httr2", "ggplot2", "tidyterra", "lubridate",
          "glue", "sf", "scales", "data.table", "dplyr", "jsonlite", "PHENTHAUproc")
renv::install(pkgs)

# type = "all" : inscrit aussi les packages appeles uniquement via library()
# dans les scripts de R/ (qui ne sont pas un package)
renv::snapshot(type = "all", prompt = FALSE)

cat("\nrenv.lock genere. A commiter : renv.lock, renv/activate.R, .Rprofile, renv/settings.json\n")

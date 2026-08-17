# =============================================================================
# PIPELINE STENO — reutilise le jeu de donnees SAFRAN partage
# -----------------------------------------------------------------------------
# R/steno_functions.R est deja allege. Le script original produisait 3 GeoJSON
# (generations, essaimage, grille) : ici grille = FALSE pour rester a 2
# fichiers, comme les 3 autres modeles.
# =============================================================================

source("R/safran_common.R")
source("R/steno_functions.R")   # garde uniquement le modele + les exports GeoJSON

meteo_wrap <- readRDS("data/cache_meteo/meteo_national.rds")
meteo      <- lapply(meteo_wrap, terra::unwrap)

steno <- modeliser_steno(meteo, exposition = "sunny", k_ecorce = 0.2)

exporter_geojson(steno, pas_essaimage = 15, grille = FALSE)

copier_vers_nom_stable <- function(pattern, nom_stable) {
  fichiers <- list.files(DIR_OUT, pattern = pattern, full.names = TRUE)
  if (length(fichiers) == 0) return(invisible(NULL))
  plus_recent <- fichiers[which.max(file.mtime(fichiers))]
  file.copy(plus_recent, file.path(DIR_OUT, nom_stable), overwrite = TRUE)
  cat(glue::glue("  Copie stable : {basename(plus_recent)} -> {nom_stable}\n"))
}
copier_vers_nom_stable("^steno_generations_.*\\.geojson$", "steno_generations_latest.geojson")
copier_vers_nom_stable("^steno_essaimage_.*\\.geojson$",   "steno_essaimage_latest.geojson")

cat("\n>> STENO termine.\n\n")

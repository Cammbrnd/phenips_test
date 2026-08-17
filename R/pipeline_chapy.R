# =============================================================================
# PIPELINE CHAPY — reutilise le jeu de donnees SAFRAN partage
# -----------------------------------------------------------------------------
# R/chapy_functions.R est deja allege (fonctions communes retirees, lancement
# automatique neutralise). Garde modeliser_chapy, produire_geojson_generations_chapy,
# produire_geojson_envol_chapy et tous les helpers/cartes/courbes du fichier
# original (non utilises ici mais inoffensifs a laisser).
# =============================================================================

source("R/safran_common.R")
source("R/chapy_functions.R")   # garde uniquement le modele + les exports GeoJSON

meteo_wrap <- readRDS("data/cache_meteo/meteo_national.rds")
meteo      <- lapply(meteo_wrap, terra::unwrap)

pheno <- modeliser_chapy(meteo, mode = "max")

produire_geojson_generations_chapy(pheno)
produire_geojson_envol_chapy(pheno)

# Copies a nom stable pour lecture par l'app Shiny (meme logique que phenips)
copier_vers_nom_stable <- function(pattern, nom_stable) {
  fichiers <- list.files(DIR_OUT, pattern = pattern, full.names = TRUE)
  if (length(fichiers) == 0) return(invisible(NULL))
  plus_recent <- fichiers[which.max(file.mtime(fichiers))]
  file.copy(plus_recent, file.path(DIR_OUT, nom_stable), overwrite = TRUE)
  cat(glue::glue("  Copie stable : {basename(plus_recent)} -> {nom_stable}\n"))
}
copier_vers_nom_stable("^chapy_generations_.*\\.geojson$", "chapy_generations_latest.geojson")
copier_vers_nom_stable("^chapy_onset_.*\\.geojson$",        "chapy_onset_latest.geojson")

cat("\n>> CHAPY termine.\n\n")

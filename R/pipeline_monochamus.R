# =============================================================================
# PIPELINE MONOCHAMUS — reutilise le jeu de donnees SAFRAN partage
# -----------------------------------------------------------------------------
# R/monochamus_functions.R est deja allege ET patche : le fichier [1]
# (mono_stade) porte desormais, en plus du stade gen1, le bivoltisme [A]
# (risque_bivoltA 0/1 + dj10_projete), extrait par centroide depuis
# pheno$risque_bivolt et pheno$dj_projete. Le fichier [2] (mono_bivoltB_gen2)
# reste le bivoltisme [B]. Toujours 2 GeoJSON au total.
# =============================================================================

source("R/safran_common.R")
source("R/monochamus_functions.R")   # garde uniquement le modele + exports GeoJSON (patchee)

meteo_wrap <- readRDS("data/cache_meteo/meteo_national.rds")
meteo      <- lapply(meteo_wrap, terra::unwrap)

pheno <- calculer_dj(meteo)

exporter_geojson(pheno)

copier_vers_nom_stable <- function(pattern, nom_stable) {
  fichiers <- list.files(DIR_OUT, pattern = pattern, full.names = TRUE)
  if (length(fichiers) == 0) return(invisible(NULL))
  plus_recent <- fichiers[which.max(file.mtime(fichiers))]
  file.copy(plus_recent, file.path(DIR_OUT, nom_stable), overwrite = TRUE)
  cat(glue::glue("  Copie stable : {basename(plus_recent)} -> {nom_stable}\n"))
}
copier_vers_nom_stable("^mono_stade_.*\\.geojson$",       "mono_stade_latest.geojson")
copier_vers_nom_stable("^mono_bivoltB_gen2_.*\\.geojson$", "mono_bivoltB_gen2_latest.geojson")

cat("\n>> Monochamus (avec bivoltisme A et B) termine.\n\n")

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

# [v3] On ne garde que les mailles SAFRAN qui croisent l'essence-hote (BD Foret IGN)
source("R/masques.R")
meteo <- appliquer_masque(meteo, charger_masque("epicea", meteo$tmean), nom = "epicea")

pheno <- modeliser_chapy(meteo, mode = "max")

produire_geojson_generations_chapy(pheno)
produire_geojson_envol_chapy(pheno)

# Copies a nom stable pour lecture par l'app Shiny (le plus recent = date du nom),
# puis purge : on garde les 4 dernieres dates (R/fichiers_dates.R)
source("R/fichiers_dates.R")
copier_vers_nom_stable("^chapy_generations_[0-9]{8}\\.geojson$", "chapy_generations_latest.geojson")
copier_vers_nom_stable("^chapy_onset_[0-9]{8}\\.geojson$",        "chapy_onset_latest.geojson")
purger_fichiers_dates_modele("chapy", garder = 4L)

cat("\n>> CHAPY termine.\n\n")

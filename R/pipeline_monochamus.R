# =============================================================================
# PIPELINE MONOCHAMUS — reutilise le jeu de donnees SAFRAN partage
# -----------------------------------------------------------------------------
# Le fichier [1] (mono_stade) porte le stade gen1 et, en plus, le bivoltisme [A]
# (risque_bivoltA 0/1 + dj10_projete), extrait au CENTROIDE de chaque polygone de
# stade depuis pheno$risque_bivolt et pheno$dj_projete. ATTENTION : un polygone
# regroupe plusieurs mailles, donc une seule valeur par classe de stade (4 valeurs
# en tout), pas une valeur par maille : correction a faire (voir README, « Points
# de vigilance »). Le fichier [2] (mono_bivoltB_gen2) est le bivoltisme [B].
# Toujours 2 GeoJSON au total.
# =============================================================================

source("R/safran_common.R")
source("R/monochamus_functions.R")   # garde uniquement le modele + exports GeoJSON (patchee)

meteo_wrap <- readRDS("data/cache_meteo/meteo_national.rds")
meteo      <- lapply(meteo_wrap, terra::unwrap)

pheno <- calculer_dj(meteo)

exporter_geojson(pheno)

# Copies a nom stable (le plus recent = date du nom), puis purge : 4 dernieres dates
source("R/fichiers_dates.R")
copier_vers_nom_stable("^mono_stade_[0-9]{8}_WGS84\\.geojson$",        "mono_stade_latest.geojson")
copier_vers_nom_stable("^mono_bivoltB_gen2_[0-9]{8}_WGS84\\.geojson$", "mono_bivoltB_gen2_latest.geojson")
purger_fichiers_dates_modele("mono", garder = 4L)

cat("\n>> Monochamus (avec bivoltisme A et B) termine.\n\n")

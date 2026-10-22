# =============================================================================
# PIPELINE STENO — reutilise le jeu de donnees SAFRAN partage
# -----------------------------------------------------------------------------
# R/steno_functions.R est deja allege. Le script original produisait 3 GeoJSON
# (generations, essaimage, grille) : ici grille = FALSE pour rester a 2
# fichiers, comme les 3 autres modeles.
#
# v2 : regles de Pineau et al. (2016) — temperature de l'AIR, 506 DJ au-dessus de
# 10.9 C (moyenne (Tmax+Tmin)/2) en scenario moyen, debut d'activite = Tmax > 18 C 2 jours de suite,
# aucune maturation ni arret photoperiodique. Les generations affichees sont des
# MAXIMA THERMIQUES. Parametres : bloc STENO de R/steno_functions.R.
# =============================================================================

source("R/safran_common.R")
source("R/steno_functions.R")   # garde uniquement le modele + les exports GeoJSON

meteo_wrap <- readRDS("data/cache_meteo/meteo_national.rds")
meteo      <- lapply(meteo_wrap, terra::unwrap)

# [v3] On ne garde que les mailles SAFRAN qui croisent l'essence-hote (BD Foret IGN)
source("R/masques.R")
meteo <- appliquer_masque(meteo, charger_masque("pin", meteo$tmean), nom = "pin")

# Trois scenarios = marges d'erreur de la these (voir SCENARIOS dans
# R/steno_functions.R) : pire (developpement rapide), moyen, moins_pire (lent).
# Temperature de l'air, comme la these. Le scenario "moyen" garde les noms de
# fichiers historiques (steno_*_latest.geojson) ; les autres portent un suffixe.
#
# Scenarios de sensibilite (a lancer a part, sans ecraser les sorties ci-dessous) :
#   sous ecorce : modeliser_steno(meteo, exposition = "sunny", k_ecorce = 0.2)
#   maturation  : modeliser_steno(meteo, par = modifyList(STENO, list(DJ_MATURATION = 100)))
#   photoperiode: modeliser_steno(meteo, par = modifyList(STENO, list(DL_INIT_MIN = 12)))

for (sc in names(SCENARIOS)) {
  cat(glue::glue("\n>>> STENO — scenario {sc} : base {SCENARIOS[[sc]]$T_MIN_DEV} C, ",
                 "{SCENARIOS[[sc]]$DJ_GEN} DJ\n"))
  steno <- modeliser_steno(meteo, exposition = "shaded",
                           par = modifyList(STENO, SCENARIOS[[sc]]))
  exporter_geojson(steno, pas_essaimage = 15, grille = FALSE, scenario = sc)
  rm(steno); gc()
}

copier_vers_nom_stable <- function(pattern, nom_stable) {
  fichiers <- list.files(DIR_OUT, pattern = pattern, full.names = TRUE)
  if (length(fichiers) == 0) return(invisible(NULL))
  plus_recent <- fichiers[which.max(file.mtime(fichiers))]
  file.copy(plus_recent, file.path(DIR_OUT, nom_stable), overwrite = TRUE)
  cat(glue::glue("  Copie stable : {basename(plus_recent)} -> {nom_stable}\n"))
}
# Motifs stricts : "moyen" = fichiers sans suffixe (date AAAAMMJJ / annee AAAA)
copier_vers_nom_stable("^steno_generations_[0-9]{8}\\.geojson$", "steno_generations_latest.geojson")
copier_vers_nom_stable("^steno_essaimage_[0-9]{4}\\.geojson$",   "steno_essaimage_latest.geojson")
for (sc in setdiff(names(SCENARIOS), "moyen")) {
  copier_vers_nom_stable(paste0("^steno_generations_", sc, "_[0-9]{8}\\.geojson$"),
                         paste0("steno_generations_", sc, "_latest.geojson"))
  copier_vers_nom_stable(paste0("^steno_essaimage_", sc, "_[0-9]{4}\\.geojson$"),
                         paste0("steno_essaimage_", sc, "_latest.geojson"))
}

cat("\n>> STENO termine.\n\n")

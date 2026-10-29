# =============================================================================
# FETCH SAFRAN AMONT -- 1er septembre N-1 -> 31 decembre N-1 (PHENTHAU uniquement)
# -----------------------------------------------------------------------------
# PHENTHAUproc compte les jours de froid a partir du 1er septembre de l'annee
# precedente : le package echoue sans cette periode. Le fetch partage
# (R/fetch_safran.R) ne remonte qu'au 1er janvier N et reste INCHANGE pour les 4
# modeles existants. Ce script recupere seulement la periode manquante, avec les
# memes fonctions et la meme chaine de controle qualite (safran_common.R).
#
# Les donnees de septembre a decembre N-1 sont figees : le cache (fichiers
# safran_v2_*_<N-1>.csv dans data/cache_safran, distincts de ceux de l'annee N)
# n'est retelecharge qu'une fois par an. Sortie : data/cache_meteo/meteo_amont.rds
# (tmean, tmin, tmax ; meme grille que meteo_national.rds, avec terra::wrap).
# =============================================================================

source("R/safran_common.R")

# On pilote les fonctions de la brique par ses variables globales (periode + annee des
# fichiers de cache), comme le fait R/fetch_safran.R pour l'annee courante.
ANNEE_N         <- ANNEE
ANNEE           <- ANNEE_N - 1L
DATE_DEB        <- as.Date(glue("{ANNEE}-09-01"))
DATE_FIN_SAFRAN <- min(as.Date(glue("{ANNEE}-12-31")), Sys.Date() - 5)
cat(glue(">> SAFRAN amont pour PHENTHAU : {DATE_DEB} -> {DATE_FIN_SAFRAN}\n\n"))

recuperer_safran()

dt <- lire_cache_complet()
dt <- controler_et_reparer_couverture(dt)
meteo <- assembler_rasters_nationaux(dt)
rm(dt); gc()

comblement <- tryCatch(
  combler_trous_raster(meteo),
  error = function(e) {
    message("::warning::comblement impossible : ", conditionMessage(e))
    list(meteo = meteo, n = NA_integer_)
  }
)
meteo <- comblement$meteo
verifier_couverture_finale(meteo)

# PHENTHAU n'utilise que les 3 temperatures
meteo <- meteo[c("tmean", "tmin", "tmax")]
dates_serie <- as.Date(time(meteo$tmean))

dir.create("data/cache_meteo", showWarnings = FALSE, recursive = TRUE)
saveRDS(lapply(meteo, terra::wrap), "data/cache_meteo/meteo_amont.rds")

cat(glue(
  "\n>> Meteo amont sauvegardee : data/cache_meteo/meteo_amont.rds\n",
  "   ({length(dates_serie)} jours, {min(dates_serie)} -> {max(dates_serie)}, ",
  "{comblement$n} mailles comblees)\n\n"
))

# =============================================================================
# FETCH SAFRAN — etape UNIQUE et PARTAGEE, executee une fois par run
# -----------------------------------------------------------------------------
# Telecharge (ou reutilise le cache), controle la couverture, assemble les
# rasters nationaux et calcule la daylength, puis sauvegarde le tout dans un
# seul fichier que les 4 scripts modeles (phenips, chapy, steno, monochamus)
# reutilisent directement, sans refaire ni requetes API ni reassemblage.
#
# Chaine de controle qualite (voir safran_common.R) :
#   1. recuperer_safran()                  dalles en erreur -> echec (strict)
#   2. controler_et_reparer_couverture()   mailles manquantes -> bandes etroites
#   3. combler_trous_raster()              <= 2 mailles de large -> interpolation
#   4. verifier_couverture_finale()        trous restants -> echec (strict)
#
# NOTE TECHNIQUE : un SpatRaster ne se serialise pas avec saveRDS() tel quel
# (le pointeur C++ ne survit pas a la session R). On utilise donc
# terra::wrap()/unwrap() -- obligatoire, sinon le fichier .rds est illisible
# dans le script suivant.
# =============================================================================

source("R/safran_common.R")

recuperer_safran()

dt <- lire_cache_complet()
dt <- controler_et_reparer_couverture(dt)
meteo <- assembler_rasters_nationaux(dt)
rm(dt); gc()

# Comblement des petits trous residuels. Si cette etape plante, on ne perd
# rien : verifier_couverture_finale() decide ensuite sur les donnees brutes.
comblement <- tryCatch(
  combler_trous_raster(meteo),
  error = function(e) {
    message("::warning::comblement impossible : ", conditionMessage(e))
    list(meteo = meteo, n = NA_integer_)
  }
)
meteo <- comblement$meteo
verifier_couverture_finale(meteo)

dates_serie     <- as.Date(time(meteo$tmean))
meteo$daylength <- calculer_daylength(RST_FRANCE, dates_serie)

meteo_wrap <- lapply(meteo, terra::wrap)

dir.create("data/cache_meteo", showWarnings = FALSE, recursive = TRUE)
saveRDS(meteo_wrap, "data/cache_meteo/meteo_national.rds")

# Trace du controle qualite, relue par les pipelines modeles (metadonnees)
saveRDS(
  list(date_fin_serie     = max(dates_serie),
       n_jours            = length(dates_serie),
       n_cellules_comblees = comblement$n,
       dalle_km           = DALLE_KM,
       marge_km           = MARGE_DALLE / 1000),
  "data/cache_meteo/qc_safran.rds"
)

cat(glue::glue(
  "\n>> Meteo nationale partagee sauvegardee : data/cache_meteo/meteo_national.rds\n",
  "   ({length(dates_serie)} jours, {ncell(RST_FRANCE)} cellules, ",
  "variables : {paste(names(meteo), collapse=', ')})\n\n"
))

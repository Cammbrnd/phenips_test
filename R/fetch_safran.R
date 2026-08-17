# =============================================================================
# FETCH SAFRAN — etape UNIQUE et PARTAGEE, executee une fois par run
# -----------------------------------------------------------------------------
# Telecharge (ou reutilise le cache), assemble les rasters nationaux et calcule
# la daylength, puis sauvegarde le tout dans un seul fichier que les 4 scripts
# modeles (phenips, chapy, steno, monochamus) reutilisent directement, sans
# refaire ni requetes API ni reassemblage.
#
# NOTE TECHNIQUE : un SpatRaster ne se serialise pas avec saveRDS() tel quel
# (le pointeur C++ ne survit pas a la session R). On utilise donc
# terra::wrap()/unwrap() -- obligatoire, sinon le fichier .rds est illisible
# dans le script suivant.
# =============================================================================

source("R/safran_common.R")

recuperer_safran()

dt    <- lire_cache_complet()
meteo <- assembler_rasters_nationaux(dt)
rm(dt); gc()

dates_serie     <- as.Date(time(meteo$tmean))
meteo$daylength <- calculer_daylength(RST_FRANCE, dates_serie)

meteo_wrap <- lapply(meteo, terra::wrap)

dir.create("data/cache_meteo", showWarnings = FALSE, recursive = TRUE)
saveRDS(meteo_wrap, "data/cache_meteo/meteo_national.rds")

cat(glue::glue(
  "\n>> Meteo nationale partagee sauvegardee : data/cache_meteo/meteo_national.rds\n",
  "   ({length(dates_serie)} jours, {ncell(RST_FRANCE)} cellules, ",
  "variables : {paste(names(meteo), collapse=', ')})\n\n"
))

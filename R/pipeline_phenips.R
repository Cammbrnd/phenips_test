# =============================================================================
# PIPELINE PHENIPS — reutilise le jeu de donnees SAFRAN partage
# -----------------------------------------------------------------------------
# Pour CHAQUE scenario d'exposition (sunny, shaded) :
#   - classes de generations   phenips_generations_<scenario>_<date>.geojson
#   - premier essaimage        phenips_onset_<scenario>_<date>.geojson
#   - indice continu           phenips_indice_<scenario>_<date>.geojson
# puis copie vers les noms stables *_latest.geojson lus par l'app, ecriture de
# phenips_meta.json, et purge des anciens fichiers dates.
# Les fonctions sont dans R/phenips_functions.R, le telechargement et le
# controle qualite SAFRAN dans R/safran_common.R et R/fetch_safran.R.
# =============================================================================

source("R/safran_common.R")
source("R/phenips_functions.R")

meteo_wrap <- readRDS("data/cache_meteo/meteo_national.rds")
meteo      <- lapply(meteo_wrap, terra::unwrap)

date_donnees <- max(as.Date(time(meteo$tmean)))
cat(glue(">> PHENIPS : donnees SAFRAN jusqu'au {date_donnees}\n\n"))

for (v in names(VARIANTES_PHENIPS)) {
  cat(glue("\n######## Scenario {v} : {VARIANTES_PHENIPS[[v]]$libelle} ########\n\n"))
  param <- VARIANTES_PHENIPS[[v]]
  pheno <- modeliser_phenips(meteo, exposure = param$exposure, scenario = param$scenario)

  ch_gen <- produire_geojson_generations(pheno, variante = v)
  ch_ons <- produire_geojson_onset(pheno,       variante = v)
  ch_idx <- produire_geojson_indice(pheno,      variante = v)

  # [P2] on publie le fichier qu'on vient d'ecrire, ou on retire le "latest"
  publier_latest(ch_gen, glue("phenips_generations_{v}_latest.geojson"))
  publier_latest(ch_ons, glue("phenips_onset_{v}_latest.geojson"))
  publier_latest(ch_idx, glue("phenips_indice_{v}_latest.geojson"))

  # Transition : l'ancienne app lit encore les noms sans suffixe (= sunny)
  if (COMPAT_ANCIENS_NOMS && v == "sunny") {
    publier_latest(ch_gen, "phenips_generations_latest.geojson")
    publier_latest(ch_ons, "phenips_onset_latest.geojson")
  }

  # Diagnostic de log : il ne doit jamais faire echouer la publication
  tryCatch(synthese_console(pheno),
           error = function(e) message("  [INFO] synthese_console ignoree : ", conditionMessage(e)))

  rm(pheno); invisible(gc())
}

ecrire_meta_phenips(date_donnees, names(VARIANTES_PHENIPS))
purger_fichiers_dates(garder = 4L)

cat("\n>> PHENIPS termine.\n\n")

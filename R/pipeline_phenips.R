# =============================================================================
# PIPELINE PHENIPS — reutilise le jeu de donnees SAFRAN partage
# -----------------------------------------------------------------------------
# R/phenips_functions.R est deja allege : ne contient que modeliser_phenips(),
# produire_geojson_generations(), produire_geojson_onset() et les helpers
# barrks (get_btmean, get_gen1, get_gen_rst, inspect_barrks). Les fonctions
# de telechargement/cache/assemblage sont dans R/safran_common.R.
# =============================================================================

source("R/safran_common.R")
source("R/phenips_functions.R")

meteo_wrap <- readRDS("data/cache_meteo/meteo_national.rds")
meteo      <- lapply(meteo_wrap, terra::unwrap)

pheno <- modeliser_phenips(meteo, exposure = "sunny", scenario = "max")

produire_geojson_generations(pheno)
produire_geojson_onset(pheno)
synthese_console(pheno)

copier_vers_nom_stable <- function(pattern, nom_stable) {
  fichiers <- list.files(DIR_OUT, pattern = pattern, full.names = TRUE)
  if (length(fichiers) == 0) return(invisible(NULL))
  plus_recent <- fichiers[which.max(file.mtime(fichiers))]
  file.copy(plus_recent, file.path(DIR_OUT, nom_stable), overwrite = TRUE)
  cat(glue::glue("  Copie stable : {basename(plus_recent)} -> {nom_stable}\n"))
}
copier_vers_nom_stable("^phenips_generations_.*\\.geojson$", "phenips_generations_latest.geojson")
copier_vers_nom_stable("^phenips_onset_.*\\.geojson$",       "phenips_onset_latest.geojson")

cat("\n>> PHENIPS termine.\n\n")

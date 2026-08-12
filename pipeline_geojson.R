# =============================================================================
# PIPELINE ALLEGE — GeoJSON generations + onset uniquement
# -----------------------------------------------------------------------------
# Utilise pour l'automatisation GitHub Actions.
# Reprend les fonctions de phenips_functions.R (script original) mais saute
# les exports GeoTIFF et les cartes PNG, inutiles en production et plus lents.
# =============================================================================

source("R/phenips_functions.R")

pipeline_geojson_only <- function(exposure = "sunny", scenario = "max") {
  cat("\n", strrep("=", 60), "\n")
  cat("  PHENIPS-Clim — Pipeline GeoJSON (CI)\n")
  cat(strrep("=", 60), "\n\n")
  t0 <- proc.time()

  # 1. Telechargement / mise a jour cache SAFRAN (rapide si deja a jour)
  recuperer_safran()

  # 2. Lecture et assemblage
  dt    <- lire_cache_complet()
  meteo <- assembler_rasters_nationaux(dt)
  rm(dt); gc()

  # 3. Daylength
  dates_serie     <- as.Date(time(meteo$tmean))
  meteo$daylength <- calculer_daylength(RST_FRANCE, dates_serie)

  # 4. Modelisation PHENIPS-Clim
  pheno <- modeliser_phenips(meteo, exposure = exposure, scenario = scenario)
  rm(meteo); gc()

  # 5. Export GeoJSON uniquement (pas de GeoTIFF, pas de PNG)
  produire_geojson_generations(pheno)
  produire_geojson_onset(pheno)

  # 6. Copies a nom fixe pour lecture stable par l'app Shiny
  #    (les fonctions ci-dessus nomment les fichiers avec la date du jour,
  #     l'app a besoin d'une URL qui ne change jamais)
  copier_vers_nom_stable <- function(pattern, nom_stable) {
    fichiers <- list.files(DIR_OUT, pattern = pattern, full.names = TRUE)
    if (length(fichiers) == 0) {
      cat(glue("  [INFO] Aucun fichier correspondant a {pattern} — pas de copie stable.\n"))
      return(invisible(NULL))
    }
    plus_recent <- fichiers[which.max(file.mtime(fichiers))]
    file.copy(plus_recent, file.path(DIR_OUT, nom_stable), overwrite = TRUE)
    cat(glue("  Copie stable : {basename(plus_recent)} -> {nom_stable}\n"))
  }

  copier_vers_nom_stable("^phenips_generations_.*\\.geojson$", "phenips_generations_latest.geojson")
  copier_vers_nom_stable("^phenips_onset_.*\\.geojson$",       "phenips_onset_latest.geojson")

  duree <- (proc.time() - t0)[["elapsed"]]
  cat(glue("\n  Duree totale : {round(duree/60, 1)} min\n"))
  cat(glue("  Sorties      : {DIR_OUT}/\n\n"))

  invisible(pheno)
}

# --- Lancement ---
resultats <- pipeline_geojson_only(exposure = "sunny", scenario = "max")

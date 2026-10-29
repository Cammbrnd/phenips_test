# =============================================================================
# ASSEMBLAGE DES TABLES HISTORIQUES (base R uniquement)
# -----------------------------------------------------------------------------
# Lit data/historique_tmp/<modele>_<annee>.csv (un par annee, issus des artefacts du
# workflow), les fusionne avec les tables existantes dans data/historique/ (les annees
# recalculees remplacent les anciennes) et ecrit :
#   data/historique/phenips_annuel.csv.gz   annee,x,y,gen_sunny,doy_sunny,gen_shaded,doy_shaded
#   data/historique/chapy_annuel.csv.gz     annee,x,y,gen,doy_envol
#   data/historique/steno_annuel.csv.gz     annee,x,y,gen_pire,doy_pire,gen_moyen,doy_moyen,gen_moins_pire,doy_moins_pire
#   data/historique/meta.json
# Coordonnees x,y : centre de maille SAFRAN 8 km, Lambert-93 (EPSG:2154).
# =============================================================================

DIR_TMP <- "data/historique_tmp"
DIR_OUT <- "data/historique"
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

modeles <- c("phenips", "chapy", "steno")
meta <- list()
for (m in modeles) {
  fichiers <- list.files(DIR_TMP, pattern = sprintf("^%s_[0-9]{4}\\.csv$", m), full.names = TRUE)
  cible <- file.path(DIR_OUT, sprintf("%s_annuel.csv.gz", m))
  ancien <- if (file.exists(cible)) utils::read.csv(gzfile(cible)) else NULL
  if (length(fichiers) == 0 && is.null(ancien)) next
  nouveau <- if (length(fichiers)) do.call(rbind, lapply(fichiers, utils::read.csv)) else NULL
  if (!is.null(nouveau) && !is.null(ancien)) {
    ancien <- ancien[!ancien$annee %in% unique(nouveau$annee), , drop = FALSE]
    nouveau <- nouveau[, names(ancien), drop = FALSE]
  }
  tout <- rbind(ancien, nouveau)
  tout <- tout[order(tout$annee, tout$y, tout$x), ]
  con <- gzfile(cible, "w"); utils::write.csv(tout, con, row.names = FALSE, na = ""); close(con)
  cat(sprintf("%s : %d lignes, annees %d-%d (%d annees), %d mailles -> %s (%.1f Mo)\n", m, nrow(tout),
              min(tout$annee), max(tout$annee), length(unique(tout$annee)),
              nrow(unique(tout[, c("x", "y")])), cible, file.size(cible) / 1e6))
  meta[[m]] <- sprintf('"%s": {"annee_min": %d, "annee_max": %d, "n_annees": %d, "n_mailles": %d, "n_lignes": %d}',
                       m, min(tout$annee), max(tout$annee), length(unique(tout$annee)),
                       nrow(unique(tout[, c("x", "y")])), nrow(tout))
}
if (length(meta)) {
  writeLines(c("{",
    sprintf('  "date_construction": "%s",', format(Sys.Date())),
    '  "source": "Meteo-France (SAFRAN-ISBA) via GeoSAS INRAE ; masques IGN BD Foret V2",',
    '  "crs_xy": "EPSG:2154 (centre de maille 8 km)",',
    paste0("  ", paste(unlist(meta), collapse = ",\n  ")),
    "}"), file.path(DIR_OUT, "meta.json"))
}

# =============================================================================
# FICHIERS DATES : copie vers le nom stable (*_latest) et purge
# -----------------------------------------------------------------------------
# Utilise par CHAPY, STENO et Monochamus (PHENIPS et PHENTHAU ont leurs propres
# fonctions, meme logique). Base R uniquement, sans dependance.
#
#   - Le fichier "le plus recent" se lit dans la DATE DU NOM (AAAAMMJJ, ou AAAA
#     pour un fichier annuel), jamais dans la date de modification : un checkout
#     Git remet tous les fichiers a la meme date.
#   - La purge garde les `garder` dates les plus recentes par type de fichier
#     (les anciens fichiers restent consultables dans l'historique Git).
# =============================================================================

# Cle (chiffres) lue dans le nom : AAAAMMJJ ou AAAA, avec suffixe _WGS84 optionnel
cle_du_nom <- function(noms) {
  motif <- "^.*_([0-9]{4}|[0-9]{8})(_WGS84)?\\.geojson$"
  ifelse(grepl(motif, noms), sub(motif, "\\1", noms), NA_character_)
}

#' Copie le fichier le plus recent (selon la date du nom) vers `nom_stable`.
#' @param pattern expression reguliere (sur le nom de fichier) des fichiers dates
copier_vers_nom_stable <- function(pattern, nom_stable, dir_out = "data") {
  f   <- list.files(dir_out, pattern = pattern)
  cle <- cle_du_nom(f)
  f   <- f[!is.na(cle)]; cle <- cle[!is.na(cle)]
  if (length(f) == 0) {
    message("::warning::Aucun fichier date ne correspond a ", pattern, " : ", nom_stable, " inchange")
    return(invisible(NULL))
  }
  plus_recent <- f[which.max(cle)]
  file.copy(file.path(dir_out, plus_recent), file.path(dir_out, nom_stable), overwrite = TRUE)
  cat("  Copie stable : ", plus_recent, " -> ", nom_stable, "\n", sep = "")
  invisible(plus_recent)
}

#' Supprime les anciens fichiers dates d'un modele : on garde les `garder` dates
#' les plus recentes pour chaque type de fichier (ex. steno_generations,
#' steno_generations_pire, ...). Ne touche ni aux *_latest ni aux fichiers annuels.
#' @param prefixe debut du nom de fichier du modele ("chapy", "steno", "mono")
purger_fichiers_dates_modele <- function(prefixe, garder = 4L, dir_out = "data") {
  f <- list.files(dir_out, pattern = paste0("^", prefixe, "_.*_[0-9]{8}(_WGS84)?\\.geojson$"))
  if (length(f) == 0) return(invisible(character(0)))
  date   <- sub("^.*_([0-9]{8})(_WGS84)?\\.geojson$", "\\1", f)
  groupe <- sub("_[0-9]{8}(_WGS84)?\\.geojson$", "", f)
  a_supprimer <- character(0)
  for (g in unique(groupe)) {
    dates_g <- sort(unique(date[groupe == g]), decreasing = TRUE)
    if (length(dates_g) > garder)
      a_supprimer <- c(a_supprimer, f[groupe == g & date %in% dates_g[-seq_len(garder)]])
  }
  if (length(a_supprimer) > 0) {
    file.remove(file.path(dir_out, a_supprimer))
    cat("  Purge ", prefixe, " : ", length(a_supprimer), " ancien(s) fichier(s) date(s) supprime(s).\n", sep = "")
  }
  invisible(a_supprimer)
}

# =============================================================================
# MASQUES D'ESSENCE (v3) : ne garder que les mailles SAFRAN qui croisent l'essence-hote
# -----------------------------------------------------------------------------
#   epicea -> PHENIPS (Ips typographus) et CHAPY (Pityogenes chalcographus)
#   pin    -> STENO (Ips sexdentatus)
#   chene  -> PHENTHAU (Thaumetopoea processionea ; "Chenes decidus")
# Les fichiers data/masques/masque_<nom>.csv sont produits par R/construire_masques.R
# (BD Foret V2, IGN). Une maille masquee passe a NA dans tous les rasters
# meteo : le modele ne la calcule pas et elle n'apparait pas dans les GeoJSON.
#
# SEUIL_PART_MASQUE : part minimale de la maille (0-1) couverte par l'essence pour
# la garder. 0 = la maille est gardee des qu'elle croise un peuplement de l'essence.
# =============================================================================

SEUIL_PART_MASQUE <- 0
DIR_MASQUES       <- "data/masques"

charger_masque <- function(nom, modele_raster, seuil = SEUIL_PART_MASQUE) {
  f <- file.path(DIR_MASQUES, paste0("masque_", nom, ".csv"))
  if (!file.exists(f))
    stop("Masque introuvable : ", f, ". Lancer le workflow avec l'entree 'construire_masques'.", call. = FALSE)
  d <- utils::read.csv(f)
  d <- d[d$part >= seuil, , drop = FALSE]
  m <- terra::rast(modele_raster[[1]])          # meme grille que la meteo, 1 couche vide
  terra::values(m) <- NA_real_
  cel <- terra::cellFromXY(m, cbind(d$x, d$y))
  if (anyNA(cel) || length(cel) == 0)
    stop("Masque '", nom, "' incompatible avec la grille meteo (mailles hors emprise ou masque vide)", call. = FALSE)
  m[cel] <- 1
  m
}

appliquer_masque <- function(meteo, masque, nom = "essence") {
  avant <- sum(!is.na(terra::values(meteo$tmean[[terra::nlyr(meteo$tmean)]])))
  out <- lapply(meteo, function(r) {
    tps <- terra::time(r)
    m   <- terra::mask(r, masque)
    if (!all(is.na(tps))) terra::time(m) <- tps
    m
  })
  names(out) <- names(meteo)
  apres <- sum(!is.na(terra::values(out$tmean[[terra::nlyr(out$tmean)]])))
  cat(sprintf(">> Masque %s : %d mailles gardees sur %d avec donnees SAFRAN (%.0f %%)\n",
              nom, apres, avant, 100 * apres / max(avant, 1)))
  if (apres < 100)
    stop("Masque '", nom, "' : moins de 100 mailles avec donnees apres masquage", call. = FALSE)
  out
}

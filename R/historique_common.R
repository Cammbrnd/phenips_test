# =============================================================================
# HISTORIQUE 1960+ : fonctions communes aux scripts R/historique_<modele>.R
# -----------------------------------------------------------------------------
# Un script = UNE annee (variable d'environnement ANNEE_SAFRAN), UN modele, tous ses
# scenarios. Sortie : data/historique_tmp/<modele>_<annee>.csv, une ligne par maille
# du masque d'essence (v3) avec donnees SAFRAN.
#   gen_<scenario> : nombre de generations au 31/12 (ou a la derniere date si annee en cours)
#   doy_<scenario> : jour de l'annee du premier envol / essaimage (vide si aucun)
# A sourcer APRES safran_common.R et le fichier de fonctions du modele (ces derniers
# redefinissent ANNEE, d'ou la reaffectation ci-dessous).
# =============================================================================

ANNEE_H <- suppressWarnings(as.integer(Sys.getenv("ANNEE_SAFRAN", "")))
if (is.na(ANNEE_H)) stop("ANNEE_SAFRAN non definie (ex. ANNEE_SAFRAN=1960)", call. = FALSE)
ANNEE <- ANNEE_H
DIR_HISTO_TMP <- "data/historique_tmp"

source("R/masques.R")

#' Meteo nationale de l'annee, masquee sur l'essence-hote ("epicea" ou "pin").
charger_meteo_masquee <- function(nom_masque) {
  meteo <- lapply(readRDS("data/cache_meteo/meteo_national.rds"), terra::unwrap)
  appliquer_masque(meteo, charger_masque(nom_masque, meteo$tmean), nom = nom_masque)
}

#' Mailles gardees : masque d'essence ET donnees SAFRAN a la derniere date.
cellules_valides <- function(meteo) {
  v  <- terra::values(meteo$tmean[[terra::nlyr(meteo$tmean)]])[, 1]
  id <- which(!is.na(v))
  xy <- terra::xyFromCell(meteo$tmean, id)
  data.frame(id = id, x = round(xy[, 1]), y = round(xy[, 2]))
}

#' Libelles ("-", "1", "1s", "2"...) d'un raster categoriel barrks, un par cellule.
etiquettes_raster <- function(rst) {
  df <- tryCatch(terra::values(rst, dataframe = TRUE), error = function(e) NULL)
  x  <- if (!is.null(df)) df[[1]] else terra::values(rst)[, 1]
  if (is.factor(x)) return(as.character(x))
  lv <- tryCatch(terra::cats(rst)[[1]], error = function(e) NULL)
  if (!is.null(lv) && ncol(lv) >= 2) return(as.character(lv[[2]][match(x, lv[[1]])]))
  as.character(x)
}

#' Jour de l'annee du premier jour ou une pile "onset" passe au-dessus de 0 (NA si jamais).
doy_premier <- function(onset, dates) {
  n <- terra::ncell(onset)
  if (is.null(onset) || terra::nlyr(onset) == 0) return(rep(NA_real_, n))
  vals <- terra::values(onset)
  a    <- rowSums(vals > 0, na.rm = TRUE) > 0
  idx  <- rep(NA_integer_, nrow(vals))
  if (any(a)) {
    sb <- vals[a, , drop = FALSE] > 0
    sb[is.na(sb)] <- FALSE
    idx[a] <- max.col(sb, ties.method = "first")
  }
  ok  <- !is.na(idx)
  doy <- rep(NA_real_, length(idx))
  doy[ok] <- as.integer(format(dates[idx[ok]], "%j"))
  doy
}

#' Resume dans les annotations GitHub (les journaux ne sont pas lisibles a distance).
resume_colonne <- function(etiquette, x) {
  t <- table(x, useNA = "ifany")
  message("::notice::", etiquette, " : ", paste(names(t), as.integer(t), sep = "=", collapse = " | "))
}

ecrire_table_annee <- function(df, nom) {
  dir.create(DIR_HISTO_TMP, showWarnings = FALSE, recursive = TRUE)
  f <- file.path(DIR_HISTO_TMP, sprintf("%s_%d.csv", nom, ANNEE_H))
  utils::write.csv(df, f, row.names = FALSE, na = "")
  cat(sprintf(">> %s : %d mailles, %d colonnes -> %s\n", nom, nrow(df), ncol(df), f))
  invisible(f)
}

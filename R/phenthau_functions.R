# =============================================================================
# PHENTHAU -- Thaumetopoea processionea (processionnaire du chene)
# -----------------------------------------------------------------------------
# Modele  : package PHENTHAUproc (Halbig et al. 2024, Forest Ecol. Manag. 552:121525)
#           mode regional "dailymeanminmax" : tmean / tmin / tmax JOURNALIERES SAFRAN
#           (T_Q, TINF_H_Q, TSUP_H_Q) : aucune reconstruction horaire.
# Entree  : serie continue du 1er septembre N-1 au plus tard au 30 septembre N
#           (le package compte les jours de froid depuis le 1er septembre N-1).
# Sorties du package : stades 0-8 par jour (1 fev -> 30 sept), eclosion des L1 selon
#           3 modeles (Custers, Wagenhoff, Meurisse), debourrement du chene, mortalite L1.
#
# LIMITE D'APPLICABILITE : le modele est calibre sur Quercus robur en Allemagne, jusqu'a
# 800 m d'altitude. Ailleurs (sud de la France, montagne) les sorties sont une
# extrapolation : AVERTISSEMENT_PHENTHAU est ecrit dans chaque GeoJSON et le meta.
#
# Le fetch SAFRAN partage ne remonte qu'au 1er janvier N : la periode septembre-decembre
# N-1 vient de R/fetch_safran_amont.R (meteo_amont.rds).
# =============================================================================

library(PHENTHAUproc)
library(sf)

# Repertoire de sortie : "data" en production ; les tests sur petite emprise le redirigent
# (variable d'environnement PHENTHAU_DIR_OUT) pour ne rien publier par erreur.
DIR_OUT <- Sys.getenv("PHENTHAU_DIR_OUT", "data")
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

ESPECE_PHENTHAU <- "Thaumetopoea processionea"
MODELE_PHENTHAU <- "PHENTHAU"

# Definition de la date d'eclosion des L1 utilisee pour la suite du modele :
# "first" = la plus precoce des 3 modeles (defaut du package) ; "mean" = jour moyen.
DEF_ECLOSION <- Sys.getenv("PHENTHAU_DEF_ECLOSION", "first")   # variable d'environnement : tests de comparaison seulement
if (!DEF_ECLOSION %in% c("first", "mean")) stop("PHENTHAU_DEF_ECLOSION doit valoir first ou mean", call. = FALSE)

SOURCE_PHENTHAU <- paste0("Source : Meteo-France (donnees SAFRAN) via GeoSAS INRAE ; ",
                          "modele PHENTHAUproc (Halbig et al. 2024)")
AVERTISSEMENT_PHENTHAU <- paste0("Modele calibre sur Quercus robur en Allemagne jusqu'a 800 m : ",
                                 "hors de ce domaine (sud de la France, montagne) ",
                                 "les resultats sont une extrapolation")

# Stades du package (get_legend("stages")) : couleurs figees ici en hexadecimal.
# L3 = debut du risque urticant (poils urticants a partir du 3e stade larvaire).
TABLE_STADES <- data.frame(
  stade   = 0:8,
  code    = c("Oeuf", "L1", "L2", "L3", "L4", "L5", "L6", "Pp", "Ad"),
  label   = c("Oeuf (pas d'eclosion)", "L1", "L2", "L3 (debut risque urticant)",
              "L4", "L5", "L6", "Pupe", "Adulte"),
  couleur = c("#E8E8E8", "#BEBEBE", "#696969", "#FFFFE0", "#FFFF00",
              "#FFA500", "#FF0000", "#D3D3D3", "#A9A9A9"),
  stringsAsFactors = FALSE
)
STADE_L3 <- 3L

# Mortalite des L1 (famine entre eclosion et debourrement, %), 5 classes + "a determiner"
CLASSES_MORTALITE <- data.frame(
  ordre   = 0:5,
  borne   = c(NA, 0, 10, 30, 60, 90),
  label   = c("A determiner (pas encore de debourrement)", "< 10 %", "10 - 30 %",
              "30 - 60 %", "60 - 90 %", ">= 90 %"),
  couleur = c("#E8E8E8", "#FFFFCC", "#FED976", "#FD8D3C", "#E31A1C", "#800026"),
  stringsAsFactors = FALSE
)


# =============================================================================
# 1. SERIE METEO : amont (sept-dec N-1) + annee N, controles
# =============================================================================

#' Assemble tmean / tmin / tmax en une serie continue du 1er sept N-1 jusqu'a
#' min(derniere date disponible, 30 sept N). `meteo_amont` peut etre NULL si `meteo`
#' couvre deja tout (tests sur petite emprise). Echec explicite si la serie a un trou.
assembler_serie_phenthau <- function(meteo, meteo_amont = NULL, annee = ANNEE) {
  debut   <- as.Date(sprintf("%d-09-01", annee - 1L))
  fin_max <- as.Date(sprintf("%d-09-30", annee))
  vars    <- c("tmean", "tmin", "tmax")
  serie <- lapply(vars, function(v) {
    r <- if (is.null(meteo_amont)) meteo[[v]] else c(meteo_amont[[v]], meteo[[v]])
    d <- if (is.null(meteo_amont)) as.Date(terra::time(meteo[[v]]))
         else c(as.Date(terra::time(meteo_amont[[v]])), as.Date(terra::time(meteo[[v]])))
    o <- order(d); r <- r[[o]]; d <- d[o]
    k <- !duplicated(d) & d >= debut & d <= fin_max
    r <- r[[which(k)]]; terra::time(r) <- d[k]
    r
  })
  names(serie) <- vars
  d <- as.Date(terra::time(serie$tmean))
  if (length(d) == 0 || min(d) != debut)
    stop("PHENTHAU : la serie SAFRAN ne commence pas le ", debut,
         " (meteo_amont.rds absent ou incomplet ?)", call. = FALSE)
  manquants <- setdiff(as.character(seq(debut, max(d), by = "day")), as.character(d))
  if (length(manquants) > 0)
    stop("PHENTHAU : ", length(manquants), " jour(s) manquant(s) dans la serie SAFRAN, ex. ",
         paste(head(manquants, 3), collapse = ", "), call. = FALSE)
  serie
}


# =============================================================================
# 2. MODELE
# =============================================================================

modeliser_phenthau <- function(serie, annee = ANNEE, def_hatch = DEF_ECLOSION) {
  res <- suppressMessages(
    phenthau(serie[c("tmean", "tmin", "tmax")],
             params = parameter("dailymeanminmax", year = annee), def_hatch = def_hatch)
  )
  dates <- as.Date(terra::time(res$stages))
  list(res = res, dates = dates, annee = annee, date_ref = max(dates))
}

#' DOY (jour de l'annee) du premier jour ou `cond` est vraie, par maille (matrice mailles x jours)
doy_premier <- function(m, dates, cond = function(x) x > 0) {
  ok <- cond(m)
  ok[is.na(ok)] <- FALSE
  idx <- max.col(ok, ties.method = "first")
  doy <- as.integer(format(dates[idx], "%j"))
  doy[rowSums(ok) == 0L] <- NA_integer_
  doy
}

#' Tous les rasters derives (vecteurs par maille de RST_FRANCE) : un seul calcul, reutilise
#' par le GeoTIFF, les GeoJSON et la synthese.
calculer_rasters_phenthau <- function(pheno) {
  res <- pheno$res; dates <- pheno$dates
  S <- terra::values(res$stages)
  v <- list(
    doy_debourrement   = doy_premier(terra::values(res$budswelling),   dates),
    doy_feuillaison    = doy_premier(terra::values(res$leafunfolding), dates),
    doy_eclosion_custers   = doy_premier(terra::values(res$custers),   dates),
    doy_eclosion_wagenhoff = doy_premier(terra::values(res$wagenhoff), dates),
    doy_eclosion_meurisse  = doy_premier(terra::values(res$meurisse),  dates)
  )
  for (s in 1:8) v[[paste0("doy_stade_", s)]] <- doy_premier(S, dates, function(x) x >= s)
  v$doy_eclosion <- v$doy_stade_1                 # eclosion retenue par le modele (DEF_ECLOSION)

  # Mortalite L1 : definitive seulement une fois le debourrement atteint
  # (avant, le package renvoie 0, ce qui serait trompeur).
  mort <- as.numeric(terra::values(res$mortality)[, 1])
  mort[is.na(v$doy_debourrement)] <- NA_real_
  v$mortalite_pct <- mort
  v$stade_courant <- as.integer(S[, ncol(S)])
  v$risque_urticant <- as.integer(v$stade_courant >= STADE_L3)
  v
}

vers_raster <- function(vec) {
  r <- RST_FRANCE
  terra::values(r) <- vec
  r
}


# =============================================================================
# 3. EXPORT GEOJSON (dissous, champs couleur precalcules)
# =============================================================================

#' Vectorise sur la grille native L93 puis reprojette (voir phenips_functions.R [v3]).
vectoriser_phenthau <- function(vec, nom_champ) {
  v <- tryCatch(
    terra::as.polygons(vers_raster(vec), dissolve = TRUE) |> sf::st_as_sf() |> sf::st_transform(4326),
    error = function(e) { message("::warning::PHENTHAU vectorisation : ", conditionMessage(e)); NULL }
  )
  if (is.null(v) || nrow(v) == 0) return(NULL)
  names(v)[1] <- nom_champ
  v
}

ajouter_meta_phenthau <- function(v, date_ref, annee) {
  v$date_calcul    <- as.character(date_ref)
  v$annee          <- annee
  v$espece         <- ESPECE_PHENTHAU
  v$modele         <- MODELE_PHENTHAU
  v$def_eclosion   <- DEF_ECLOSION
  v$source         <- SOURCE_PHENTHAU
  v$avertissement  <- AVERTISSEMENT_PHENTHAU
  v
}

ecrire_geojson_phenthau <- function(v, nom_fichier, dir_out = DIR_OUT) {
  v <- sf::st_simplify(v, dTolerance = 0.005, preserveTopology = TRUE)
  v <- v[!sf::st_is_empty(v), ]
  chemin <- file.path(dir_out, nom_fichier)
  sf::st_write(v, chemin, delete_dsn = TRUE, quiet = TRUE,
               layer_options = "COORDINATE_PRECISION=5")
  cat(glue("  {nrow(v)} polygones | {round(file.size(chemin) / 1024, 0)} Ko -> {nom_fichier}\n\n"))
  chemin
}

nom_fichier_phenthau <- function(couche, date_ref) {
  glue("phenthau_{couche}_{format(date_ref, '%Y%m%d')}.geojson")
}

#' Stade atteint a la date des donnees (classes 0-8)
produire_geojson_stades <- function(pheno, ras, dir_out = DIR_OUT) {
  cat(">> [PHENTHAU] Stade a la date des donnees -> GeoJSON...\n")
  v <- vectoriser_phenthau(ras$stade_courant, "stade")
  if (is.null(v)) { cat("  Aucun polygone produit.\n\n"); return(invisible(NULL)) }
  i <- match(as.integer(v$stade), TABLE_STADES$stade)
  v$stade_code     <- TABLE_STADES$code[i]
  v$label          <- TABLE_STADES$label[i]
  v$couleur        <- TABLE_STADES$couleur[i]
  v$ordre          <- v$stade
  v$risque_urticant <- as.integer(v$stade >= STADE_L3)
  v <- v[order(v$ordre), ]
  v <- ajouter_meta_phenthau(v, pheno$date_ref, pheno$annee)
  invisible(ecrire_geojson_phenthau(v, nom_fichier_phenthau("stades", pheno$date_ref), dir_out))
}

#' Couche de dates (gradient jj/mm cote app : champ date_iso) : couche = "l3", "eclosion", "debourrement"
produire_geojson_date <- function(pheno, ras, nom_doy, couche, champ_date, dir_out = DIR_OUT) {
  cat(glue(">> [PHENTHAU] Date {couche} -> GeoJSON...\n"))
  doy <- ras[[nom_doy]]
  if (all(is.na(doy))) {
    cat("  Aucune maille n'a encore atteint ce stade.\n\n")
    return(invisible(NULL))
  }
  v <- vectoriser_phenthau(doy, "doy")
  if (is.null(v)) { cat("  Aucun polygone produit.\n\n"); return(invisible(NULL)) }
  origine <- as.Date(sprintf("%d-01-01", pheno$annee))
  d <- origine + as.integer(v$doy) - 1L
  v[[champ_date]] <- format(d, "%d/%m/%Y")
  v$date_iso      <- format(d, "%Y-%m-%d")
  v$semaine       <- as.integer(format(d, "%V"))
  v$doy <- NULL
  v <- ajouter_meta_phenthau(v, pheno$date_ref, pheno$annee)
  invisible(ecrire_geojson_phenthau(v, nom_fichier_phenthau(couche, pheno$date_ref), dir_out))
}

#' Mortalite des L1 par classes
produire_geojson_mortalite <- function(pheno, ras, dir_out = DIR_OUT) {
  cat(">> [PHENTHAU] Mortalite L1 -> GeoJSON...\n")
  m <- ras$mortalite_pct
  presentes <- !is.na(ras$stade_courant)                  # mailles du masque chenaie
  if (!any(presentes)) { cat("  Aucune maille.\n\n"); return(invisible(NULL)) }
  classe <- rep(NA_integer_, length(m))
  classe[presentes] <- 0L                                 # a determiner
  det <- presentes & !is.na(m)
  classe[det] <- findInterval(m[det], CLASSES_MORTALITE$borne[-1])   # 1..5 = ordre des classes
  v <- vectoriser_phenthau(classe, "ordre")
  if (is.null(v)) { cat("  Aucun polygone produit.\n\n"); return(invisible(NULL)) }
  v$ordre   <- as.integer(v$ordre)
  i <- match(v$ordre, CLASSES_MORTALITE$ordre)
  v$label   <- CLASSES_MORTALITE$label[i]
  v$couleur <- CLASSES_MORTALITE$couleur[i]
  v <- v[order(v$ordre), ]
  v <- ajouter_meta_phenthau(v, pheno$date_ref, pheno$annee)
  invisible(ecrire_geojson_phenthau(v, nom_fichier_phenthau("mortalite", pheno$date_ref), dir_out))
}


# =============================================================================
# 4. EXPORT GEOTIFF (pile des DOY par maille, L93, parite PHENIPS : LZW + TILED + BIGTIFF)
# =============================================================================

exporter_geotiff_phenthau <- function(pheno, ras, dir_out = DIR_OUT) {
  cat(">> [PHENTHAU] Export GeoTIFF (L93, LZW)...\n")
  couches <- lapply(names(ras), function(n) { r <- vers_raster(as.numeric(ras[[n]])); names(r) <- n; r })
  pile <- terra::rast(couches)
  nom  <- glue("phenthau_{pheno$annee}_L93.tif")
  terra::writeRaster(pile, file.path(dir_out, nom), overwrite = TRUE,
                     gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))
  cat(glue("  {nlyr(pile)} bandes -> {nom}\n\n"))
  invisible(file.path(dir_out, nom))
}


# =============================================================================
# 4b. CARTES NATIONALES PNG (stade a la date des donnees + date de debut du risque urticant L3)
# =============================================================================
# ggplot2 seul (deja installe par le workflow), grille L93 en tuiles : pas de fond de carte.

cartes_nationales_phenthau <- function(pheno, ras, dir_out = DIR_OUT) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    cat("  ggplot2 absent : PNG ignores.\n\n"); return(invisible(NULL))
  }
  cat(">> [PHENTHAU] Cartes nationales PNG...\n")
  date_j     <- format(pheno$date_ref, "%Y%m%d")
  date_label <- format(pheno$date_ref, "%d/%m/%Y")
  vers_df <- function(vec) {
    r <- vers_raster(vec)
    d <- as.data.frame(r, xy = TRUE, na.rm = TRUE)
    names(d)[3] <- "v"
    d
  }
  theme_carte <- ggplot2::theme_void(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   plot.caption = ggplot2::element_text(size = 8, colour = "grey40", hjust = 0),
                   plot.caption.position = "plot")
  legende_cadre <- paste(c(strwrap(AVERTISSEMENT_PHENTHAU, 90), strwrap(SOURCE_PHENTHAU, 90)), collapse = "\n")

  # 1. Stade a la date des donnees
  d1 <- vers_df(ras$stade_courant)
  if (nrow(d1) > 0) {
    d1$stade <- factor(TABLE_STADES$label[match(d1$v, TABLE_STADES$stade)], levels = TABLE_STADES$label)
    p1 <- ggplot2::ggplot(d1, ggplot2::aes(x, y, fill = stade)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_manual(values = stats::setNames(TABLE_STADES$couleur, TABLE_STADES$label),
                                 drop = TRUE, name = "Stade") +
      ggplot2::coord_equal() +
      ggplot2::labs(title = glue("Processionnaire du chene - stade au {date_label}"),
                    subtitle = "Mailles SAFRAN 8 km croisant des chenes decidus (BD Foret V2)",
                    caption = legende_cadre) + theme_carte
    f1 <- file.path(dir_out, glue("phenthau_stades_{date_j}.png"))
    ggplot2::ggsave(f1, p1, width = 8, height = 8, dpi = 150, bg = "white")
    cat(glue("  -> {basename(f1)}\n"))
  }

  # 2. Date de debut du risque urticant (L3)
  d2 <- vers_df(ras$doy_stade_3)
  if (nrow(d2) > 0) {
    d2$date <- as.Date(sprintf("%d-01-01", pheno$annee)) + d2$v - 1L
    p2 <- ggplot2::ggplot(d2, ggplot2::aes(x, y, fill = date)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradientn(colours = c("#FFF3B0", "#FDE725", "#35B779", "#31688E", "#440154"),
                                    trans = "date", name = "Date", labels = function(x) format(x, "%d/%m")) +
      ggplot2::coord_equal() +
      ggplot2::labs(title = glue("Processionnaire du chene - debut du risque urticant (L3) {pheno$annee}"),
                    subtitle = glue("Donnees au {date_label}"), caption = legende_cadre) + theme_carte
    f2 <- file.path(dir_out, glue("phenthau_l3_{date_j}.png"))
    ggplot2::ggsave(f2, p2, width = 8, height = 8, dpi = 150, bg = "white")
    cat(glue("  -> {basename(f2)}\n"))
  }
  cat("\n")
  invisible(TRUE)
}


# =============================================================================
# 5. PUBLICATION : nom stable, metadonnees, purge (memes regles que phenips_functions.R)
# =============================================================================

publier_latest_phenthau <- function(chemin, nom_stable, dir_out = DIR_OUT) {
  cible <- file.path(dir_out, nom_stable)
  if (is.null(chemin)) {
    if (file.exists(cible)) {
      file.remove(cible)
      cat(glue("  Aucune donnee : {nom_stable} supprime (pas de carte perimee).\n"))
    }
    return(invisible(FALSE))
  }
  file.copy(chemin, cible, overwrite = TRUE)
  cat(glue("  Copie stable : {basename(chemin)} -> {nom_stable}\n"))
  invisible(TRUE)
}

purger_fichiers_phenthau <- function(dir_out = DIR_OUT, garder = 4L) {
  f <- list.files(dir_out, pattern = "^phenthau_.*_[0-9]{8}\\.geojson$")
  if (length(f) == 0) return(invisible(character(0)))
  date    <- sub("^.*_([0-9]{8})\\.geojson$", "\\1", f)
  prefixe <- sub("_[0-9]{8}\\.geojson$", "", f)
  a_supprimer <- character(0)
  for (p in unique(prefixe)) {
    dates_p <- sort(unique(date[prefixe == p]), decreasing = TRUE)
    if (length(dates_p) > garder)
      a_supprimer <- c(a_supprimer, f[prefixe == p & date %in% dates_p[-seq_len(garder)]])
  }
  if (length(a_supprimer) > 0) {
    file.remove(file.path(dir_out, a_supprimer))
    cat(glue("  Purge : {length(a_supprimer)} ancien(s) fichier(s) date(s) supprime(s).\n"))
  }
  invisible(a_supprimer)
}

ecrire_meta_phenthau <- function(pheno, n_mailles, dir_out = DIR_OUT) {
  json_str <- function(x) paste0('"', gsub('"', '\\\\"', as.character(x)), '"')
  champs <- c(
    sprintf('"date_donnees": %s',    json_str(pheno$date_ref)),
    sprintf('"date_generation": %s', json_str(format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))),
    sprintf('"annee": %d',           pheno$annee),
    sprintf('"espece": %s',          json_str(ESPECE_PHENTHAU)),
    sprintf('"modele": %s',          json_str(MODELE_PHENTHAU)),
    sprintf('"phenthauproc_version": %s', json_str(as.character(utils::packageVersion("PHENTHAUproc")))),
    sprintf('"def_eclosion": %s',    json_str(DEF_ECLOSION)),
    sprintf('"mailles": %d',         as.integer(n_mailles)),
    sprintf('"avertissement": %s',   json_str(AVERTISSEMENT_PHENTHAU)),
    sprintf('"source": %s',          json_str(SOURCE_PHENTHAU))
  )
  chemin <- file.path(dir_out, "phenthau_meta.json")
  writeLines(c("{", paste0("  ", champs, collapse = ",\n"), "}"), chemin)
  cat(glue("  Metadonnees -> {basename(chemin)} (donnees au {pheno$date_ref})\n"))
  invisible(chemin)
}


# =============================================================================
# 6. SYNTHESE CONSOLE (diagnostic de log, ne doit jamais faire echouer le run)
# =============================================================================

synthese_phenthau <- function(pheno, ras) {
  n <- sum(!is.na(ras$stade_courant))
  med <- function(x) if (all(is.na(x))) "NA" else
    format(as.Date(sprintf("%d-01-01", pheno$annee)) + round(stats::median(x, na.rm = TRUE)) - 1L, "%d/%m")
  cat(glue("\n--- PHENTHAU {pheno$annee} au {pheno$date_ref} : {n} mailles ---\n"))
  cat("  Stade courant (mailles) : ",
      paste(sprintf("%s=%d", TABLE_STADES$code[match(sort(unique(stats::na.omit(ras$stade_courant))), TABLE_STADES$stade)],
                    as.integer(table(stats::na.omit(ras$stade_courant)))), collapse = " | "), "\n", sep = "")
  cat(glue("  Dates medianes : debourrement {med(ras$doy_debourrement)} | eclosion {med(ras$doy_eclosion)} | ",
           "L3 {med(ras$doy_stade_3)} | pupe {med(ras$doy_stade_7)} | adulte {med(ras$doy_stade_8)}\n"))
  cat(glue("  Mortalite L1 : {sum(!is.na(ras$mortalite_pct))} mailles determinees, ",
           "{sum(ras$mortalite_pct >= 50, na.rm = TRUE)} avec >= 50 %\n\n"))
}

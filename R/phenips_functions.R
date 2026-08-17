# =============================================================================
# PHENIPS-Clim — France metropolitaine
# v3 — SAFRAN via API GéoSAS INRAE | Architecture dalles + cache CSV
# =============================================================================


# =============================================================================
# 0. PACKAGES
# =============================================================================
# Seulement ce qui est utilise APRES le chargement du jeu SAFRAN partage
# (safran_common.R deja charge terra/httr2/glue/data.table/lubridate).
# ggplot2/tidyterra/scales ne sont plus necessaires : plus d'export PNG.

library(barrks)
library(sf)


# =============================================================================
# 1. PARAMETRES SPECIFIQUES PHENIPS
# =============================================================================
# ANNEE, DATE_DEB, RST_FRANCE, CRS_L93/CRS_WGS proviennent deja de
# safran_common.R (source avant ce fichier) : pas de redeclaration ici.
# DIR_OUT n'est PAS defini par safran_common.R (qui ne gere que le cache
# SAFRAN, commun) : chaque modele garde le sien, tous pointant vers "data".

DIR_OUT <- "data"
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# 7. HELPERS ACCES OBJET BARRKS
# =============================================================================

inspect_barrks <- function(pheno) {
  cat("\n--- Structure objet barrks (diagnostic) ---\n")
  cat("Classe     :", class(pheno), "\n")
  cat("Noms slots :", paste(names(pheno), collapse = " | "), "\n")
  if (!is.null(pheno$data))    cat("  $data    :", paste(names(pheno$data),    collapse = " | "), "\n")
  if (!is.null(pheno$storage)) cat("  $storage :", paste(names(pheno$storage), collapse = " | "), "\n")
  if (!is.null(pheno$development)) {
    cat("  $development :", paste(names(pheno$development), collapse = " | "), "\n")
  }
  cat("---\n\n")
}

get_btmean <- function(pheno) {
  if (!is.null(pheno[["data"]][["btmean"]]))    return(pheno$data$btmean)
  if (!is.null(pheno[["storage"]][["btmean"]])) return(pheno$storage$btmean)
  if (!is.null(pheno[["btmean"]]))              return(pheno$btmean)
  inspect_barrks(pheno)
  stop("btmean introuvable dans l'objet phenips. Voir structure ci-dessus.")
}

get_gen1 <- function(pheno) {
  if (!is.null(pheno[["development"]][["gen_1"]]))    return(pheno$development$gen_1)
  if (!is.null(pheno[["development"]][["generation_1"]])) return(pheno$development$generation_1)
  if (!is.null(pheno[["gen_1"]]))                     return(pheno$gen_1)
  NULL
}

get_gen_rst <- function(pheno) {
  date_str <- as.character(max(pheno$dates))
  tryCatch(
    get_generations_rst(pheno, date_str),
    error = function(e) {
      message(glue("  [INFO] get_generations_rst({date_str}) : {e$message}"))
      message("  -> Pas de generation complete disponible (normal avant ete).")
      NULL
    }
  )
}


# =============================================================================
# 8. MODELISATION PHENIPS-Clim
# =============================================================================

modeliser_phenips <- function(meteo,
                               exposure = "sunny",
                               scenario = "max") {
  cat(glue(">> Modelisation PHENIPS-Clim [exposure={exposure} | scenario={scenario}]...\n"))
  cat(glue("   {ncell(meteo$tmean)} cellules x {nlyr(meteo$tmean)} jours\n"))

  stopifnot(
    "time(tmean) non defini" = !any(is.na(as.Date(time(meteo$tmean)))),
    "time(tmin) non defini"  = !any(is.na(as.Date(time(meteo$tmin)))),
    "time(tmax) non defini"  = !any(is.na(as.Date(time(meteo$tmax)))),
    "time(rad) non defini"   = !any(is.na(as.Date(time(meteo$rad)))),
    "time(daylength) non defini" = !any(is.na(as.Date(time(meteo$daylength))))
  )

  pheno <- phenology(
    "phenips-clim",
    tmin          = meteo$tmin,
    tmean         = meteo$tmean,
    tmax          = meteo$tmax,
    rad           = meteo$rad,
    daylength     = meteo$daylength,
    exposure      = exposure,
    scenario      = scenario,
    sister_broods = TRUE
  )

  pheno$dates <- as.Date(time(meteo$tmean))
  cat(glue("  Modelisation terminee ({length(pheno$dates)} jours).\n\n"))
  if (getOption("phenips.debug", default = FALSE)) inspect_barrks(pheno)
  pheno
}


# =============================================================================
# 10b. EXPORT GEOJSON DES GENERATIONS
# =============================================================================

produire_geojson_generations <- function(pheno, dir_out = DIR_OUT) {
  cat(">> Vectorisation generations -> GeoJSON...\n")
  rst_gen <- get_gen_rst(pheno)
  if (is.null(rst_gen) || all(is.na(values(rst_gen)))) {
    cat("  Pas de generation complete disponible (normal avant juillet).\n\n")
    return(invisible(NULL))
  }
  date_j      <- format(max(pheno$dates), "%Y%m%d")
  date_calcul <- as.character(max(pheno$dates))
  rst_wgs <- project(rst_gen, "EPSG:4326", method = "near")
  vect_gen <- tryCatch(
    as.polygons(rst_wgs, dissolve = TRUE) |> st_as_sf(),
    error = function(e) { message("  Erreur vectorisation : ", e$message); NULL }
  )
  if (is.null(vect_gen) || nrow(vect_gen) == 0) {
    cat("  Aucun polygone produit.\n\n")
    return(invisible(NULL))
  }
  names(vect_gen)[1] <- "generation"
  vect_gen$label <- dplyr::case_when(
    vect_gen$generation == 1   ~ "Generation 1 complete",
    vect_gen$generation == 1.5 ~ "Brood soeur (gen 1.5)",
    vect_gen$generation == 2   ~ "Generation 2 complete",
    vect_gen$generation == 2.5 ~ "Brood soeur gen 2",
    vect_gen$generation == 3   ~ "Generation 3 complete",
    TRUE ~ paste0("Gen ", vect_gen$generation)
  )
  vect_gen$couleur <- dplyr::case_when(
    vect_gen$generation == 1   ~ "#35B779FF",
    vect_gen$generation == 1.5 ~ "#FDE725FF",
    vect_gen$generation == 2   ~ "#31688EFF",
    vect_gen$generation == 2.5 ~ "#FFA500",
    vect_gen$generation == 3   ~ "#440154FF",
    TRUE ~ "#cccccc"
  )
  vect_gen$date_calcul <- date_calcul
  vect_gen$annee       <- ANNEE
  vect_gen <- sf::st_simplify(vect_gen, dTolerance = 0.005, preserveTopology = TRUE)
  vect_gen <- vect_gen[!sf::st_is_empty(vect_gen), ]
  chemin <- file.path(dir_out, glue("phenips_generations_{date_j}.geojson"))
  sf::st_write(vect_gen, chemin, delete_dsn = TRUE, quiet = TRUE,
               layer_options = "COORDINATE_PRECISION=5")
  n_poly <- nrow(vect_gen)
  taille <- round(file.size(chemin) / 1024, 0)
  cat(glue("  {n_poly} polygones | {taille} Ko -> {basename(chemin)}\n\n"))
  invisible(vect_gen)
}


# =============================================================================
# 10c. EXPORT GEOJSON DATE DE PREMIER ESSAIMAGE (ONSET)
# =============================================================================

produire_geojson_onset <- function(pheno, dir_out = DIR_OUT) {
  cat(">> Vectorisation onset (premier essaimage) -> GeoJSON...\n")
  onset  <- pheno$onset
  dates  <- pheno$dates
  if (is.null(onset) || nlyr(onset) == 0) {
    cat("  Raster onset absent.\n\n")
    return(invisible(NULL))
  }
  date_j <- format(max(dates), "%Y%m%d")
  onset_vals <- values(onset)
  has_onset  <- rowSums(onset_vals > 0, na.rm = TRUE) > 0
  idx_onset <- rep(NA_integer_, nrow(onset_vals))
  if (any(has_onset, na.rm = TRUE)) {
    sub_mat <- onset_vals[has_onset, , drop = FALSE]
    idx_onset[has_onset] <- max.col(sub_mat > 0, ties.method = "first")
  }
  doy_onset_vec <- rep(NA_real_, length(idx_onset))
  valides <- !is.na(idx_onset)
  doy_onset_vec[valides] <- as.integer(format(dates[idx_onset[valides]], "%j"))
  rst_onset_doy <- RST_FRANCE
  values(rst_onset_doy) <- doy_onset_vec
  rst_onset_doy[is.na(rst_onset_doy)] <- NA
  rst_wgs <- project(rst_onset_doy, "EPSG:4326", method = "near")
  vect_onset <- tryCatch(
    as.polygons(rst_wgs, dissolve = TRUE) |> st_as_sf(),
    error = function(e) { message("  Erreur vectorisation onset : ", e$message); NULL }
  )
  if (is.null(vect_onset) || nrow(vect_onset) == 0) {
    cat("  Aucun polygone onset produit (pas encore d essaimage).\n\n")
    return(invisible(NULL))
  }
  names(vect_onset)[1] <- "doy"
  origine <- as.Date(glue("{ANNEE}-01-01"))
  vect_onset$date_onset   <- format(origine + vect_onset$doy - 1L, "%d/%m/%Y")
  vect_onset$date_iso     <- format(origine + vect_onset$doy - 1L, "%Y-%m-%d")
  vect_onset$semaine      <- as.integer(format(origine + vect_onset$doy - 1L, "%V"))
  vect_onset$date_calcul  <- as.character(max(dates))
  vect_onset$annee        <- ANNEE
  vect_onset$doy <- NULL
  vect_onset <- sf::st_simplify(vect_onset, dTolerance = 0.005, preserveTopology = TRUE)
  vect_onset <- vect_onset[!sf::st_is_empty(vect_onset), ]
  chemin <- file.path(dir_out, glue("phenips_onset_{date_j}.geojson"))
  sf::st_write(vect_onset, chemin, delete_dsn = TRUE, quiet = TRUE,
               layer_options = "COORDINATE_PRECISION=5")
  n_poly <- nrow(vect_onset)
  taille <- round(file.size(chemin) / 1024, 0)
  n_dates <- length(unique(vect_onset$date_onset))
  cat(glue("  {n_poly} polygones | {n_dates} dates uniques | {taille} Ko -> {basename(chemin)}\n\n"))
  invisible(vect_onset)
}


# =============================================================================
# 11. SYNTHESE CONSOLE (diagnostic dans les logs GitHub Actions)
# =============================================================================

synthese_console <- function(pheno) {
  cat("\n", strrep("=", 60), "\n")
  cat("  SYNTHESE PHENIPS-Clim — France metropolitaine\n")
  cat(strrep("=", 60), "\n")

  date_j  <- max(pheno$dates)
  btmean  <- get_btmean(pheno)
  idx_j   <- nlyr(btmean)

  bt_vals  <- values(btmean[[idx_j]])
  n_valide <- sum(!is.na(bt_vals))
  n_total  <- ncell(RST_FRANCE)

  cat(glue("  Date      : {format(date_j, '%d/%m/%Y')}\n"))
  cat(glue("  Cellules  : {n_valide} valides / {n_total} total\n\n"))

  onset_vals  <- values(pheno$onset)
  n_essaime   <- sum(rowSums(onset_vals > 0, na.rm = TRUE) > 0, na.rm = TRUE)
  pct_essaime <- round(100 * n_essaime / n_valide, 1)
  cat(glue("  Essaimage : {n_essaime} cellules ({pct_essaime}%) ont essaime\n\n"))

  cat("  --- Temperature sous ecorce (btmean) ---\n")
  cat(glue("  Moyenne : {round(mean(bt_vals, na.rm=TRUE), 1)} C\n"))
  cat(glue("  Min     : {round(min(bt_vals,  na.rm=TRUE), 1)} C\n"))
  cat(glue("  Max     : {round(max(bt_vals,  na.rm=TRUE), 1)} C\n"))
  n_au_dessus <- sum(bt_vals >= 16.5, na.rm = TRUE)
  cat(glue("  >= seuil essaimage (16.5C) : {n_au_dessus} cellules",
           " ({round(100*n_au_dessus/n_valide,1)}%)\n\n"))

  gen1 <- get_gen1(pheno)
  if (!is.null(gen1)) {
    idx_g   <- nlyr(gen1)
    g1_vals <- values(gen1[[idx_g]])
    g1_mean <- mean(g1_vals, na.rm = TRUE) * 100
    cat(glue("  Dev. gen1 moyen : {round(g1_mean, 1)}%\n"))
    n_complete <- sum(g1_vals >= 1.0, na.rm = TRUE)
    cat(glue("  Gen1 complete   : {n_complete} cellules",
             " ({round(100*n_complete/n_valide,1)}%)\n\n"))
  } else {
    cat("  [INFO] Gen1 non disponible (normal avant mi-saison).\n\n")
  }

  cat(strrep("=", 60), "\n\n")
}

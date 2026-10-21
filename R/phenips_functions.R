# =============================================================================
# PHENIPS-Clim — France metropolitaine
# v4 — SAFRAN via API GeoSAS INRAE | scenarios sunny + shaded | exports GeoJSON
# =============================================================================
# Changements v4 (par rapport a v3) :
#   [P1] codes de generation barrks ("1s", "2s", "-") correctement traduits :
#        avant, ils tombaient tous dans le gris #cccccc du case_when.
#   [P3] max.col() : les NA sont forces a FALSE (une maille avec un jour SAFRAN
#        manquant perdait sa date d'essaimage).
#   [P7] deux scenarios publies : exposition sud (sunny) et a l'ombre (shaded).
#   [P8] indice continu de generations cumulees (en plus de la classe discrete).
#   [P2] "latest" = le fichier REELLEMENT ecrit, jamais "le plus recent par
#        date de modification" ; supprime s'il n'y a rien a publier.
#   [P15] purge des fichiers dates (on garde les 4 derniers par couche).
#   [P16] champ `source` dans les GeoJSON + fichier phenips_meta.json.
# =============================================================================


# =============================================================================
# 0. PACKAGES
# =============================================================================
# Seulement ce qui est utilise APRES le chargement du jeu SAFRAN partage
# (safran_common.R deja charge terra/httr2/glue/data.table/lubridate).

library(barrks)
library(sf)


# =============================================================================
# 1. PARAMETRES SPECIFIQUES PHENIPS
# =============================================================================
# ANNEE, DATE_DEB, RST_FRANCE, CRS_L93/CRS_WGS proviennent deja de
# safran_common.R (source avant ce fichier) : pas de redeclaration ici.

DIR_OUT <- "data"
dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)

SOURCE_DONNEES <- paste0("Source : Meteo-France (donnees SAFRAN) via GeoSAS INRAE ; ",
                         "modele PHENIPS-Clim (package barrks)")

# Deux scenarios d'exposition. Aucun n'est "le pire cas" partout : en region
# mediterraneenne, l'ecorce exposee au sud depasse souvent 38.9 C (seuil
# au-dela duquel PHENIPS arrete le developpement), si bien que "shaded" peut
# donner PLUS de generations que "sunny" autour d'Avignon.
VARIANTES_PHENIPS <- list(
  sunny  = list(exposure = "sunny",  scenario = "max", libelle = "Ecorce exposee au sud"),
  shaded = list(exposure = "shaded", scenario = "max", libelle = "Ecorce a l'ombre")
)

# Les anciens noms de fichiers (sans suffixe de scenario) continuent d'etre
# ecrits, = scenario "sunny", pour que l'ancienne version de l'app Shiny deja
# deployee ne fige pas pendant la transition. A passer a FALSE une fois la
# nouvelle app deployee.
COMPAT_ANCIENS_NOMS <- TRUE

# Table des generations barrks. `code` : numerique, la couvee soeur vaut +0.5.
# Rampe MONOTONE jaune -> rouge sombre (plus c'est fonce, plus c'est avance).
TABLE_GENERATIONS <- data.frame(
  code    = c(0, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5),
  texte   = c("-", "1", "1s", "2", "2s", "3", "3s", "4", "4s"),
  label   = c("Aucune generation achevee",
              "Gen. 1 achevee", "Gen. 1 + couvee soeur",
              "Gen. 2 achevee", "Gen. 2 + couvee soeur",
              "Gen. 3 achevee", "Gen. 3 + couvee soeur",
              "Gen. 4 achevee", "Gen. 4 + couvee soeur"),
  couleur = c("#E8E8E8", "#FFFFB2", "#FED976", "#FD8D3C", "#E31A1C",
              "#BD0026", "#800026", "#4D0015", "#2B000D"),
  stringsAsFactors = FALSE
)

# Indice continu : classes de 0.25 generation, de 0 a 4, couleurs FIXES
# (la meme valeur a toujours la meme couleur, d'une semaine a l'autre).
PAS_INDICE     <- 0.25
INDICE_MAX     <- 4
PALETTE_INDICE <- grDevices::colorRampPalette(
  c("#FFFFCC", "#FED976", "#FD8D3C", "#E31A1C", "#800026")
)(INDICE_MAX / PAS_INDICE + 1)


# =============================================================================
# 2. HELPERS ACCES OBJET BARRKS
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

#' Raster NUMERIQUE des generations a une date (0, 1, 1.5, 2, 2.5...).
#' [P1] categorical = FALSE : la couvee soeur vaut +0.5 (doc barrks). Si la
#' version installee ne connait pas cet argument, repli sur le raster
#' categoriel ("1s"...), que normaliser_generation() sait aussi lire.
get_gen_rst <- function(pheno, date_ref = max(pheno$dates)) {
  date_str <- as.character(date_ref)
  tryCatch(
    tryCatch(get_generations_rst(pheno, date_str, categorical = FALSE),
             error = function(e) get_generations_rst(pheno, date_str)),
    error = function(e) {
      message(glue("  [INFO] get_generations_rst({date_str}) : {conditionMessage(e)}"))
      message("  -> Pas de generation complete disponible (normal avant ete).")
      NULL
    }
  )
}


# =============================================================================
# 3. NOMENCLATURE DES GENERATIONS  [P1]
# =============================================================================

#' "-" -> 0 ; "1" / 1 -> 1 ; "1s" / 1.5 -> 1.5 ; "2s" -> 2.5 ...
#' Accepte aussi bien les codes texte de barrks que les nombres.
normaliser_generation <- function(x) {
  x_chr <- trimws(as.character(x))
  soeur <- grepl("s$", x_chr, ignore.case = TRUE)
  base  <- suppressWarnings(as.numeric(sub("s$", "", x_chr, ignore.case = TRUE)))
  out   <- ifelse(soeur, base + 0.5, base)
  out[x_chr %in% "-"] <- 0
  out
}

#' Libelle, couleur, ordre et code texte pour un vecteur de codes numeriques.
#' Un code absent de TABLE_GENERATIONS ressort en MAGENTA "Code inconnu" :
#' une future version de barrks qui ajouterait un code se verra tout de suite.
habiller_generations <- function(code) {
  i <- match(round(code * 2) / 2, TABLE_GENERATIONS$code)
  inconnus <- is.na(i)
  if (any(inconnus)) {
    message("::warning::Codes de generation barrks inconnus : ",
            paste(unique(code[inconnus]), collapse = ", "))
  }
  list(
    generation = ifelse(inconnus, as.character(code), TABLE_GENERATIONS$texte[i]),
    label      = ifelse(inconnus, paste0("Code inconnu : ", code), TABLE_GENERATIONS$label[i]),
    couleur    = ifelse(inconnus, "#FF00FF", TABLE_GENERATIONS$couleur[i]),
    ordre      = ifelse(inconnus, NA_integer_, i - 1L)
  )
}


# =============================================================================
# 4. MODELISATION PHENIPS-Clim
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
# 5. OUTILS D'EXPORT GEOJSON
# =============================================================================

#' Raster -> polygones dissous (WGS84), premiere colonne renommee `nom_champ`.
vectoriser_raster <- function(rst, nom_champ) {
  rst_wgs <- project(rst, "EPSG:4326", method = "near")
  v <- tryCatch(
    as.polygons(rst_wgs, dissolve = TRUE) |> st_as_sf(),
    error = function(e) { message("  Erreur vectorisation : ", conditionMessage(e)); NULL }
  )
  if (is.null(v) || nrow(v) == 0) return(NULL)
  names(v)[1] <- nom_champ
  v
}

#' Metadonnees communes a tous les GeoJSON PHENIPS. [P13][P16]
ajouter_meta <- function(v, variante, date_ref) {
  v$date_calcul <- as.character(date_ref)
  v$annee       <- ANNEE
  v$scenario    <- variante
  v$exposition  <- VARIANTES_PHENIPS[[variante]]$libelle
  v$source      <- SOURCE_DONNEES
  v
}

#' Simplifie, ecrit, et RETOURNE LE CHEMIN ecrit (utilise par publier_latest).
ecrire_geojson <- function(v, nom_fichier, dir_out = DIR_OUT) {
  v <- sf::st_simplify(v, dTolerance = 0.005, preserveTopology = TRUE)
  v <- v[!sf::st_is_empty(v), ]
  chemin <- file.path(dir_out, nom_fichier)
  sf::st_write(v, chemin, delete_dsn = TRUE, quiet = TRUE,
               layer_options = "COORDINATE_PRECISION=5")
  cat(glue("  {nrow(v)} polygones | {round(file.size(chemin) / 1024, 0)} Ko -> {nom_fichier}\n\n"))
  chemin
}


# =============================================================================
# 6. EXPORT : GENERATIONS (classes)  [P1]
# =============================================================================

produire_geojson_generations <- function(pheno, variante = "sunny", dir_out = DIR_OUT) {
  cat(glue(">> [{variante}] Vectorisation generations -> GeoJSON...\n"))
  date_ref <- max(pheno$dates)
  rst_gen  <- get_gen_rst(pheno, date_ref)
  if (is.null(rst_gen) || all(is.na(values(rst_gen)))) {
    cat("  Pas de generation complete disponible (normal avant juillet).\n\n")
    return(invisible(NULL))
  }
  v <- vectoriser_raster(rst_gen, "generation_brute")
  if (is.null(v)) { cat("  Aucun polygone produit.\n\n"); return(invisible(NULL)) }

  hab <- habiller_generations(normaliser_generation(v$generation_brute))
  v$generation_brute <- NULL
  v$generation <- hab$generation
  v$label      <- hab$label
  v$couleur    <- hab$couleur
  v$ordre      <- hab$ordre
  v <- v[order(v$ordre), ]
  v <- ajouter_meta(v, variante, date_ref)

  invisible(ecrire_geojson(v, glue("phenips_generations_{variante}_{format(date_ref, '%Y%m%d')}.geojson"),
                           dir_out))
}


# =============================================================================
# 7. EXPORT : DATE DE PREMIER ESSAIMAGE (onset)  [P3]
# =============================================================================

produire_geojson_onset <- function(pheno, variante = "sunny", dir_out = DIR_OUT) {
  cat(glue(">> [{variante}] Vectorisation onset (premier essaimage) -> GeoJSON...\n"))
  onset <- pheno$onset
  dates <- pheno$dates
  if (is.null(onset) || nlyr(onset) == 0) {
    cat("  Raster onset absent.\n\n")
    return(invisible(NULL))
  }
  date_ref   <- max(dates)
  onset_vals <- values(onset)
  has_onset  <- rowSums(onset_vals > 0, na.rm = TRUE) > 0
  idx_onset  <- rep(NA_integer_, nrow(onset_vals))
  if (any(has_onset, na.rm = TRUE)) {
    sous_bin <- onset_vals[has_onset, , drop = FALSE] > 0
    # [P3] Sans cette ligne, une maille avec un seul jour SAFRAN manquant
    # donnait max.col() = NA et disparaissait de la carte.
    sous_bin[is.na(sous_bin)] <- FALSE
    idx_onset[has_onset] <- max.col(sous_bin, ties.method = "first")
  }
  doy_onset_vec <- rep(NA_real_, length(idx_onset))
  valides <- !is.na(idx_onset)
  doy_onset_vec[valides] <- as.integer(format(dates[idx_onset[valides]], "%j"))
  rst_onset_doy <- RST_FRANCE
  values(rst_onset_doy) <- doy_onset_vec

  v <- vectoriser_raster(rst_onset_doy, "doy")
  if (is.null(v)) {
    cat("  Aucun polygone onset produit (pas encore d'essaimage).\n\n")
    return(invisible(NULL))
  }
  origine <- as.Date(glue("{ANNEE}-01-01"))
  v$date_onset <- format(origine + v$doy - 1L, "%d/%m/%Y")
  v$date_iso   <- format(origine + v$doy - 1L, "%Y-%m-%d")
  v$semaine    <- as.integer(format(origine + v$doy - 1L, "%V"))
  v$doy <- NULL
  v <- ajouter_meta(v, variante, date_ref)

  invisible(ecrire_geojson(v, glue("phenips_onset_{variante}_{format(date_ref, '%Y%m%d')}.geojson"),
                           dir_out))
}


# =============================================================================
# 8. EXPORT : INDICE CONTINU DE GENERATIONS CUMULEES  [P8]
# =============================================================================
# Indice = somme, sur les generations filiales 1..n, du developpement de chaque
# generation ramene a [0, 1]. Exemple : 2.40 = 2 generations achevees + 40 %
# de la 3e. Contrairement a la classe discrete, il distingue deux mailles
# situees du meme cote d'un seuil. Contrairement au seul developpement de la
# generation 1, il ne sature pas a 100 % en milieu d'ete.
# Les couvees soeurs ne sont pas comptees (cohorte distincte de la filiale).

produire_geojson_indice <- function(pheno, variante = "sunny", dir_out = DIR_OUT) {
  cat(glue(">> [{variante}] Indice continu de generations -> GeoJSON...\n"))
  date_ref <- max(pheno$dates)

  gens <- tryCatch(as.numeric(unlist(prop_filial_generations(pheno))), error = function(e) NULL)
  if (is.null(gens) || length(gens) == 0 || anyNA(gens)) gens <- 1:3
  gens <- sort(unique(gens[gens == floor(gens)]))

  couches <- lapply(gens, function(g) {
    tryCatch({
      r <- get_development_rst(pheno, generation = g, dates = as.character(date_ref))
      # -1 = "generation pas encore presente" (doc barrks) -> 0 ; plafond a 1
      clamp(r[[1]], lower = 0, upper = 1, values = TRUE)
    }, error = function(e) {
      message(glue("  [INFO] developpement generation {g} indisponible : {conditionMessage(e)}"))
      NULL
    })
  })
  couches <- Filter(Negate(is.null), couches)
  if (length(couches) == 0) {
    cat("  Aucun developpement disponible.\n\n")
    return(invisible(NULL))
  }
  indice <- Reduce(`+`, couches)
  classe <- floor(indice / PAS_INDICE) * PAS_INDICE

  v <- vectoriser_raster(classe, "borne_inf")
  if (is.null(v)) { cat("  Aucun polygone produit.\n\n"); return(invisible(NULL)) }

  n_max <- length(couches)
  v$borne_inf <- as.numeric(v$borne_inf)
  v$borne_sup <- pmin(v$borne_inf + PAS_INDICE, n_max)
  v$label <- ifelse(v$borne_inf >= n_max,
                    sprintf("%d generations achevees", n_max),
                    sprintf("%.2f - %.2f generations", v$borne_inf, v$borne_sup))
  v$ordre   <- as.integer(round(v$borne_inf / PAS_INDICE))
  v$couleur <- PALETTE_INDICE[pmin(v$ordre + 1L, length(PALETTE_INDICE))]
  v <- v[order(v$ordre), ]
  v <- ajouter_meta(v, variante, date_ref)

  invisible(ecrire_geojson(v, glue("phenips_indice_{variante}_{format(date_ref, '%Y%m%d')}.geojson"),
                           dir_out))
}


# =============================================================================
# 9. PUBLICATION : copie "latest", metadonnees, purge  [P2][P15][P16]
# =============================================================================

#' Copie `chemin` (le fichier qu'on vient d'ecrire) vers le nom stable lu par
#' l'app. [P2] Si rien n'a ete produit (`chemin` NULL), le nom stable est
#' SUPPRIME : avant, on recopiait "le plus recent fichier date du dossier",
#' donc en janvier l'app affichait la carte d'aout de l'annee precedente.
publier_latest <- function(chemin, nom_stable, dir_out = DIR_OUT) {
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

#' Garde les `garder` dates les plus recentes pour chaque couche PHENIPS
#' (phenips_<couche>_<AAAAMMJJ>.geojson), supprime les autres. Les fichiers
#' *_latest.geojson ne correspondent pas au motif : ils ne sont jamais touches.
purger_fichiers_dates <- function(dir_out = DIR_OUT, garder = 4L) {
  f <- list.files(dir_out, pattern = "^phenips_.*_[0-9]{8}\\.geojson$")
  if (length(f) == 0) return(invisible(character(0)))
  date    <- sub("^.*_([0-9]{8})\\.geojson$", "\\1", f)
  prefixe <- sub("_[0-9]{8}\\.geojson$", "", f)
  a_supprimer <- character(0)
  for (p in unique(prefixe)) {
    dates_p <- sort(unique(date[prefixe == p]), decreasing = TRUE)
    if (length(dates_p) > garder) {
      vieux <- dates_p[-seq_len(garder)]
      a_supprimer <- c(a_supprimer, f[prefixe == p & date %in% vieux])
    }
  }
  if (length(a_supprimer) > 0) {
    file.remove(file.path(dir_out, a_supprimer))
    cat(glue("  Purge : {length(a_supprimer)} ancien(s) fichier(s) date(s) supprime(s).\n"))
  }
  invisible(a_supprimer)
}

#' phenips_meta.json : date des donnees, version de barrks, etc. Lu par l'app
#' pour afficher la fraicheur meme quand une couche est absente. [P13]
ecrire_meta_phenips <- function(date_donnees, variantes, dir_out = DIR_OUT) {
  qc <- tryCatch(readRDS("data/cache_meteo/qc_safran.rds"), error = function(e) NULL)
  comblees <- if (!is.null(qc) && !is.null(qc$n_cellules_comblees) &&
                  !is.na(qc$n_cellules_comblees)) qc$n_cellules_comblees else NA
  json_str <- function(x) paste0('"', gsub('"', '\\\\"', as.character(x)), '"')
  champs <- c(
    sprintf('"date_donnees": %s',     json_str(date_donnees)),
    sprintf('"date_generation": %s',  json_str(format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"))),
    sprintf('"barrks_version": %s',   json_str(as.character(utils::packageVersion("barrks")))),
    sprintf('"variantes": [%s]',      paste(json_str(variantes), collapse = ", ")),
    sprintf('"mailles_comblees": %s', if (is.na(comblees)) "null" else as.character(comblees)),
    sprintf('"source": %s',           json_str(SOURCE_DONNEES))
  )
  chemin <- file.path(dir_out, "phenips_meta.json")
  writeLines(c("{", paste0("  ", champs, collapse = ",\n"), "}"), chemin)
  cat(glue("  Metadonnees -> {basename(chemin)} (donnees au {date_donnees})\n"))
  invisible(chemin)
}


# =============================================================================
# 10. SYNTHESE CONSOLE (diagnostic dans les logs GitHub Actions)
# =============================================================================

synthese_console <- function(pheno) {
  cat("\n", strrep("=", 60), "\n")
  cat("  SYNTHESE PHENIPS-Clim - France metropolitaine\n")
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

# =============================================================================
# PIPELINE PHENTHAU (Thaumetopoea processionea) -- reutilise le jeu SAFRAN partage
# -----------------------------------------------------------------------------
# Entrees  : data/cache_meteo/meteo_national.rds   (fetch partage, 1er janvier N -> J-5)
#            data/cache_meteo/meteo_amont.rds      (R/fetch_safran_amont.R, sept-dec N-1)
#            data/masques/masque_chene.csv         (BD Foret V2, "Chenes decidus")
# Sorties  : phenthau_{stades,l3,eclosion,debourrement,mortalite}_<date>.geojson
#            -> copie vers *_latest.geojson, phenthau_meta.json, phenthau_<annee>_L93.tif
# Fonctions: R/phenthau_functions.R
#
# MODE TEST (petite emprise, plusieurs annees), sans le fetch national ni l'amont :
#   PHENTHAU_TEST_EMPRISE="xmin,ymin,xmax,ymax;xmin,ymin,xmax,ymax"  (L93)
#   PHENTHAU_TEST_ANNEES="2023,2024,2025"
#   PHENTHAU_DIR_OUT="data/diag/phenthau_test"   (sinon "data")
# Les donnees SAFRAN de l'emprise sont alors telechargees avec les fonctions de la brique.
# =============================================================================

source("R/safran_common.R")
source("R/phenthau_functions.R")
source("R/masques.R")

EMPRISES  <- Sys.getenv("PHENTHAU_TEST_EMPRISE", "")
MODE_TEST <- nzchar(EMPRISES)

# --- MODE TEST : SAFRAN d'une petite emprise sur plusieurs annees, via la brique partagee ---
meteo_test_emprise <- function(emprises, annees) {
  deb <- as.Date(sprintf("%d-09-01", min(annees) - 1L))
  fin <- as.Date(sprintf("%d-09-30", max(annees)))
  jan <- if (lubridate::year(deb) < lubridate::year(fin))
    as.Date(sprintf("%d-01-01", (lubridate::year(deb) + 1L):lubridate::year(fin))) else as.Date(character(0))
  bornes <- c(deb, jan, fin + 1L)
  rects <- lapply(strsplit(emprises, ";")[[1]], function(s) as.numeric(strsplit(s, ",")[[1]]))

  # les fonctions de la brique lisent la periode et l'annee de cache dans des variables globales
  noms <- c("DATE_DEB", "DATE_FIN_SAFRAN", "ANNEE", "DIR_CACHE")
  sauve <- mget(noms, envir = globalenv())
  on.exit(list2env(sauve, envir = globalenv()), add = TRUE)
  assign("DATE_DEB", deb, envir = globalenv())
  assign("DATE_FIN_SAFRAN", fin, envir = globalenv())
  assign("ANNEE", as.integer(max(annees)), envir = globalenv())
  assign("DIR_CACHE", file.path(tempdir(), "safran_test_phenthau"), envir = globalenv())
  dir.create(DIR_CACHE, showWarnings = FALSE, recursive = TRUE)

  k <- 0L
  for (r in rects) for (i in seq_len(length(bornes) - 1L)) {
    k <- k + 1L
    f <- file.path(DIR_CACHE, sprintf("safran_v2_patch_%03d_%d.csv", k, ANNEE))
    res <- telecharger_bbox(r[1], r[2], r[3], r[4], f, date_deb = bornes[i], date_fin = bornes[i + 1L] - 1L)
    cat(glue("  [test] emprise {k}: {bornes[i]} -> {bornes[i + 1L] - 1L} : {res$statut} (HTTP {res$code})\n"))
    if (res$statut != "ok") stop("Telechargement SAFRAN de test en echec", call. = FALSE)
  }
  dt <- lire_cache_complet(dalles = data.frame(id = 999L))
  assembler_rasters_nationaux(dt)[c("tmean", "tmin", "tmax")]
}

if (MODE_TEST) {
  annees <- as.integer(strsplit(Sys.getenv("PHENTHAU_TEST_ANNEES", as.character(ANNEE - 1L)), ",")[[1]])
  cat(glue(">> PHENTHAU MODE TEST : annees {paste(annees, collapse = ', ')}, sortie dans {DIR_OUT}\n\n"))
  meteo <- meteo_test_emprise(EMPRISES, annees)
  amont <- NULL
} else {
  annees <- ANNEE
  meteo <- lapply(readRDS("data/cache_meteo/meteo_national.rds"), terra::unwrap)
  amont <- lapply(readRDS("data/cache_meteo/meteo_amont.rds"),    terra::unwrap)
}

NOMS_COUCHES <- c("stades", "l3", "eclosion", "debourrement", "mortalite")

# --- Avant le 1er fevrier : le package ne produit rien. On retire les cartes de l'annee precedente. ---
date_brute <- max(as.Date(time(meteo$tmean)))
if (!MODE_TEST && date_brute < as.Date(sprintf("%d-02-01", ANNEE))) {
  cat(glue(">> PHENTHAU : donnees au {date_brute}, avant le 1er fevrier : rien a calculer.\n"))
  for (n in NOMS_COUCHES) publier_latest_phenthau(NULL, glue("phenthau_{n}_latest.geojson"))
  quit(save = "no", status = 0)
}

# --- Masque chenaie (v3) : on ne garde que les mailles qui croisent des Chenes decidus ---
masque <- charger_masque("chene", meteo$tmean)
mini   <- if (MODE_TEST) 1L else 100L
meteo  <- appliquer_masque(meteo, masque, nom = "chene", min_mailles = mini)
if (!is.null(amont)) amont <- appliquer_masque(amont, masque, nom = "chene (amont)", min_mailles = mini)

for (a in annees) {
  cat(glue("\n######## PHENTHAU {a} ########\n\n"))
  serie <- assembler_serie_phenthau(meteo, amont, annee = a)
  cat(glue(">> Serie SAFRAN : {min(as.Date(time(serie$tmean)))} -> {max(as.Date(time(serie$tmean)))} ",
           "({nlyr(serie$tmean)} jours)\n"))

  pheno <- modeliser_phenthau(serie, annee = a)
  cat(glue(">> Modele PHENTHAUproc : stades du {min(pheno$dates)} au {pheno$date_ref}\n\n"))
  ras <- calculer_rasters_phenthau(pheno)

  ch <- list(
    stades       = produire_geojson_stades(pheno, ras),
    l3           = produire_geojson_date(pheno, ras, "doy_stade_3",      "l3",           "date_l3"),
    eclosion     = produire_geojson_date(pheno, ras, "doy_eclosion",     "eclosion",     "date_eclosion"),
    debourrement = produire_geojson_date(pheno, ras, "doy_debourrement", "debourrement", "date_debourrement"),
    mortalite    = produire_geojson_mortalite(pheno, ras)
  )
  for (n in NOMS_COUCHES) publier_latest_phenthau(ch[[n]], glue("phenthau_{n}_latest.geojson"))

  # GeoTIFF et synthese : jamais bloquants pour la publication
  tryCatch(exporter_geotiff_phenthau(pheno, ras),
           error = function(e) message("  [INFO] GeoTIFF ignore : ", conditionMessage(e)))
  tryCatch(synthese_phenthau(pheno, ras),
           error = function(e) message("  [INFO] synthese ignoree : ", conditionMessage(e)))

  n_mailles <- sum(!is.na(ras$stade_courant))
  dernier   <- list(date_ref = pheno$date_ref, annee = a)
  rm(serie, pheno, ras); invisible(gc())
}

# meta et purge : sur la derniere annee traitee
ecrire_meta_phenthau(dernier, n_mailles)
purger_fichiers_phenthau(garder = 4L)

cat("\n>> PHENTHAU termine.\n\n")

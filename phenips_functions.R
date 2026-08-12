# =============================================================================
# PHENIPS-Clim — France metropolitaine
# v3 — SAFRAN via API GéoSAS INRAE | Architecture dalles + cache CSV
# -----------------------------------------------------------------------------
# CORRECTIONS v2 (bugs logiques barrks) :
#   [BUG1] pheno$dates absent de l'objet barrks -> ajoute dans modeliser_phenips()
#   [BUG2] pheno$data$btmean : acces securise via helper get_btmean()
#   [BUG3] tryCatch silencieux sur get_generations_rst() -> message() ajoute
#   [BUG4] apply() ligne par ligne sur onset -> remplace par max.col() vectorise
#   [BUG5] Verification SSI_Q range avant pipeline (commentaire + arret si aberrant)
#   [BONUS] pheno$development$gen_1 securise via helper get_gen1()
#   [BONUS] synthese_console() robustifiee (ncell depuis btmean, pas RST_FRANCE)
#   [BONUS] Ajout d'un helper inspect_barrks() pour diagnostiquer la structure
#           de l'objet retourne par phenology() selon la version de barrks
#
# CORRECTIONS v3 (bugs API GéoSAS — source : doc officielle geosas.fr) :
#   [BUG6] TMIN_Q n'existe PAS dans SAFRAN -> remplace par TINF_H_Q
#          (temp. min des 24h horaires, ajoutee le 23/01/2024)
#          -> corrige l'erreur de syntaxe API qui corrompait tous les CSV cache
#   [BUG7] SSI_Q est en J/cm2 (PAS en W/m2 comme suppose en v1/v2)
#          Conversion vers Wh/m2 (attendu par barrks) :
#          1 J/cm2 = 10000 J/m2 = 2.7778 Wh/m2
#          SSI_CHECK_MAX ajuste en consequence (mediane ete ~15 J/cm2)
# -----------------------------------------------------------------------------
# Source meteo : API OGC EDR GéoSAS INRAE
#   URL        : https://api.geosas.fr/edr/collections/safran-isba/cube
#   Format     : CSV (bbox en L93, pas de requete point par point)
#   Variables  : T_Q (Tmoy), TINF_H_Q (Tmin), TSUP_H_Q (Tmax), SSI_Q (J/cm2)
#   Resolution : grille SAFRAN 8x8 km, France metropolitaine
#   Disponible : jusqu'a J-5 environ — mise a jour mensuelle (debut de mois)
#
# Modele pheno : PHENIPS-Clim via barrks >= 1.1.2
#   - Calcul temperature sous ecorce (btmean/btmax) via SSI_Q
#   - Onset (premier essaimage), developpement larvaire, generations
#   - Fix terra::time() obligatoire pour barrks (maintenu)
#
# References  : Baier et al. 2007 (PHENIPS), Ogris et al. 2019 (PHENIPS-Clim)
#               GéoSAS INRAE : https://geosas.fr
# Auteur      : CNPF / CRPF AuRA — avril 2026
# =============================================================================


# =============================================================================
# 0. PACKAGES
# =============================================================================

# install.packages(c("barrks", "terra", "httr2", "ggplot2", "tidyterra",
#                    "lubridate", "glue", "sf", "scales", "data.table"))

library(barrks)
library(terra)
library(httr2)
library(ggplot2)
library(tidyterra)
library(lubridate)
library(glue)
library(sf)
library(scales)
library(data.table)


# =============================================================================
# 1. PARAMETRES GLOBAUX
# =============================================================================

# --- Periode ---
ANNEE    <- as.integer(format(Sys.Date(), "%Y"))
DATE_DEB <- as.Date(glue("{ANNEE}-01-01"))

# SAFRAN disponible jusqu'a ~J-5 (mise a jour mensuelle)
DATE_FIN_SAFRAN <- Sys.Date() - 5
if (DATE_FIN_SAFRAN < DATE_DEB) DATE_FIN_SAFRAN <- DATE_DEB

# --- API GéoSAS ---
GEOSAS_URL  <- "https://api.geosas.fr/edr/collections/safran-isba/cube"
GEOSAS_CRS  <- "EPSG:2154"
# Variables PHENIPS-Clim : Tmoy, Tmin, Tmax, Rayonnement global
# [BUG6] La variable Tmin dans SAFRAN est TINF_H_Q (ajoutee 23/01/2024)
#        TMIN_Q n'existe PAS -> l'API retournait une erreur JSON corrompant le cache
GEOSAS_VARS <- "T_Q,TINF_H_Q,TSUP_H_Q,SSI_Q"

# --- Grille SAFRAN France metropolitaine (L93) ---
FRANCE_XMIN <- 100000
FRANCE_XMAX <- 1200000
FRANCE_YMIN <- 6050000
FRANCE_YMAX <- 7150000
SAFRAN_RES  <- 8000

# Decoupage en dalles 200x200 km pour eviter les timeouts API
DALLE_KM <- 200
DALLE_M  <- DALLE_KM * 1000

CRS_L93 <- "EPSG:2154"
CRS_WGS <- "EPSG:4326"

# Conversion rayonnement :
#   SAFRAN SSI_Q en J/cm2 (cumul quotidien) — confirme par physique Baier 2007
#   Valeurs brutes observees : ~627 J/cm2 (jan-mars), ~1765 J/cm2 (avril),
#                              ~2000-2500 J/cm2 (juillet)
#   barrks phenips-clim attend Wh/m2 -> conversion : 1 J/cm2 = 10000/3600 Wh/m2
#   Verification coefficient Baier (btmax = tmax + 0.000814 * I) :
#     Juillet avec factor=2.778 : 2300 J/cm2 -> 6389 Wh/m2 -> deltaT = +5.2 C (coherent)
#     Juillet avec factor=1     : 2500 Wh/m2 seulement     -> deltaT = +2.0 C (trop faible)
#   [BUG8] v3 avait mis WHm2_FACTOR=1 par erreur -> mode shaded de facto
WHm2_FACTOR <- 10000 / 3600   # J/cm2 -> Wh/m2 = 2.7778

# Seuil de controle SSI_Q brut (J/cm2) : mediane ete ~2000-2500 J/cm2
SSI_CHECK_MAX <- 4000   # J/cm2 — au-dessus = unite suspecte

# --- Repertoires ---
DIR_OUT   <- "output_phenips_france"
DIR_CACHE <- file.path(DIR_OUT, "cache_safran")
dir.create(DIR_OUT,   showWarnings = FALSE, recursive = TRUE)
dir.create(DIR_CACHE, showWarnings = FALSE, recursive = TRUE)

cat("=== PHENIPS-Clim — France metropolitaine (SAFRAN) v2 ===\n")
cat(glue("  Periode    : {DATE_DEB} -> {DATE_FIN_SAFRAN}\n"))
cat(glue("  Variables  : {GEOSAS_VARS}\n"))
cat(glue("  Resolution : {SAFRAN_RES/1000} km | Dalles : {DALLE_KM}x{DALLE_KM} km\n"))
cat(glue("  Cache      : {DIR_CACHE}/\n\n"))


# =============================================================================
# 2. GRILLE DE DALLES FRANCE
# =============================================================================

creer_dalles <- function() {
  xs <- seq(FRANCE_XMIN, FRANCE_XMAX, by = DALLE_M)
  ys <- seq(FRANCE_YMIN, FRANCE_YMAX, by = DALLE_M)
  dalles      <- expand.grid(xmin = xs, ymin = ys)
  dalles$xmax <- pmin(dalles$xmin + DALLE_M, FRANCE_XMAX)
  dalles$ymax <- pmin(dalles$ymin + DALLE_M, FRANCE_YMAX)
  dalles$id   <- seq_len(nrow(dalles))
  dalles
}

DALLES <- creer_dalles()
cat(glue("  Dalles France : {nrow(DALLES)}\n\n"))


# =============================================================================
# 3. RASTER DE REFERENCE FRANCE (grille SAFRAN native, L93)
# =============================================================================

creer_raster_france <- function() {
  rast(
    xmin = FRANCE_XMIN, xmax = FRANCE_XMAX,
    ymin = FRANCE_YMIN, ymax = FRANCE_YMAX,
    resolution = SAFRAN_RES, crs = CRS_L93
  )
}

RST_FRANCE <- creer_raster_france()
cat(glue("  Raster France : {nrow(RST_FRANCE)}L x {ncol(RST_FRANCE)}C",
         " = {ncell(RST_FRANCE)} cellules\n\n"))


# =============================================================================
# 4. TELECHARGEMENT SAFRAN PAR DALLE (avec cache)
# =============================================================================
# Format CSV retourne par GéoSAS /cube :
#   colonnes typiques : datetime (ou date), x, y, T_Q, TMIN_Q, TSUP_H_Q, SSI_Q
#   x, y : coordonnees centroide tuile SAFRAN en L93
#
# Cache : 1 fichier CSV par dalle et par annee
#   Valide si la derniere date du fichier >= DATE_FIN_SAFRAN
#   Sinon : re-telecharge toute la saison depuis DATE_DEB

chemin_cache <- function(dalle_id) {
  file.path(DIR_CACHE,
            glue("safran_dalle{sprintf('%03d', dalle_id)}_{ANNEE}.csv"))
}

cache_est_valide <- function(dalle_id) {
  f <- chemin_cache(dalle_id)
  if (!file.exists(f)) return(FALSE)
  tryCatch({
    entete <- names(fread(f, nrows = 0))
    if (length(entete) == 0) return(TRUE)   # dalle hors France (vide), OK
    col_t <- entete[grep("^time$|^date$|datetime", tolower(entete))][1]
    if (is.na(col_t)) return(TRUE)
    df <- fread(f, select = col_t)
    if (nrow(df) == 0) return(TRUE)
    max_date <- as.Date(max(df[[col_t]], na.rm = TRUE))
    if (is.na(max_date) || is.infinite(as.numeric(max_date))) return(TRUE)
    max_date >= DATE_FIN_SAFRAN
  }, error = function(e) FALSE)
}

#' Telecharge une dalle SAFRAN et sauvegarde en CSV cache
telecharger_dalle <- function(dalle,
                              vars     = GEOSAS_VARS,
                              date_deb = DATE_DEB,
                              date_fin = DATE_FIN_SAFRAN) {
  bbox    <- glue("{dalle$xmin},{dalle$ymin},{dalle$xmax},{dalle$ymax}")
  periode <- glue("{date_deb}/{date_fin}")

  req <- request(GEOSAS_URL) |>
    req_url_query(
      bbox             = bbox,
      crs              = GEOSAS_CRS,
      `parameter-name` = vars,
      f                = "CSV",
      datetime         = periode
    ) |>
    req_timeout(180) |>
    req_retry(
      max_tries    = 5,
      backoff      = \(i) 30 * i,
      is_transient = \(resp) resp_status(resp) %in% c(429, 503, 504)
    )

  resp <- tryCatch(
    req_perform(req),
    error = function(e) {
      message("  ERREUR dalle ", dalle$id, ": ", e$message)
      NULL
    }
  )

  if (is.null(resp)) return(NULL)
  if (resp_status(resp) != 200) {
    message("  HTTP ", resp_status(resp), " dalle ", dalle$id)
    return(NULL)
  }

  writeBin(resp_body_raw(resp), chemin_cache(dalle$id))
  invisible(chemin_cache(dalle$id))
}

#' Boucle de telechargement sur toutes les dalles (avec cache)
recuperer_safran <- function(dalles = DALLES) {
  cat(glue(">> Telechargement SAFRAN ({DATE_DEB} -> {DATE_FIN_SAFRAN})...\n"))
  n       <- nrow(dalles)
  n_cache <- 0L; n_dl <- 0L; n_err <- 0L

  for (i in seq_len(n)) {
    cat(glue("\r  Dalle {i}/{n} | cache={n_cache} dl={n_dl} err={n_err}  "))
    d <- dalles[i, ]

    if (cache_est_valide(d$id)) {
      n_cache <- n_cache + 1L
      next
    }

    res <- telecharger_dalle(d)
    if (is.null(res)) {
      n_err <- n_err + 1L
    } else {
      n_dl <- n_dl + 1L
    }
    Sys.sleep(1.0)  # rate-limit GéoSAS
  }
  cat(glue("\n  Bilan : {n_cache} cache | {n_dl} telecharges | {n_err} erreurs\n\n"))
}


# =============================================================================
# 5. LECTURE ET ASSEMBLAGE RASTERS NATIONAUX
# =============================================================================

#' Verifie si un fichier CSV cache est corrompu (reponse JSON d'erreur API)
#' Retourne TRUE si le fichier est un vrai CSV SAFRAN, FALSE si corrompu
csv_cache_valide <- function(f) {
  # Lire les premieres lignes brutes pour detecter un JSON {"response":"erreur"...}
  premieres_lignes <- tryCatch(
    readLines(f, n = 2, warn = FALSE, encoding = "UTF-8"),
    error = function(e) ""
  )
  if (length(premieres_lignes) == 0) return(FALSE)   # fichier vide -> dalle hors France OK
  # Si la 1ere ligne commence par '{' c'est du JSON = erreur API = corrompu
  if (grepl("^\\s*\\{", premieres_lignes[1])) return(FALSE)
  TRUE
}

#' Lit tous les CSV du cache et retourne un data.table consolide
lire_cache_complet <- function(dalles = DALLES) {
  cat(">> Lecture du cache SAFRAN...\n")
  fichiers <- chemin_cache(dalles$id)
  fichiers <- fichiers[file.exists(fichiers)]

  if (length(fichiers) == 0) stop("Aucun fichier cache. Lancer recuperer_safran() d'abord.")

  # --- Detection et purge des fichiers corrompus (JSON d'erreur API) ---
  # Cause : requetes v1/v2 utilisaient TMIN_Q (inexistante) -> API retournait JSON
  valides   <- vapply(fichiers, csv_cache_valide, logical(1))
  corrompus <- fichiers[!valides]
  fichiers  <- fichiers[valides]

  if (length(corrompus) > 0) {
    cat(glue("  !! {length(corrompus)} fichier(s) cache corrompus (JSON erreur API) -> suppression\n"))
    file.remove(corrompus)
    cat("  -> Ces dalles seront re-telechargees au prochain appel a recuperer_safran().\n")
    cat("  -> Relancer pipeline_france() pour completer le cache.\n\n")
    if (length(fichiers) == 0) {
      stop("Tous les fichiers cache etaient corrompus. Relancer pipeline_france(force_dl=TRUE).")
    }
  }

  # --- Fichier temoin : premier fichier NON VIDE ---
  # (les dalles hors France ont 0 octets ou juste l'entete -> chercher un vrai)
  temoin <- NULL
  for (f in fichiers) {
    dt_test <- tryCatch(fread(f, nrows = 5), error = function(e) NULL)
    if (!is.null(dt_test) && nrow(dt_test) > 0) {
      temoin <- dt_test
      cat(glue("  Fichier temoin : {basename(f)}\n"))
      break
    }
  }
  if (is.null(temoin)) stop("Aucun fichier cache non-vide trouve. Verifier le telechargement.")

  cat(glue("  Colonnes detectees : {paste(names(temoin), collapse=' | ')}\n"))

  # --- Detection flexible des colonnes (insensible a la casse) ---
  detecter_col <- function(dt, patterns) {
    noms <- tolower(names(dt))
    for (p in patterns) {
      idx <- grep(p, noms, ignore.case = TRUE)
      if (length(idx) > 0) return(names(dt)[idx[1]])
    }
    NULL
  }

  # GéoSAS retourne "time" (pas "date") comme colonne temporelle
  col_date <- detecter_col(temoin, c("^time$", "^date$", "datetime"))
  col_x    <- detecter_col(temoin, c("^x$", "coord_x", "x_l93"))
  col_y    <- detecter_col(temoin, c("^y$", "coord_y", "y_l93"))
  col_tmoy <- detecter_col(temoin, c("t_q", "tmoy", "tmean"))
  # [BUG6] SAFRAN Tmin = TINF_H_Q (TMIN_Q n'existe pas dans l'API)
  col_tmin <- detecter_col(temoin, c("tinf_h_q", "tmin_q", "tmin"))
  col_tmax <- detecter_col(temoin, c("tsup_h_q", "tmax"))
  col_rad  <- detecter_col(temoin, c("ssi_q", "rad", "rayonnement"))

  manquants <- c(
    if (is.null(col_date)) "time/date",
    if (is.null(col_x))    "x (L93)",
    if (is.null(col_y))    "y (L93)",
    if (is.null(col_tmoy)) "T_Q (tmoy)",
    if (is.null(col_tmin)) "TINF_H_Q (tmin)",
    if (is.null(col_tmax)) "TSUP_H_Q (tmax)",
    if (is.null(col_rad))  "SSI_Q (rayonnement)"
  )

  if (length(manquants) > 0) {
    cat(glue("\n  Colonnes disponibles : {paste(names(temoin), collapse=', ')}\n"))
    stop(glue("Colonnes non identifiees : {paste(manquants, collapse=', ')}"))
  }

  cat(glue("  Mapping : time={col_date} | x={col_x} | y={col_y} | ",
           "tmoy={col_tmoy} | tmin={col_tmin} | tmax={col_tmax} | rad={col_rad}\n"))

  # --- Lecture de tous les fichiers valides ---
  liste <- lapply(seq_along(fichiers), function(i) {
    if (i %% 5 == 0) cat(glue("\r  Lecture {i}/{length(fichiers)}  "))
    tryCatch({
      dt <- fread(fichiers[i])
      if (nrow(dt) == 0) return(NULL)
      # Verifier que les colonnes attendues sont presentes dans ce fichier
      cols_presentes <- all(c(col_date, col_x, col_y, col_tmoy,
                               col_tmin, col_tmax, col_rad) %in% names(dt))
      if (!cols_presentes) {
        message("\n  Colonnes manquantes dans ", basename(fichiers[i]), " -> ignore")
        return(NULL)
      }
      setnames(dt,
               c(col_date, col_x, col_y, col_tmoy, col_tmin, col_tmax, col_rad),
               c("date",   "x",   "y",   "t_q",    "tmin_q", "tsup_h_q", "ssi_q"),
               skip_absent = FALSE)
      dt[, .(date, x, y, t_q, tmin_q, tsup_h_q, ssi_q)]
    }, error = function(e) {
      message("\n  Erreur lecture ", basename(fichiers[i]), ": ", e$message)
      NULL
    })
  })
  cat("\n")

  liste <- Filter(Negate(is.null), liste)
  dt    <- rbindlist(liste, use.names = TRUE, fill = TRUE)

  dt[, date := as.Date(date)]
  dt <- dt[date >= DATE_DEB & date <= DATE_FIN_SAFRAN]
  dt <- unique(dt, by = c("date", "x", "y"))

  # Verification du range SSI_Q (Wh/m2 directement, WHm2_FACTOR=1)
  # Valeurs attendues : ~200-800 Wh/m2 hiver, ~1500-5000 Wh/m2 ete
  # Guard : si toutes les valeurs sont NA (ex. debut annee sans donnees) -> skip
  ssi_vals_check <- dt$ssi_q[!is.na(dt$ssi_q)]
  if (length(ssi_vals_check) == 0) {
    warning("SSI_Q : toutes les valeurs sont NA apres filtre date. Verifier DATE_DEB/DATE_FIN_SAFRAN.")
  } else {
    ssi_med <- median(ssi_vals_check)
    cat(glue("  SSI_Q median brut : {round(ssi_med, 0)} J/cm2 (attendu ~200-2500 J/cm2)\n"))
    if (!is.na(ssi_med) && ssi_med > SSI_CHECK_MAX) {
      warning(glue(
        "SSI_Q median ({round(ssi_med,0)}) depasse {SSI_CHECK_MAX} J/cm2.\n",
        "Valeur suspecte. Verifier SSI_Q dans l'API GéoSAS."
      ))
    }
  }

  # Conversion SSI_Q J/cm2 -> Wh/m2 (attendu par barrks phenips-clim)
  # 1 J/cm2 = 10000 J/m2 = 10000/3600 Wh/m2 = 2.7778 Wh/m2
  dt[, ssi_q := ssi_q * WHm2_FACTOR]
  cat(glue("  SSI_Q apres conversion (*{round(WHm2_FACTOR,4)}) : {round(median(dt$ssi_q, na.rm=TRUE), 0)} Wh/m2\n"))

  cat(glue("  {nrow(dt)} lignes | {uniqueN(dt$date)} jours | ",
           "{uniqueN(interaction(dt$x, dt$y))} tuiles SAFRAN\n\n"))
  dt
}

#' Assemble les rasters nationaux depuis le data.table
#' Retourne une liste de SpatRasters avec time() defini (fix barrks)
assembler_rasters_nationaux <- function(dt) {
  cat(">> Assemblage rasters nationaux (France entiere)...\n")
  dates <- sort(unique(dt$date))
  n_j   <- length(dates)

  vars_info <- list(
    list(col = "t_q",      nom = "tmean"),
    list(col = "tmin_q",   nom = "tmin"),
    list(col = "tsup_h_q", nom = "tmax"),
    list(col = "ssi_q",    nom = "rad")
  )

  rasters <- setNames(vector("list", 4), c("tmean", "tmin", "tmax", "rad"))

  for (vi in vars_info) {
    cat(glue("  Variable {vi$nom}...\n"))
    couches <- lapply(seq_along(dates), function(i) {
      if (i %% 20 == 0) cat(glue("\r    Jour {i}/{n_j}  "))
      d   <- dates[i]
      sub <- dt[date == d, .(x, y, v = get(vi$col))]
      sub <- sub[!is.na(v)]

      r <- RST_FRANCE
      if (nrow(sub) < 3) {
        values(r) <- NA_real_
      } else {
        vpts <- vect(as.data.frame(sub[, .(x, y, v)]),
                     geom = c("x", "y"), crs = CRS_L93)
        r    <- rasterize(vpts, RST_FRANCE, field = "v",
                          fun = mean, na.rm = TRUE)
      }
      names(r) <- as.character(d)
      r
    })
    cat("\n")
    rst             <- rast(couches)
    time(rst)       <- dates   # FIX BARRKS : time() obligatoire
    rasters[[vi$nom]] <- rst
  }

  cat(glue("  Rasters assembles : {n_j} jours x {ncell(RST_FRANCE)} cellules\n\n"))
  rasters
}


# =============================================================================
# 6. DUREE DU JOUR (daylength)
# =============================================================================

#' Calcule la duree du jour (Brock 1981) pour chaque cellule et date
#' Necessite de convertir les centroïdes L93 -> WGS84 pour les latitudes
calculer_daylength <- function(rst_ref, dates) {
  cat(">> Calcul daylength (France entiere)...\n")

  # Latitudes WGS84 des centroïdes L93
  lats <- crds(rst_ref) |>
    as.data.frame() |>
    vect(geom = c("x", "y"), crs = CRS_L93) |>
    project(CRS_WGS) |>
    crds() |>
    (\(m) m[, 2])()

  duree <- function(lat_deg, doy) {
    decl   <- -asin(0.39779 * cos(pi/180 * (0.98563*(doy+10) +
                1.914*sin(pi/180 * 0.98563*(doy-2)))))
    cos_ha <- pmin(pmax(-tan(lat_deg*pi/180) * tan(decl), -1), 1)
    2 * acos(cos_ha) * 12 / pi
  }

  cat(glue("  {ncell(rst_ref)} cellules x {length(dates)} jours...\n"))
  couches <- lapply(seq_along(dates), function(i) {
    if (i %% 30 == 0) cat(glue("\r  Jour {i}/{length(dates)}  "))
    r <- rst_ref
    values(r) <- duree(lats, yday(dates[i]))
    names(r)  <- as.character(dates[i])
    r
  })
  cat("\n")

  rst       <- rast(couches)
  time(rst) <- dates  # FIX BARRKS
  rst
}


# =============================================================================
# 7. HELPERS ACCES OBJET BARRKS
# =============================================================================
# [BUG2] L'objet retourne par phenology() peut avoir des structures differentes
# selon la version de barrks (1.x vs futures versions).
# Ces helpers centralisent l'acces et levent une erreur claire si introuvable.

#' Inspecte la structure de l'objet barrks (diagnostic)
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

#' Recupere le raster btmean quel que soit le slot barrks [BUG2]
get_btmean <- function(pheno) {
  # Essayer les slots connus selon les versions de barrks
  if (!is.null(pheno[["data"]][["btmean"]]))    return(pheno$data$btmean)
  if (!is.null(pheno[["storage"]][["btmean"]])) return(pheno$storage$btmean)
  if (!is.null(pheno[["btmean"]]))              return(pheno$btmean)
  # Introuvable : afficher la structure pour diagnostic
  inspect_barrks(pheno)
  stop("btmean introuvable dans l'objet phenips. Voir structure ci-dessus.")
}

#' Recupere le raster de developpement generation 1 [BONUS]
get_gen1 <- function(pheno) {
  if (!is.null(pheno[["development"]][["gen_1"]]))    return(pheno$development$gen_1)
  if (!is.null(pheno[["development"]][["generation_1"]])) return(pheno$development$generation_1)
  if (!is.null(pheno[["gen_1"]]))                     return(pheno$gen_1)
  NULL   # gen1 absente = normal en debut de saison, pas une erreur
}

#' Recupere le raster des generations completes via get_generations_rst()
#' [BUG3] ajoute un message() informatif au lieu de silencieusement retourner NULL
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

#' Lance PHENIPS-Clim via barrks sur la France entiere
#' [BUG1] Ajoute pheno$dates explicitement (absent de l'objet barrks natif)
#' Prerequis : time() defini sur tous les SpatRasters
modeliser_phenips <- function(meteo,
                               exposure = "sunny",
                               scenario = "max") {
  cat(glue(">> Modelisation PHENIPS-Clim [exposure={exposure} | scenario={scenario}]...\n"))
  cat(glue("   {ncell(meteo$tmean)} cellules x {nlyr(meteo$tmean)} jours\n"))

  # Verification que time() est bien defini (prerequis strict de barrks)
  stopifnot(
    "time(tmean) non defini — verifier assembler_rasters_nationaux()" =
      !any(is.na(as.Date(time(meteo$tmean)))),
    "time(tmin) non defini" =
      !any(is.na(as.Date(time(meteo$tmin)))),
    "time(tmax) non defini" =
      !any(is.na(as.Date(time(meteo$tmax)))),
    "time(rad) non defini" =
      !any(is.na(as.Date(time(meteo$rad)))),
    "time(daylength) non defini" =
      !any(is.na(as.Date(time(meteo$daylength))))
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

  # [BUG1] FIX CRITIQUE : ajouter $dates a l'objet barrks
  # phenology() ne stocke pas les dates directement — on les reconstruit
  # depuis le raster d'entree (source de verite)
  pheno$dates <- as.Date(time(meteo$tmean))

  cat(glue("  Modelisation terminee ({length(pheno$dates)} jours).\n\n"))

  # Diagnostic de structure (utile si bug sur un poste)
  if (getOption("phenips.debug", default = FALSE)) inspect_barrks(pheno)

  pheno
}


# =============================================================================
# 9. EXPORT RASTERS
# =============================================================================

exporter_rasters <- function(pheno, dir_out = DIR_OUT) {
  cat(">> Export GeoTIFF (L93, LZW)...\n")

  # [BUG1] pheno$dates maintenant garanti par modeliser_phenips()
  date_j  <- format(max(pheno$dates), "%Y%m%d")
  btmean  <- get_btmean(pheno)   # [BUG2] acces securise

  ecrire <- function(rst, nom) {
    ch <- file.path(dir_out, nom)
    writeRaster(rst, ch, overwrite = TRUE,
                gdal = c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))
    cat(glue("  {nom}\n"))
  }

  # Onset (premier essaimage)
  ecrire(pheno$onset,
         glue("phenips_onset_{ANNEE}_L93.tif"))

  # Temperature sous ecorce (derniere date)
  idx_j <- nlyr(btmean)
  ecrire(btmean[[idx_j]],
         glue("phenips_btmean_{date_j}_L93.tif"))

  # Developpement generation 1 (derniere date)
  gen1 <- get_gen1(pheno)   # [BONUS] helper securise
  if (!is.null(gen1)) {
    idx_g      <- nlyr(gen1)
    gen1_export <- gen1[[idx_g]] * 100
    gen1_export[gen1_export <= 0] <- NA   # masquer absence de developpement
    ecrire(gen1_export,
           glue("phenips_gen1_dev_{date_j}_L93.tif"))
  } else {
    cat("  [INFO] gen1 non disponible (normal avant mi-saison).\n")
  }

  # Generations completes (si disponibles)
  rst_gen <- get_gen_rst(pheno)   # [BUG3] avec message informatif
  if (!is.null(rst_gen) && !all(is.na(values(rst_gen)))) {
    ecrire(rst_gen, glue("phenips_generations_{date_j}_L93.tif"))
  }

  cat(glue("\n  Rasters exportes dans {dir_out}/\n\n"))
}



# =============================================================================
# 10b. EXPORT GEOJSON DES GENERATIONS
# =============================================================================
# Vectorise le raster categoriel des generations en polygones dissous par classe.
# Sortie : GeoJSON WGS84 directement consommable par MapLibre / Leaflet / deck.gl
# Disponible uniquement juillet-octobre (quand au moins une generation est complete)

produire_geojson_generations <- function(pheno, dir_out = DIR_OUT) {
  cat(">> Vectorisation generations -> GeoJSON...\n")

  rst_gen <- get_gen_rst(pheno)
  if (is.null(rst_gen) || all(is.na(values(rst_gen)))) {
    cat("  Pas de generation complete disponible (normal avant juillet).\n\n")
    return(invisible(NULL))
  }

  date_j      <- format(max(pheno$dates), "%Y%m%d")
  date_calcul <- as.character(max(pheno$dates))

  # Reprojection WGS84 (EPSG:4326) pour usage web
  # method = "near" indispensable pour conserver les valeurs entieres/categories
  rst_wgs <- project(rst_gen, "EPSG:4326", method = "near")

  # Vectorisation : pixels contigus de meme valeur -> polygone dissous
  vect_gen <- tryCatch(
    as.polygons(rst_wgs, dissolve = TRUE) |> st_as_sf(),
    error = function(e) {
      message("  Erreur vectorisation : ", e$message)
      NULL
    }
  )
  if (is.null(vect_gen) || nrow(vect_gen) == 0) {
    cat("  Aucun polygone produit.\n\n")
    return(invisible(NULL))
  }

  # Renommer la colonne valeur raster
  names(vect_gen)[1] <- "generation"

  # Labels lisibles pour popup / legende web
  vect_gen$label <- dplyr::case_when(
    vect_gen$generation == 1   ~ "Generation 1 complete",
    vect_gen$generation == 1.5 ~ "Brood soeur (gen 1.5)",
    vect_gen$generation == 2   ~ "Generation 2 complete",
    vect_gen$generation == 2.5 ~ "Brood soeur gen 2",
    vect_gen$generation == 3   ~ "Generation 3 complete",
    TRUE ~ paste0("Gen ", vect_gen$generation)
  )

  # Couleurs MapLibre par generation (coherentes avec barrks_colors())
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

  # Simplification geometrique legere pour alleger le GeoJSON web
  # tolerance en degres (~500m a la latitude France)
  vect_gen <- sf::st_simplify(vect_gen, dTolerance = 0.005, preserveTopology = TRUE)

  # Supprimer les geometries vides issues de la simplification
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
# Vectorise la date de premier essaimage par pixel en polygones dissous par date.
# Le raster onset (nlyr = nb jours) contient 1 par jour d'essaimage, 0 sinon.
# On extrait le premier jour > 0 par pixel, on le convertit en date calendaire
# jj/mm/aaaa, puis on vectorise par date unique.
# Sortie : GeoJSON WGS84 avec champ date_onset "jj/mm/aaaa" pour popup web.

produire_geojson_onset <- function(pheno, dir_out = DIR_OUT) {
  cat(">> Vectorisation onset (premier essaimage) -> GeoJSON...\n")

  onset  <- pheno$onset
  dates  <- pheno$dates   # vecteur Date de longueur nlyr(onset)

  if (is.null(onset) || nlyr(onset) == 0) {
    cat("  Raster onset absent.\n\n")
    return(invisible(NULL))
  }

  date_j <- format(max(dates), "%Y%m%d")

  # --- Extraire l index du premier jour d essaimage par pixel ---
  # max.col sur la matrice (nlyr colonnes) avec ties.method="first"
  # Uniquement sur les pixels ayant au moins un 1
  onset_vals <- values(onset)   # matrice ncell x nlyr
  has_onset  <- rowSums(onset_vals > 0, na.rm = TRUE) > 0

  idx_onset <- rep(NA_integer_, nrow(onset_vals))
  if (any(has_onset, na.rm = TRUE)) {
    sub_mat <- onset_vals[has_onset, , drop = FALSE]
    idx_onset[has_onset] <- max.col(sub_mat > 0, ties.method = "first")
  }

  # --- Convertir l index en date calendaire ---
  # idx_onset est un entier 1..nlyr -> dates[idx_onset]
  date_onset_vec <- rep(NA_character_, length(idx_onset))
  valides <- !is.na(idx_onset)
  date_onset_vec[valides] <- format(dates[idx_onset[valides]], "%d/%m/%Y")

  # --- Construire un raster avec la date encodee en jours depuis origine ---
  # On stocke le numero de jour dans l annee (DOY) pour pouvoir vectoriser
  # puis on reattache la date lisible en attribut apres vectorisation
  doy_onset_vec <- rep(NA_real_, length(idx_onset))
  doy_onset_vec[valides] <- as.integer(format(dates[idx_onset[valides]], "%j"))

  rst_onset_doy <- RST_FRANCE
  values(rst_onset_doy) <- doy_onset_vec

  # Pixels sans essaimage = NA (zones froides, hors saison)
  rst_onset_doy[is.na(rst_onset_doy)] <- NA

  # --- Reprojection WGS84 ---
  rst_wgs <- project(rst_onset_doy, "EPSG:4326", method = "near")

  # --- Vectorisation par DOY unique ---
  vect_onset <- tryCatch(
    as.polygons(rst_wgs, dissolve = TRUE) |> st_as_sf(),
    error = function(e) {
      message("  Erreur vectorisation onset : ", e$message)
      NULL
    }
  )
  if (is.null(vect_onset) || nrow(vect_onset) == 0) {
    cat("  Aucun polygone onset produit (pas encore d essaimage).\n\n")
    return(invisible(NULL))
  }

  names(vect_onset)[1] <- "doy"

  # Reattacher la date calendaire lisible depuis le DOY
  origine <- as.Date(glue("{ANNEE}-01-01"))
  vect_onset$date_onset   <- format(origine + vect_onset$doy - 1L, "%d/%m/%Y")
  vect_onset$date_iso     <- format(origine + vect_onset$doy - 1L, "%Y-%m-%d")
  vect_onset$semaine      <- as.integer(format(origine + vect_onset$doy - 1L, "%V"))
  vect_onset$date_calcul  <- as.character(max(dates))
  vect_onset$annee        <- ANNEE

  # Supprimer la colonne DOY brute (non necessaire cote web)
  vect_onset$doy <- NULL

  # --- Simplification geometrique ---
  vect_onset <- sf::st_simplify(vect_onset, dTolerance = 0.005,
                                 preserveTopology = TRUE)
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
# 10. VISUALISATIONS NATIONALES
# =============================================================================

#' Carte nationale generique (gradient ou categorielle)
carte_france <- function(rst_l93, titre, legende, palette,
                          breaks     = NULL,
                          type       = c("gradient", "categoriel"),
                          date_label = NULL,
                          fichier    = NULL) {
  type <- match.arg(type)
  if (is.null(date_label)) date_label <- format(DATE_FIN_SAFRAN, "%d/%m/%Y")
  rst_wgs <- project(rst_l93, CRS_WGS)

  p <- ggplot() +
    geom_spatraster(data = rst_wgs) +
    {
      if (type == "gradient") {
        scale_fill_gradientn(
          colours  = palette,
          values   = if (!is.null(breaks)) rescale(breaks) else NULL,
          na.value = "lightblue",
          name     = legende,
          guide    = guide_colorbar(barwidth = 14, barheight = 0.8,
                                    title.position = "top")
        )
      } else {
        scale_fill_manual(
          values   = palette,
          na.value = "lightblue",
          name     = legende,
          drop     = FALSE
        )
      }
    } +
    labs(
      title    = titre,
      subtitle = glue("PHENIPS-Clim | SAFRAN INRAE | {date_label} | {ANNEE}"),
      caption  = "GéoSAS INRAE — CNPF/CRPF AuRA",
      x = NULL, y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold"),
      plot.caption    = element_text(size = 8, color = "grey50"),
      legend.position = "bottom"
    )

  if (!is.null(fichier)) {
    ggsave(fichier, plot = p, width = 10, height = 9, dpi = 150)
    cat(glue("  Carte : {fichier}\n"))
  }
  invisible(p)
}


#' Genere les cartes nationales de synthese
generer_cartes <- function(pheno, dir_out = DIR_OUT) {
  cat(">> Generation des cartes nationales...\n")

  # [BUG1] pheno$dates garanti
  date_j     <- format(max(pheno$dates), "%Y%m%d")
  date_label <- format(max(pheno$dates), "%d/%m/%Y")
  btmean     <- get_btmean(pheno)   # [BUG2]

  # 1. Carte onset (premier essaimage)
  # [BUG4] Remplacement de apply() ligne par ligne par max.col() vectorise
  # apply() sur ~14000 lignes x n jours = tres lent ; max.col() instantane
  onset_vals <- values(pheno$onset)
  has_onset  <- rowSums(onset_vals > 0, na.rm = TRUE) > 0

  date_onset_idx <- rep(NA_integer_, nrow(onset_vals))
  if (any(has_onset, na.rm = TRUE)) {
    # max.col avec ties.method="first" donne le premier TRUE par ligne
    # on applique seulement sur les lignes qui ont au moins un 1
    sub_mat <- onset_vals[has_onset, , drop = FALSE]
    date_onset_idx[has_onset] <- max.col(sub_mat > 0, ties.method = "first")
  }

  rst_onset_date <- RST_FRANCE
  values(rst_onset_date) <- as.integer(date_onset_idx)

  carte_france(
    rst_onset_date,
    titre      = "Date de premier essaimage — Ips typographus — France",
    legende    = glue("Jour de l'annee (1=01/01/{ANNEE})"),
    palette    = c("grey95", "#FFF3B0", "#FDE725FF", "#35B779FF",
                   "#31688EFF", "#440154FF"),
    type       = "gradient",
    date_label = date_label,
    fichier    = file.path(dir_out, glue("phenips_onset_{date_j}.png"))
  )

  # 2. Temperature sous ecorce (derniere date)
  idx_j  <- nlyr(btmean)
  rst_bt <- btmean[[idx_j]]
  carte_france(
    rst_bt,
    titre      = glue("Temperature sous ecorce — Ips typographus — {date_label}"),
    legende    = "Temp. sous ecorce (C)",
    palette    = c("grey95", "#FFF3B0", "#FDE725FF", "#35B779FF",
                   "#F39C12", "#C0392B"),
    breaks     = c(0, 8, 16.5, 20, 25, 35),
    type       = "gradient",
    date_label = date_label,
    fichier    = file.path(dir_out, glue("phenips_btmean_{date_j}.png"))
  )

  # 3. Developpement generation 1 (si disponible)
  gen1 <- get_gen1(pheno)   # [BONUS]
  if (!is.null(gen1)) {
    idx_g  <- nlyr(gen1)
    rst_g1 <- gen1[[idx_g]] * 100  # -> pourcentage
    rst_g1[rst_g1 <= 0] <- NA       # masquer absence de developpement (montagnes, hiver)
    carte_france(
      rst_g1,
      titre      = glue("Developpement generation 1 — Ips typographus — {date_label}"),
      legende    = "Developpement gen1 (%)",
      palette    = c("grey95", "#FFF7BC", "#FDE725FF", "#35B779FF", "#31688EFF", "#440154FF"),
      breaks     = c(0, 10, 25, 50, 75, 100),
      type       = "gradient",
      date_label = date_label,
      fichier    = file.path(dir_out, glue("phenips_gen1_{date_j}.png"))
    )
  }

  # 4. Generations completes (si disponibles en ete/automne)
  rst_gen  <- get_gen_rst(pheno)   # [BUG3]
  gen_vide <- is.null(rst_gen) || all(is.na(values(rst_gen)))

  if (!gen_vide) {
    carte_france(
      rst_gen,
      titre      = glue("Generations Ips typographus — France — {date_label}"),
      legende    = "Generation",
      palette    = barrks_colors(),
      type       = "categoriel",
      date_label = date_label,
      fichier    = file.path(dir_out, glue("phenips_generations_{date_j}.png"))
    )
    cat("  Carte generations produite.\n")
  } else {
    cat("  Carte generations : pas encore de generation complete (normal avant ete).\n")
  }

  cat(glue("  Cartes exportees dans {dir_out}/\n\n"))
}


# =============================================================================
# 11. SYNTHESE CONSOLE
# =============================================================================

synthese_console <- function(pheno) {
  cat("\n", strrep("=", 60), "\n")
  cat("  SYNTHESE PHENIPS-Clim — France metropolitaine\n")
  cat(strrep("=", 60), "\n")

  # [BUG1] pheno$dates garanti
  date_j  <- max(pheno$dates)
  btmean  <- get_btmean(pheno)   # [BUG2]
  idx_j   <- nlyr(btmean)

  bt_vals  <- values(btmean[[idx_j]])
  # [BONUS] n_valide depuis btmean (pas RST_FRANCE qui inclut mer/hors-France)
  n_valide <- sum(!is.na(bt_vals))
  n_total  <- ncell(RST_FRANCE)

  cat(glue("  Date      : {format(date_j, '%d/%m/%Y')}\n"))
  cat(glue("  Cellules  : {n_valide} valides / {n_total} total\n\n"))

  # Onset
  onset_vals  <- values(pheno$onset)
  n_essaime   <- sum(rowSums(onset_vals > 0, na.rm = TRUE) > 0, na.rm = TRUE)
  pct_essaime <- round(100 * n_essaime / n_valide, 1)
  cat(glue("  Essaimage : {n_essaime} cellules ({pct_essaime}%) ont essaime\n\n"))

  # Temperature sous ecorce
  cat("  --- Temperature sous ecorce (btmean) ---\n")
  cat(glue("  Moyenne : {round(mean(bt_vals, na.rm=TRUE), 1)} C\n"))
  cat(glue("  Min     : {round(min(bt_vals,  na.rm=TRUE), 1)} C\n"))
  cat(glue("  Max     : {round(max(bt_vals,  na.rm=TRUE), 1)} C\n"))
  n_au_dessus <- sum(bt_vals >= 16.5, na.rm = TRUE)
  cat(glue("  >= seuil essaimage (16.5C) : {n_au_dessus} cellules",
           " ({round(100*n_au_dessus/n_valide,1)}%)\n\n"))

  # Generation 1
  gen1 <- get_gen1(pheno)   # [BONUS]
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


# =============================================================================
# 12. PIPELINE PRINCIPAL
# =============================================================================

#' Pipeline complet France entiere — mise a jour mensuelle
#' Etapes : cache SAFRAN -> assemblage -> daylength -> PHENIPS -> export -> cartes
#'
#' @param exposure  "sunny" (defaut) ou "shaded"
#' @param scenario  "max" (defaut) ou "mean"
#' @param force_dl  TRUE = re-telecharge meme si cache valide
#' @return Objet pheno barrks enrichi (avec $dates)
pipeline_france <- function(exposure = "sunny",
                             scenario = "max",
                             force_dl = FALSE) {
  cat("\n", strrep("=", 60), "\n")
  cat("  PHENIPS-Clim — France metropolitaine — SAFRAN v2\n")
  cat(strrep("=", 60), "\n\n")
  t0 <- proc.time()

  # IMPORTANT v3 : si tu viens d'une v1/v2, le cache est corrompu (JSON erreur API)
  # car les requetes utilisaient TMIN_Q (inexistante). Purger avant le premier run :
  #   file.remove(list.files(DIR_CACHE, pattern="\\.csv$", full.names=TRUE))
  # OU utiliser force_dl = TRUE au premier lancement

  # 1. Telechargement / mise a jour cache SAFRAN
  if (force_dl) {
    fichiers_cache <- list.files(DIR_CACHE, pattern = "\\.csv$", full.names = TRUE)
    if (length(fichiers_cache) > 0) {
      file.remove(fichiers_cache)
      cat(glue("  {length(fichiers_cache)} fichiers cache supprimes (force_dl=TRUE)\n"))
    }
  }
  recuperer_safran()

  # 2. Lecture et assemblage
  dt    <- lire_cache_complet()
  meteo <- assembler_rasters_nationaux(dt)
  rm(dt); gc()

  # 3. Daylength
  dates_serie     <- as.Date(time(meteo$tmean))
  meteo$daylength <- calculer_daylength(RST_FRANCE, dates_serie)

  # 4. Modelisation PHENIPS-Clim
  # [BUG1] pheno$dates sera ajoute dans modeliser_phenips()
  pheno <- modeliser_phenips(meteo, exposure = exposure, scenario = scenario)
  rm(meteo); gc()

  # 5. Export rasters GeoTIFF
  exporter_rasters(pheno)

  # 5b. Export GeoJSON generations (vectoriel web)
  produire_geojson_generations(pheno)

  # 5c. Export GeoJSON onset (date premier essaimage)
  produire_geojson_onset(pheno)

  # 6. Cartes PNG
  generer_cartes(pheno)

  # 7. Synthese console
  synthese_console(pheno)

  duree <- (proc.time() - t0)[["elapsed"]]
  cat(glue("  Duree totale : {round(duree/60, 1)} min\n"))
  cat(glue("  Sorties      : {DIR_OUT}/\n\n"))

  invisible(pheno)
}


# =============================================================================
# 13. LANCEMENT
# =============================================================================

# --- Test connexion API + verification variables AVANT tout pipeline ---
# A lancer systematiquement sur un nouveau poste ou apres changement d'API
# Utilise une dalle centree sur la France (Massif Central) qui contient des donnees
# La dalle 1 (coin SO) est hors France -> 0 lignes = normal, inutile pour le test
#if (FALSE) {
  # Dalle Massif Central (~centre France, garantie non vide)
# dalle_test <- DALLES[DALLES$xmin == 600000 & DALLES$ymin == 6200000, ]
# if (nrow(dalle_test) == 0) {
    # Fallback : prendre la dalle la plus centrale disponible
#   centre_x <- (FRANCE_XMIN + FRANCE_XMAX) / 2
#   centre_y <- (FRANCE_YMIN + FRANCE_YMAX) / 2
#   dalle_test <- DALLES[
#     which.min(abs(DALLES$xmin - centre_x) + abs(DALLES$ymin - centre_y)), ]
# }
# cat(glue("Test dalle Massif Central : bbox={dalle_test$xmin},{dalle_test$ymin},",
#          "{dalle_test$xmax},{dalle_test$ymax}\n"))

# resp_test <- request(GEOSAS_URL) |>
#   req_url_query(
#     bbox             = glue("{dalle_test$xmin},{dalle_test$ymin},",
#                             "{dalle_test$xmax},{dalle_test$ymax}"),
#     crs              = GEOSAS_CRS,
#     `parameter-name` = GEOSAS_VARS,
#     f                = "CSV",
#     datetime         = glue("{DATE_DEB}/{DATE_DEB + 7}")   # 7 jours seulement
#   ) |>
#   req_timeout(60) |>
#   req_perform()
#
# cat("Status HTTP :", resp_status(resp_test), "\n")
# csv_test <- fread(text = resp_body_string(resp_test))
# cat("Colonnes    :", paste(names(csv_test), collapse = " | "), "\n")
# cat("Lignes      :", nrow(csv_test), "\n")
#
# # Si 0 lignes = dalle hors France ou periode sans donnees
# if (nrow(csv_test) == 0) {
#   cat("!! Aucune donnee retournee pour cette dalle/periode.\n")
#   cat("   Verifier que DATE_DEB est dans la plage SAFRAN disponible.\n")
# } else {
#   print(head(csv_test, 3))
#
    # Verification TINF_H_Q presente [BUG6]
#   col_tmin_test <- names(csv_test)[grep("tinf_h_q", tolower(names(csv_test)))]
#   if (length(col_tmin_test) > 0) {
#     cat("TINF_H_Q OK — Tmin presente\n")
#   } else {
#     cat("!! ATTENTION : TINF_H_Q non trouvee dans la reponse !\n")
#     cat("   Colonnes disponibles :", paste(names(csv_test), collapse = ", "), "\n")
#   }

    # Verification SSI_Q en J/cm2 [BUG7]
    # Guard NA : range/median sur 0 valeurs retourne Inf/-Inf et NA
#   col_rad <- names(csv_test)[grep("ssi", tolower(names(csv_test)))]
#   if (length(col_rad) > 0) {
#     ssi_vals  <- csv_test[[col_rad]]
#     ssi_vals  <- ssi_vals[!is.na(ssi_vals)]
#     if (length(ssi_vals) == 0) {
#       cat("!! SSI_Q presente mais toutes les valeurs sont NA.\n")
#     } else {
#       ssi_range <- range(ssi_vals)
#       ssi_med   <- median(ssi_vals)
#       cat(glue("SSI_Q : range={round(ssi_range[1],2)}-{round(ssi_range[2],2)}",
#                " | mediane={round(ssi_med,2)} J/cm2\n"))
#       if (is.na(ssi_med)) {
#         cat("!! SSI_Q mediane NA — impossible de valider l'unite.\n")
#       } else if (ssi_med > SSI_CHECK_MAX) {
#         cat(glue("!! ATTENTION : mediane ({round(ssi_med,1)}) > {SSI_CHECK_MAX} J/cm2",
#                  " = unite suspecte, verifier WHm2_FACTOR\n"))
#       } else {
#         cat(glue("OK : SSI_Q en J/cm2 (attendu ~0-25) ->",
#                  " WHm2_FACTOR={round(WHm2_FACTOR,4)} correct\n"))
#       }
#     }
#   } else {
#     cat("!! ATTENTION : SSI_Q non trouvee dans la reponse CSV !\n")
#   }
# }
#}

# --- Test sur une dalle unique (Massif Central) sans pipeline France entiere ---
#if (FALSE) {
# dalle_mc <- DALLES[DALLES$xmin == 600000 & DALLES$ymin == 6200000, ]
# if (nrow(dalle_mc) == 0) dalle_mc <- DALLES[8, ]
#
# cat("Test sur dalle Massif Central...\n")
# telecharger_dalle(dalle_mc,
#                   date_deb = as.Date(glue("{ANNEE}-01-01")),
#                   date_fin = as.Date(glue("{ANNEE}-03-31")))

# dt_test    <- lire_cache_complet(dalles = dalle_mc)
# meteo_test <- assembler_rasters_nationaux(dt_test)
#
# dates_t             <- as.Date(time(meteo_test$tmean))
# meteo_test$daylength <- calculer_daylength(RST_FRANCE, dates_t)
#
  # Activer le diagnostic de structure barrks
# options(phenips.debug = TRUE)
# pheno_test <- modeliser_phenips(meteo_test)
# options(phenips.debug = FALSE)

# synthese_console(pheno_test)
#}

# --- Pipeline complet France entiere ---
# Mise a jour mensuelle — lancer le 6 ou 7 du mois
# Duree estimee : 15-30 min selon la connexion et le cache
#if (TRUE) {
# pheno <- pipeline_france(exposure = "sunny", scenario = "max")
#}

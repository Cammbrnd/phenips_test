# =============================================================================
# SAFRAN COMMON — fonctions partagees par PHENIPS, CHAPY, STENO, Monochamus
# -----------------------------------------------------------------------------
# Regroupe tout ce qui etait duplique a l'identique dans les 4 scripts :
# grille de dalles, telechargement, cache CSV, lecture, assemblage rasters,
# daylength. Chaque modele n'a plus qu'a sourcer ce fichier puis definir SES
# propres parametres biologiques et sa fonction modeliser_xxx().
#
# IMPORTANT : on telecharge systematiquement le SUPERSET des 4 variables
# (T_Q, TINF_H_Q, TSUP_H_Q, SSI_Q), meme pour les modeles qui n'utilisent pas
# le rayonnement (Monochamus). Cela permet un cache 100% partage entre les
# 4 pipelines : un seul jeu de requetes API GeoSAS pour tous les modeles.
# =============================================================================

library(terra)
library(httr2)
library(glue)
library(data.table)
library(lubridate)

# --- Periode (identique pour tous les modeles) ---
ANNEE    <- as.integer(format(Sys.Date(), "%Y"))
DATE_DEB <- as.Date(glue("{ANNEE}-01-01"))
DATE_FIN_SAFRAN <- Sys.Date() - 5
if (DATE_FIN_SAFRAN < DATE_DEB) DATE_FIN_SAFRAN <- DATE_DEB

# --- API GeoSAS ---
GEOSAS_URL  <- "https://api.geosas.fr/edr/collections/safran-isba/cube"
GEOSAS_CRS  <- "EPSG:2154"
GEOSAS_VARS <- "T_Q,TINF_H_Q,TSUP_H_Q,SSI_Q"   # superset partage par les 4 modeles

# --- Grille SAFRAN France metropolitaine (L93) ---
FRANCE_XMIN <- 100000
FRANCE_XMAX <- 1200000
FRANCE_YMIN <- 6050000
FRANCE_YMAX <- 7150000
SAFRAN_RES  <- 8000
DALLE_KM    <- 200
DALLE_M     <- DALLE_KM * 1000
CRS_L93 <- "EPSG:2154"
CRS_WGS <- "EPSG:4326"

# --- Repertoire de cache PARTAGE (un seul, pour les 4 modeles) ---
DIR_CACHE <- "data/cache_safran"
dir.create(DIR_CACHE, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# GRILLE DE DALLES
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

creer_raster_france <- function() {
  rast(
    xmin = FRANCE_XMIN, xmax = FRANCE_XMAX,
    ymin = FRANCE_YMIN, ymax = FRANCE_YMAX,
    resolution = SAFRAN_RES, crs = CRS_L93
  )
}

RST_FRANCE <- creer_raster_france()


# =============================================================================
# TELECHARGEMENT + CACHE (identique aux 4 scripts d'origine)
# =============================================================================

chemin_cache <- function(dalle_id) {
  file.path(DIR_CACHE, glue("safran_dalle{sprintf('%03d', dalle_id)}_{ANNEE}.csv"))
}

csv_cache_valide <- function(f) {
  premieres_lignes <- tryCatch(readLines(f, n = 2, warn = FALSE, encoding = "UTF-8"),
                                error = function(e) "")
  if (length(premieres_lignes) == 0) return(FALSE)
  if (grepl("^\\s*\\{", premieres_lignes[1])) return(FALSE)
  TRUE
}

cache_est_valide <- function(dalle_id) {
  f <- chemin_cache(dalle_id)
  if (!file.exists(f)) return(FALSE)
  tryCatch({
    entete <- names(fread(f, nrows = 0))
    if (length(entete) == 0) return(TRUE)
    col_t <- entete[grep("^time$|^date$|datetime", tolower(entete))][1]
    if (is.na(col_t)) return(TRUE)
    df <- fread(f, select = col_t)
    if (nrow(df) == 0) return(TRUE)
    max_date <- as.Date(max(df[[col_t]], na.rm = TRUE))
    if (is.na(max_date) || is.infinite(as.numeric(max_date))) return(TRUE)
    max_date >= DATE_FIN_SAFRAN
  }, error = function(e) FALSE)
}

telecharger_dalle <- function(dalle, vars = GEOSAS_VARS,
                              date_deb = DATE_DEB, date_fin = DATE_FIN_SAFRAN) {
  bbox    <- glue("{dalle$xmin},{dalle$ymin},{dalle$xmax},{dalle$ymax}")
  periode <- glue("{date_deb}/{date_fin}")

  req <- request(GEOSAS_URL) |>
    req_url_query(bbox = bbox, crs = GEOSAS_CRS, `parameter-name` = vars,
                  f = "CSV", datetime = periode) |>
    req_timeout(180) |>
    req_retry(max_tries = 5, backoff = \(i) 30 * i,
              is_transient = \(resp) resp_status(resp) %in% c(429, 503, 504))

  resp <- tryCatch(req_perform(req), error = function(e) {
    message("  ERREUR dalle ", dalle$id, ": ", e$message)
    NULL
  })
  if (is.null(resp)) return(NULL)
  if (resp_status(resp) != 200) {
    message("  HTTP ", resp_status(resp), " dalle ", dalle$id)
    return(NULL)
  }
  writeBin(resp_body_raw(resp), chemin_cache(dalle$id))
  invisible(chemin_cache(dalle$id))
}

recuperer_safran <- function(dalles = DALLES) {
  cat(glue(">> Telechargement SAFRAN partage ({DATE_DEB} -> {DATE_FIN_SAFRAN})...\n"))
  n <- nrow(dalles)
  n_cache <- 0L; n_dl <- 0L; n_err <- 0L
  for (i in seq_len(n)) {
    cat(glue("\r  Dalle {i}/{n} | cache={n_cache} dl={n_dl} err={n_err}  "))
    d <- dalles[i, ]
    if (cache_est_valide(d$id)) { n_cache <- n_cache + 1L; next }
    res <- telecharger_dalle(d)
    if (is.null(res)) n_err <- n_err + 1L else n_dl <- n_dl + 1L
    Sys.sleep(1.0)
  }
  cat(glue("\n  Bilan : {n_cache} cache | {n_dl} telecharges | {n_err} erreurs\n\n"))
}


# =============================================================================
# LECTURE + ASSEMBLAGE (superset des 4 variables, chaque modele pioche ce qu'il utilise)
# =============================================================================

lire_cache_complet <- function(dalles = DALLES) {
  cat(">> Lecture du cache SAFRAN partage...\n")
  fichiers <- chemin_cache(dalles$id)
  fichiers <- fichiers[file.exists(fichiers)]
  if (length(fichiers) == 0) stop("Aucun fichier cache. Lancer recuperer_safran() d'abord.")

  valides   <- vapply(fichiers, csv_cache_valide, logical(1))
  corrompus <- fichiers[!valides]
  fichiers  <- fichiers[valides]
  if (length(corrompus) > 0) {
    cat(glue("  !! {length(corrompus)} fichier(s) corrompus -> suppression\n"))
    file.remove(corrompus)
    if (length(fichiers) == 0) stop("Tous les fichiers cache etaient corrompus.")
  }

  temoin <- NULL
  for (f in fichiers) {
    dt_test <- tryCatch(fread(f, nrows = 5), error = function(e) NULL)
    if (!is.null(dt_test) && nrow(dt_test) > 0) { temoin <- dt_test; break }
  }
  if (is.null(temoin)) stop("Aucun fichier cache non-vide trouve.")

  detecter_col <- function(dt, patterns) {
    noms <- tolower(names(dt))
    for (p in patterns) {
      idx <- grep(p, noms, ignore.case = TRUE)
      if (length(idx) > 0) return(names(dt)[idx[1]])
    }
    NULL
  }

  col_date <- detecter_col(temoin, c("^time$", "^date$", "datetime"))
  col_x    <- detecter_col(temoin, c("^x$", "coord_x", "x_l93"))
  col_y    <- detecter_col(temoin, c("^y$", "coord_y", "y_l93"))
  col_tmoy <- detecter_col(temoin, c("t_q", "tmoy", "tmean"))
  col_tmin <- detecter_col(temoin, c("tinf_h_q", "tmin_q", "tmin"))
  col_tmax <- detecter_col(temoin, c("tsup_h_q", "tmax"))
  col_rad  <- detecter_col(temoin, c("ssi_q", "rad", "rayonnement"))

  liste <- lapply(seq_along(fichiers), function(i) {
    if (i %% 5 == 0) cat(glue("\r  Lecture {i}/{length(fichiers)}  "))
    tryCatch({
      dt <- fread(fichiers[i])
      if (nrow(dt) == 0) return(NULL)
      requis <- c(col_date, col_x, col_y, col_tmoy, col_tmin, col_tmax, col_rad)
      if (!all(requis %in% names(dt))) return(NULL)
      setnames(dt, requis, c("date","x","y","t_q","tmin_q","tsup_h_q","ssi_q"),
               skip_absent = FALSE)
      dt[, .(date, x, y, t_q, tmin_q, tsup_h_q, ssi_q)]
    }, error = function(e) NULL)
  })
  cat("\n")
  liste <- Filter(Negate(is.null), liste)
  dt <- rbindlist(liste, use.names = TRUE, fill = TRUE)
  dt[, date := as.Date(date)]
  dt <- dt[date >= DATE_DEB & date <= DATE_FIN_SAFRAN]
  dt <- unique(dt, by = c("date", "x", "y"))

  # Conversion rayonnement J/cm2 -> Wh/m2 (voir diagnostiquer_rayonnement() dans
  # steno_france_safran_v1.R si tu veux revalider cette hypothese d'unite)
  dt[, ssi_q := ssi_q * (10000 / 3600)]

  cat(glue("  {nrow(dt)} lignes | {uniqueN(dt$date)} jours | ",
           "{uniqueN(interaction(dt$x, dt$y))} tuiles SAFRAN\n\n"))
  dt
}

assembler_rasters_nationaux <- function(dt) {
  cat(">> Assemblage rasters nationaux (partage)...\n")
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
      d   <- dates[i]
      sub <- dt[date == d, .(x, y, v = get(vi$col))]
      sub <- sub[!is.na(v)]
      r <- RST_FRANCE
      if (nrow(sub) < 3) {
        values(r) <- NA_real_
      } else {
        vpts <- vect(as.data.frame(sub[, .(x, y, v)]), geom = c("x","y"), crs = CRS_L93)
        r <- rasterize(vpts, RST_FRANCE, field = "v", fun = mean, na.rm = TRUE)
      }
      names(r) <- as.character(d)
      r
    })
    rst <- rast(couches)
    time(rst) <- dates
    rasters[[vi$nom]] <- rst
  }
  cat(glue("  Rasters assembles : {n_j} jours x {ncell(RST_FRANCE)} cellules\n\n"))
  rasters
}

calculer_daylength <- function(rst_ref, dates) {
  cat(">> Calcul daylength (partage)...\n")
  lats <- crds(rst_ref) |> as.data.frame() |>
    vect(geom = c("x","y"), crs = CRS_L93) |> project(CRS_WGS) |>
    crds() |> (\(m) m[, 2])()
  duree <- function(lat_deg, doy) {
    decl <- -asin(0.39779 * cos(pi/180 * (0.98563*(doy+10) + 1.914*sin(pi/180 * 0.98563*(doy-2)))))
    cos_ha <- pmin(pmax(-tan(lat_deg*pi/180) * tan(decl), -1), 1)
    2 * acos(cos_ha) * 12 / pi
  }
  couches <- lapply(seq_along(dates), function(i) {
    r <- rst_ref
    values(r) <- duree(lats, yday(dates[i]))
    names(r) <- as.character(dates[i])
    r
  })
  rst <- rast(couches)
  time(rst) <- dates
  rst
}
